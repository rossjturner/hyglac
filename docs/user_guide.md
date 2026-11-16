# Hyglac user guide

Hyglac is a computational framework for simulating water flow through ice channels and predicting the resulting seismic signals. The method is described in Turner et al. (2026), Hyglac: a computational framework to inform interpretation of hydrological processes in glacier environments through simulated seismic signals. This guide provides instructions for using the software, descriptions of its functions and parameters, and guidance on interpreting the results.

The guide describes [`hyglac.py`](../hyglac.py) as provided in the accompanying software release.

**Contents**

1. [Overview](#1-overview)
2. [Installation](#2-installation)
3. [Conventions](#3-conventions)
4. [Walkthrough: the verification test](#4-walkthrough-the-verification-test)
5. [API reference](#5-api-reference)
6. [Choosing resolution and parameters](#6-choosing-resolution-and-parameters)
7. [Using output from another SPH simulation](#7-using-output-from-another-sph-simulation)
8. [Limitations](#8-limitations)

## 1. Overview

Hyglac has three modular stages. Each stage can be used independently, provided the necessary inputs are available:

| Stage | Function | Reads | Produces | Paper |
| --- | --- | --- | --- | --- |
| Hydrodynamic simulation | `run_seismic_sph` | Domain dimensions and geometry functions | SPH snapshots (HDF5) | §2 |
| Seismic force vectors | `seismic_sph_interactions` | SPH snapshots (HDF5) | Interface forces (HDF5) | §3.1–3.2 |
| | `seismic_sph_forces` or `seismic_sph_forces_defaults` | Interface forces (HDF5) | Force vectors (DataFrame) | §3.3 |
| Seismic wave propagation | `earth_model`, `seismic_sph_signals` | Force vectors (DataFrame or CSV) and Earth model | Seismograms (DataFrame) | §4 |

The hydrodynamic simulation is by far the most expensive stage. Once it has run, the forces can be recomputed with different grids, and seismograms can be calculated for any number of seismometer positions and Earth models, without rerunning the hydrodynamic simulation.

## 2. Installation

The repository includes an [`environment.yml`](../environment.yml) file recording the Conda environment used to develop and run Hyglac, including Python 3.9, PySPH 1.0a6, pyprop8 1.1.5, NumPy, SciPy, pandas, h5py and Matplotlib.

```bash
conda env create -f environment.yml
conda activate hyglac
```

Hyglac is a single module rather than an installable package. Run Python or Jupyter from the repository root and import it as:

```python
import hyglac as ss
```

Some practical points:

- PySPH and pyprop8 are both imported when the module loads, so both must be installed even if you only use one stage.

## 3. Conventions

### Coordinate frames

Three coordinate frames are used throughout Hyglac: the hydrodynamic simulation frame, the seismic source frame, and the seismometer frame.

**Simulation frame** (stages 1 and 2). `x` runs along the tank from `0` to its length. `y` runs across the tank, centred on `0`, from `−width/2` to `+width/2`. `z` points upward from the tank floor at `0`. The geometry functions, the `xrange`/`yrange`/`zrange` windows and the coordinates in the interactions file all use this frame.

**Source frame** (output of `seismic_sph_forces`). The horizontal origin is moved to `(x0, y0)`, so `loc_1 = x − x0` and `loc_2 = y − y0`, rounded to the force grid. This corresponds to the paper's coordinate system, whose origin is the centre of the channel bends. The vertical coordinate `loc_3` is negative depth: `−loc_3` is the depth of the source below the surface of the Earth model. In practice `z0` is the depth of the simulation's floor: 500 m in the paper, at the base of the ice layer.

**Seismometer position** (input to `seismic_sph_signals`). The seismometer's `x` and `y` are measured in the source frame, i.e. relative to `(x0, y0)`. In this case, its `z` is a positive depth consistent with `pyprop8`. The paper writes the same seismometer position as `z = −0.5 m`; in the code it is `z=0.5`.

### Force sign

The unit normal to the interface points out of the ice and into the water (paper §3.1). The interactions file and the force DataFrame store `pressure × area × normal`, which points into the water. The force the water exerts on the ice is the negative of this (paper Eq. 8). The `seismic_sph_signals` function applies that minus sign internally, so pass the forces in exactly as Hyglac produces them.

### Units

All inputs and outputs are expressed in SI units (m, s, N, Pa, kg m⁻³, m s⁻¹). The `earth_model` and `seismic_sph_signals` functions automatically convert to and from pyprop8's internal units (km, km s⁻¹, g cm⁻³).

## 4. Walkthrough: the verification test

This section steps through [`seismicsph_testcase.ipynb`](../verification_test/seismicsph_testcase.ipynb). The test is a small dam break: a 1 m long, 0.5 m high column of water collapses and strikes the far wall of a 3 m tank. The resolution is deliberately coarse to minimise computational cost, and the results are not intended for quantitative scientific interpretation.

### Step 1: Define the geometry and run the simulation

The simulation domain is a box with specified length, width and height. Candidate fluid and ice particles are generated throughout the domain. Two functions select the particles to retain: one defining the initial water volume and the other defining the ice boundaries. In this example, the ice forms the tank walls and floor, while the water initially occupies the first metre of the tank:

```python
import numpy as np
import hyglac as ss

def set_fluid(x, y, z, tol):
    xtol, ytol, ztol = tol, tol, 0.01 * tol
    return np.where((0.0 + xtol < x) & (x <= 1.0 + xtol)
                    & (-0.1 + ytol < y) & (y < 0.1 - ytol)
                    & (0.0 + ztol < z) & (z <= 0.5 + ztol))[0]

def set_obstacle(x, y, z, tol):
    xtol, ytol, ztol = tol, tol, 0.01 * tol
    return np.where((x <= 0.0 + xtol) | (0.1 - ytol <= y) | (y <= -0.1 + ytol)
                    | (x >= 3.0 + xtol) | (z <= 0.0 + ztol))[0]

ss.run_seismic_sph(
    3.0, 0.2, 1.0,          # tank length, width, height (m)
    1.5,                    # duration (s)
    0.025, 5e-5,           # particle spacing (m), initial timestep (s)
    set_fluid, set_obstacle,
    output_times=np.arange(0, 1.5 + 1. / 400, 1. / 200),  # 200 Hz output
    filename="scurve_testcase",
)
```

Snapshots are written to `scurve_testcase_output/`. Individual snapshots can be visualised using `ss.plot_seismic_sph(3.0, 0.2, 1.0, step=<iteration>, filename="scurve_testcase")`, where `<iteration>` corresponds to an existing output iteration.

### Step 2: extract the interface forces

```python
ss.seismic_sph_interactions(filename="scurve_testcase", delta=0.025,
                            dx=0.05, dy=0.05, dz=0.05, xrange=[2, 4])
```

This identifies the water–ice interface between `x = 2` and `x = 4` m, including the far wall at `x = 3` m. The force on each 5 cm cell of the interface is calculated at every output time and written to `scurve_testcase_interactions.hdf5`.

### Step 3: build the seismic force vectors

```python
forces_df = ss.seismic_sph_forces_defaults(3.0, 0.0, 500.0, filename="scurve_testcase")
forces_df.to_csv("scurve_testcase_forces.csv")
```

The interface forces are summed onto a 10 cm source grid, with the horizontal origin shifted to the far wall at `(x, y) = (3, 0)` m. The sources are placed 500 m below the surface, at the base of the ice layer in the Earth model defined below, and the vertical force component is set to zero.

### Step 4: compute seismograms

The Earth model consists of 500 m of ice overlying bedrock, as in Table 3 of the paper. Seismometers are positioned 0.5 m below the surface at four locations relative to the source-frame origin:

```python
import os

model = ss.earth_model([(0, 3600, 1800, 916), (500, 6000, 3500, 2700)])

os.makedirs("seismograms", exist_ok=True)
for xs, ys in [(500, 0), (1000, 0), (0, 500), (0, 1000)]:
    seismogram_df = ss.seismic_sph_signals(xs, ys, 0.5, earth_model=model, forces_df=forces_df,
                                           duration=5.001, sampling_rate=100.0)
    seismogram_df.to_csv(f"seismograms/scurve_testcase_seismogram_{xs}_{ys}_0.csv")
```

Each seismogram has three displacement components (`disp_1`-`disp_3`, m) and three velocity components (`vel_1`-`vel_3`, m s⁻¹). The notebook plots the velocity components and their power spectral densities, with the latter calculated using `scipy.signal.welch`.

## 5. API reference

### `run_seismic_sph`

```python
run_seismic_sph(x, y, z, t, delta, dt, set_fluid, set_obstacle, output_times=None, filename="seismic_sph")
```

Runs a three-dimensional weakly compressible SPH simulation of water in a rectangular tank. The user specifies the initial distribution of water, which starts from rest and evolves under gravity, and the geometry of the fixed ice boundaries.

| Parameter | Description |
| --- | --- |
| `x`, `y`, `z` | Tank length, width and height (m). |
| `t` | Simulation duration (s). |
| `delta` | Initial particle spacing δl (m). The SPH smoothing length is `h = 1.3 δl`. |
| `dt` | Initial timestep (s). The solver adapts the timestep as the simulation runs. |
| `set_fluid`, `set_obstacle` | Functions `f(x, y, z, tol)` that receive coordinate arrays of candidate particles and return the integer indices of particles to retain as initial water and ice, respectively. `tol = δl/4` is used to handle particles lying on geometry boundaries. |
| `output_times` | List or array of times (s) at which snapshots are written. Default: the final time only. |
| `filename` | Prefix for the output folder and files. Default `"seismic_sph"`. |

**Returns** `None`. Writes snapshots `<filename>_output/<filename>_<iteration>.hdf5` and the PySPH log `<filename>_output/<filename>.log`.

**Notes.** The tank's outer walls are generated as `boundary` particles. The `set_obstacle` function only needs to set the `obstacle` particles within the tank as both types are treated as ice. The fixed simulation settings are listed in [§7](#fixed-sph-settings).

### `fastload`

```python
fastload(path)
```

Loads particle data and simulation metadata from a PySPH output snapshot. The function provides access to the particle positions, smoothing lengths, densities and pressures for each particle type, together with information about the simulation state.

| Parameter | Description |
| --- | --- |
| `path` | Path to the PySPH output file containing the particle data and simulation metadata at a particular time. |

**Returns** a dictionary with keys `arrays` and `solver_data` for HDF5 snapshots. `arrays` contains an entry for each particle type (`fluid`, `boundary`, `obstacle`), with NumPy attributes `x`, `y`, `z`, `h`, `rho` and `p` (`p` is `None` if absent). For other file formats, the function calls PySPH's `load`.

### `plot_seismic_sph`

```python
plot_seismic_sph(x, y, z, step, particle_type="fluid", filename=None)
```

Plots the particle distribution from a selected simulation output, showing plan (`x-y`) and side (`x-z`) views. Fluid particles are coloured by density on a fixed scale of 990-1010 kg m⁻³.

| Parameter | Description |
| --- | --- |
| `x`, `y`, `z` | Tank length, width and height (m), used to set the plot limits. |
| `step` | Output iteration number identifying the simulation snapshot. |
| `particle_type` | Particle type to display: `"fluid"`, `"boundary"` or `"obstacle"`. Default `"fluid"`. |
| `filename` | Simulation prefix. Default `"seismic_sph"`. |

**Returns** `None`. Displays the particle distribution using Matplotlib.

**Notes.** The function modifies Matplotlib's global settings to use LaTeX text rendering and serif fonts, requiring a working LaTeX installation.

### `seismic_sph_interactions`

```python
seismic_sph_interactions(filename=None, delta=0.01, dx=0.05, dy=0.05, dz=0.05,
                         xrange=[4, 6], yrange=None, zrange=None, seed=0)
```

Identifies the water–ice interface from the SPH simulation output and calculates the pressure forces acting on the interface at each output time (paper §3.1–3.2).

| Parameter | Description |
| --- | --- |
| `filename` | Simulation prefix. Default (`None`): `"seismic_sph"`. |
| `delta` | Grid spacing δl (m) onto which ice particles are mapped to detect the interface. This should match the simulation's particle spacing; see [§7](#choosing-the-interface-parameters). Default `0.01`. |
| `dx`, `dy`, `dz` | Dimensions of the coarse grid cells Δl (m) used to sample the pressure and calculate the interface forces. Default `0.05`. |
| `xrange`, `yrange`, `zrange` | `[min, max]` windows in the simulation frame (m), used to restrict the calculation to the region of interest; `None` means no restriction. `xrange` defaults to `[4, 6]`, so `xrange=None` is required to include the whole tank. |
| `seed` | Seed for selecting one representative interface particle per cell. Default `0`. This sets NumPy's global random seed. |

**Returns** `None`. Writes or overwrites `<filename>_interactions.hdf5` in the working directory, containing the following datasets:

| Dataset | Contents |
| --- | --- |
| `t` | Output time (s). |
| `x`, `y`, `z` | Coarse-cell coordinates in the simulation frame (m), as multiples of `dx`, `dy`, `dz`. |
| `fx`, `fy`, `fz` | Force components (N), signed as described in [§3](#force-sign). |

The file contains one row per interface cell per output time. The initial output at `t = 0` is excluded.

### `seismic_sph_forces`

```python
seismic_sph_forces(x0, y0, z0, filename=None, dx=0.1, dy=0.1, dz=0.1, average_height=True, crit_force=1e-3)
```

Combines the calculated water–ice interface forces onto a coarser spatial grid to produce a time series of seismic force vectors. The source coordinates are transformed into the seismic source frame (paper §3.3).

| Parameter | Description |
| --- | --- |
| `x0`, `y0` | Simulation-frame position (m) that becomes the horizontal origin. |
| `z0` | Depth (m) of the simulation floor below the Earth-model surface. |
| `filename` | Simulation prefix. Default (`None`): `"seismic_sph"`. |
| `dx`, `dy`, `dz` | Source grid spacing (m). Default `0.1`. The grid can be finer than, equal to or coarser than the interactions grid. |
| `average_height` | If `True` (default), sum over all heights in each horizontal cell and place every source at depth `z0` (`loc_3 = −z0`). If `False`, retain separate vertical cells, with `loc_3 ≈ z − z0`. |
| `crit_force` | Exclude contributions weaker than this fraction of the maximum force. Default `1e-3`. The threshold is applied both before and after gridding. |

**Returns** a `pandas.DataFrame` with columns `t` (s), `loc_1`, `loc_2`, `loc_3` (m, source frame) and `force_1`, `force_2`, `force_3` (N). The DataFrame is not saved automatically but can be exported using `to_csv`.

**Notes.** The output interval is taken as the median spacing of the output times, and contributions more than 10⁻⁹ s from a regular time are dropped.

### `seismic_sph_forces_defaults`

```python
seismic_sph_forces_defaults(x0, y0, z0, filename=None)
```

Calculates seismic force vectors using the settings adopted for the channel geometries in Turner et al. (2026), including spatial restrictions and the exclusion of vertical forces.

| Parameter | Description |
| --- | --- |
| `x0`, `y0` | Simulation-frame position (m) that becomes the horizontal origin. |
| `z0` | Depth (m) of the simulation floor below the Earth-model surface. |
| `filename` | Simulation prefix. Default (`None`): `"seismic_sph"`. |

The function uses a source grid spacing of `0.1` m, with `average_height=True` and `crit_force=1e-3`. Only sources within `−1.001 ≤ loc_1 ≤ 1.001` m are retained, restricting the calculation to the channel bends, and the vertical force component (`force_3`) is set to zero.

**Returns** a `pandas.DataFrame` containing the seismic force vectors, with the same columns as `seismic_sph_forces`. These settings are specific to the channel geometries considered in the paper and may not be appropriate for other geometries.

### `earth_model`

```python
earth_model(layers)
```

Constructs a plane-layered Earth model for calculating seismic wave propagation.

| Parameter | Description |
| --- | --- |
| `layers` | List of `(top_depth, vp, vs, density)` tuples defining the seismic properties of each layer. Depths are in metres, seismic velocities in m s⁻¹ and densities in kg m⁻³. The first layer must begin at a depth of 0 m. |

**Returns** a `pyprop8.LayeredStructureModel`.

### `seismic_sph_signals`

```python
seismic_sph_signals(x, y, z, earth_model=None, forces_df=None, filename=None,
                    duration=10., sampling_rate=500.)
```

Calculates synthetic three-component seismograms at a specified seismometer location from a time series of seismic force vectors (paper §4).

| Parameter | Description |
| --- | --- |
| `x`, `y` | Seismometer position in the source frame (m), relative to `(x0, y0)`. |
| `z` | Seismometer depth (m, positive downward). |
| `earth_model` | Layered seismic velocity model generated using `earth_model`. Default: a uniform ice half-space (`vp = 3600` m s⁻¹, `vs = 1800` m s⁻¹, density 916 kg m⁻³). |
| `forces_df` | DataFrame containing the seismic force vectors. If omitted, the data are read from a CSV file specified by `filename`. |
| `filename` | Name of the force CSV file without the `.csv` extension. Used only when `forces_df` is not provided. |
| `duration` | Seismogram duration (s). Default `10`. Must exceed the final force time. |
| `sampling_rate` | Seismogram sampling frequency (Hz). Default `500`. Must not exceed the sampling frequency of the force vectors (see [§7](#choosing-the-seismic-parameters)). |

**Returns** a `pandas.DataFrame` with columns `t` (s), `disp_1`, `disp_2`, `disp_3` (m) and `vel_1`, `vel_2`, `vel_3` (m s⁻¹).

## 6. Choosing resolution and parameters

### Fixed SPH settings

These are module-level constants in `hyglac.py`. Edit them there if you need to change them.

| Setting | Value |
| --- | --- |
| Kernel, integrator | Wendland quintic, EPEC, adaptive timestep (`n_damp = 50`) |
| Smoothing length | `h = 1.3 δl` |
| Water density ρ₀ | 1000 kg m⁻³ |
| Equation-of-state exponent γ | 7 |
| Artificial viscosity α, β | 0.25, 0 |
| Sound speed c₀ | `10 √(2 × 9.81 × 0.55)` ≈ 32.8 m s⁻¹ |
| Gravity | −9.81 m s⁻² (vertical) |
| Ice (boundary) layers | 1 |

The sound speed is deliberately much lower than in real water so that larger timesteps can be used (paper §2.1). The value is inherited from PySPH's 3D dam-break example, which has a 0.55 m water column. Weakly compressible SPH usually needs c₀ to be about ten times the fastest flow speed to keep density changes to around 1%. The paper's flows reach about 5.5 m s⁻¹, so check the density spread with `plot_seismic_sph`, which shows 990–1010 kg m⁻³, when modelling faster or deeper flows.

### Choosing the simulation resolution

The paper uses a particle spacing of δl = 0.01 m. Runtime grows roughly as 1/δl³, because the number of particles does, so resolution is the main cost lever. The pressure-averaging radius (`7.5 h ≈ 9.75 δl`) should stay smaller than half the channel width. Narrow channels therefore need finer resolution (paper §6.2).

Write snapshots often enough to capture an impact in a single sample. The paper uses 200 Hz, i.e. `output_times` spaced 0.005 s apart.

### Choosing the interface parameters

`delta` should match the simulation's particle spacing. If `delta` is smaller than the spacing, neighbouring ice particles no longer land on adjacent grid points, so most ice particles look like interface particles and their normals become unreliable.

The coarse cell size (`dx`, `dy`, `dz`) only needs to be fine enough to tell apart different sections of the interface; the paper uses 0.04 m for 0.2 m wide channels. Wider channels can use coarser cells, which run faster.

### Choosing the seismic parameters

**Sample rate.** `sampling_rate` must not be higher than the force output rate (1/`dt`). Above about 2/`dt` the calculation fails outright. Ideally, choose a rate for which `1 / (sampling_rate × dt)` is a whole number, since the code rounds it.

**Source grid.** The source grid can be much coarser than the interface grid: at 200 Hz, an S wave in ice travels about 9 m between samples.

**Duration.** `duration` must be longer than the last force time plus the travel time to the seismometer.

**Surface waves.** These are well established only beyond roughly 0.5–0.75 km from the source in the paper's model (paper §4.1).

## 7. Using output from another SPH simulation

`seismic_sph_interactions` reads files directly, so another simulation's output must look like Hyglac's:

- **Snapshots:** a folder `<name>_output/` holding snapshots `<name>_<iteration>.hdf5`. HDF5 is required.
- **Particle groups:** in each snapshot, the groups `particles/fluid/arrays`, `particles/boundary/arrays` and `particles/obstacle/arrays`, each with datasets `x`, `y`, `z`, `h` and `rho`. The fluid group also needs `p`. The group names are fixed in the code, and the obstacle group must contain at least one particle.
- **Output times:** ideally the PySPH log `<name>_output/<name>.log`, which is the tested route. Without it, the code looks for a `solver_data/t` dataset in each snapshot; this fallback is untested.
- **Regular output times:** snapshots at irregular times are dropped.
- **Fixed ice:** the ice must not move, since only the first snapshot's ice positions are used.

Anything else requires adapting the file-reading functions in `hyglac.py`.

## 8. Limitations

**Physics**

- **Rigid ice:** the ice is fixed in place, so elastic response and resonance in the ice are not modelled.
- **No acoustics:** the reduced sound speed rules out acoustic waves in the water and conduit resonance.
- **Turbulence:** only turbulent structure larger than the SPH kernel is resolved.
- **Finite water volume:** the dam-break set-up gives a fixed amount of water and a short event. In the paper, the long-lasting forcing in the two-90°-bend case was extended beyond the simulation by repeating and damping it; that step is not part of `hyglac.py`.

**Implementation**

- **Single source depth:** the seismic stage uses the first source's depth for every source. Use `average_height=True`, or adapt the code, if you need several depths.
- **Spectra:** the `spec=True` output ignores force timing (see [§5](#seismic_sph_signals)).
- **Deep seismometers:** results for a seismometer deeper than the source are approximate.
- **Default `xrange`:** the default `[4, 6]` m in `seismic_sph_interactions` is easy to overlook and may miss the region of interest entirely.
- **Global settings:** `seismic_sph_interactions` resets NumPy's global random seed, and `plot_seismic_sph` changes Matplotlib's global settings.
- **Overwriting:** rerunning `seismic_sph_interactions` overwrites the interactions file. The force and seismogram functions never write files themselves.
