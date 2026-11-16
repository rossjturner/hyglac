# Hyglac

**A computational framework to inform interpretation of hydrological processes in glacier environments through simulated seismic signals.**

Hyglac ([`hyglac.py`](hyglac.py)) simulates the seismic signals generated when moving water interacts with ice, as described in Turner et al. (2026). The modular framework has three key stages:

1. **Hydrodynamic simulation:** smoothed particle hydrodynamics (SPH) simulation of water flowing through an ice channel, using [PySPH](https://github.com/pypr/pysph).
2. **Seismic force vectors:** post-processing of the SPH output to find the water–ice interface, recover the smoothed dynamic pressure acting on it, and convert this into a time series of force vectors.
3. **Seismic wave propagation:** synthetic three-component seismograms for those forces in a layered ice and bedrock model, using [pyprop8](https://github.com/valentineap/pyprop8).

The stages are modular. Each one only reads the previous stage's output, so you can explore different seismometer positions or Earth models without repeating the expensive hydrodynamic simulation. You can also apply the force extraction to an existing SPH simulation, or compute seismograms from a force history produced elsewhere.

**Scope.** Hyglac targets relatively low-amplitude signals from fluid impacts, wave breaking and sloshing, and large-scale turbulent structures. The ice is treated as rigid and fixed, and the SPH sound speed is artificially reduced, so the framework does not model conduit resonance or acoustic waves in the water. Turbulence is resolved only on scales larger than the SPH kernel.

## Documentation

The [user guide](docs/user_guide.md) contains a walkthrough of the verification test, the full API reference, coordinate and sign conventions, advice on choosing resolution, and known limitations. The method is described in full in Turner et al. (2026).

## Installation

Hyglac is a single module, not an installable package. Create the Conda environment from `environment.yml` (Python 3.9, PySPH 1.0a6, pyprop8 1.1.5, NumPy, SciPy, pandas, h5py, Matplotlib), then run Python or Jupyter from the repository root:

```bash
conda env create -f environment.yml
conda activate pysph_env
```

```python
import hyglac as ss
```

`environment.yml` is a platform-specific export and may need adapting for other systems. Both PySPH and pyprop8 must be installed to import Hyglac, even if you only use one stage.

## Quick start

This runs the full pipeline with the verification-test settings. The geometry callbacks `set_fluid` and `set_obstacle` are defined in the [verification notebook](verification_test/seismicsph_testcase.ipynb).

```python
import numpy as np
import hyglac as ss

# 1. Hydrodynamic simulation: dam break in a 3 m x 0.2 m x 1 m tank, 1.5 s, outputs at 200 Hz
ss.run_seismic_sph(3.0, 0.2, 1.0, 1.5, 1/40, 5e-5, set_fluid, set_obstacle, output_times=np.arange(0, 1.5 + 0.0025, 0.005), filename="scurve_testcase")

# 2. Seismic force vectors
ss.seismic_sph_interactions(filename="scurve_testcase", delta=0.01,
                            dx=0.05, dy=0.05, dz=0.05, xrange=[2, 4])
forces_df = ss.seismic_sph_forces_defaults(3.0, 0.0, 500.0, filename="scurve_testcase")

# 3. Seismic wave propagation: 500 m of ice over bedrock, seismometer 1 km away at 0.5 m depth
model = ss.earth_model([(0, 3600, 1800, 916), (500, 6000, 3500, 2700)])
seismogram_df = ss.seismic_sph_signals(0, 1000, 0.5, earth_model=model, forces_df=forces_df, duration=5.001, sampling_rate=100.0)
```

| Stage | Function | Output |
| --- | --- | --- |
| Hydrodynamic simulation | `run_seismic_sph` | PySPH snapshots in `<name>_output/` |
| Interface forces | `seismic_sph_interactions` | `<name>_interactions.hdf5` |
| Seismic force vectors | `seismic_sph_forces`, `seismic_sph_forces_defaults` | Force-history DataFrame |
| Earth model | `earth_model` | `pyprop8.LayeredStructureModel` |
| Seismic wave propagation | `seismic_sph_signals` | Displacement and velocity DataFrame |

## Examples and verification test

Worked examples for the four subglacial channel geometries (two 45° bends, two 90° bends, a gentle 45° s-curve and a sharp 90° s-curve) and the Lobovský et al. (2014) dam-break validation case are provided in the corresponding Jupyter notebooks.

The [`seismicsph_testcase.ipynb`](verification_test/seismicsph_testcase.ipynb) notebook runs the complete workflow at low cost to check that the code works. Its parameters minimise runtime and output size so that reference outputs can be included in [`verification_test/`](verification_test/) for comparison. They are **not** physically realistic, and the outputs should not be interpreted as scientific predictions.

## Citation

If you use Hyglac, please cite the paper and the software archive:

- Turner, R. J., Magyar, J. C., Cook, S., and Reading, A. M. (2026). Hyglac: a computational framework to inform interpretation of hydrological processes in glacier environments through simulated seismic signals. In press.
- Turner, R. J. (2026). hyglac. Zenodo. https://doi.org/10.5281/zenodo.24161091

The Zenodo archive holds the exact version used in the paper, together with the scripts for all four channel geometries.

## Licence

MIT Licence.

## Contact

Ross Turner <turner.rj@icloud.com>
