# Changelog

Notable changes to this project will be documented in this file. The GitHub
Release Notes will also be useful.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

The versioning scheme followed is `YYYY.MAJOR.MINOR`:

- `YYYY`: year of release
- `MAJOR`: major feature (= change in API, no backward compatibility)
- `MINOR`: minor feature (= bug fixes, minor changes not affecting the API, backward compatible)

## [2026.2.0 - not released yet]

### `fdmnes`

    - Split `fdmnes.py` into a `fdmnes` package (`fdmnes/input.py`) with refactored input handling; `FdmnesXasInput` and `struct2fdmnes` remain importable from `larixite.fdmnes` as before
    - Direct `Cif_file` structure type support: a CIF can now be passed straight to FDMNES instead of only an expanded `Crystal` description
    - Added convolution parameters (`E_cut`, `Dec`, `Ecent`, `Elarg`, `Gamma_hole`, `Gamma_max`, `Gaussian`, `Estart`) and automatic `Conv_out` section in the generated input
    - `green`/`scf` are now plain attributes kept in sync with `params["Green"]`/`params["SCF"]` (previously computed properties)
    - Added `run_sbatch()` to submit jobs via `sbatch --parsable` and track status in `status.yaml`
    - Added `dump_params()` / `from_yaml()` to store and restore a `FdmnesXasInput` configuration as YAML
    - Extended `write_sbatch()` with `nnodes`, `mem_per_cpu`, `walltime`, `constraint`, `partition` options (default `ncpus` raised from 8 to 12)
    - Removed the `get_vmax()`/`get_rself()` helper methods; `vmax`/`rself` are now rendered directly from their attribute values (also fixes `Vmax` previously always being written as `-6` regardless of the value set)

### `struct`

    - Migration of coordination environment analysis and visualization from `struct2xas` into separate modules `analyze` and `visualize` (new `py3Dmol`-based 3D visualization, new `pyproject.toml` dependency)

### Bug fixes

    - Fix `FdmnesXasInput(optimize=True)` not actually enabling `SCF` in the generated input despite logging "SCF enabled" (regression introduced during the `fdmnes` refactor, caught before release; affects the web app's FDMNES "optimize" option)
    - Restore `from larixite.fdmnes import logger`, broken by the module split (used in `examples/fdmnes_xas_workflow.ipynb`)

### Tests

    - Added tests for analysis, visualization, and fdmnes input modules

## [2026.1.0 - 18-May-2026]

### web application

    - Separate FEFF and FDMNES sections for XAS input generation
    - FDMNES: added optimization checkbox
    - FDMNES: added Green formalism checkbox
    - Added more tests

### `fdmnes`

    - Update to FDMNES version March 2026
    - Add `make_sbatch` with a template for the ESRF SLURM scheduler
    - Add `optimize : bool = False` option in `FdmnesXasInput` class to control if the input parameters are optimized (= more accurate calcuation, but slower) or not
    - `write_input()` defaults to `$HOME/.larixite/fdmnes/{job_prefix}_{YYMMDD_HHMMSS}`

### `feff`

    - add more commented-out Feff cars with parameter hints
    - add atom count and default valus for `l` in the ipot section
    - fix exchange card to include `vr` and `vi`.
    - better handling of clusters larger than 8 Ang.

### `struc2xas`

    - Fully migrated from `larch.xrd.struc2xas` (and removed there)
    - Make `struc2xas` *deprecated* and start refactoring
