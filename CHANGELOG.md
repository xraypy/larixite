# Changelog

Notable changes to this project will be documented in this file. The GitHub
Release Notes will also be useful.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

The versioning scheme followed is `YYYY.MAJOR.MINOR`:

- `YYYY`: year of release
- `MAJOR`: major feature (= change in API, no backward compatibility)
- `MINOR`: minor feature (= bug fixes, minor changes not affecting the API, backward compatible)

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
