#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Search FDMNES job directories.
"""

from dataclasses import dataclass
from pathlib import Path

import yaml
import pandas as pd
from larixite.fdmnes.input import FdmnesXasInput
from larixite.utils import get_logger

logger = get_logger("larixite.fdmnes.output")


@dataclass
class FdmnesXasSim:
    """A FDMNES XAS simulation associated with a job directory."""

    jobdir: Path
    input: FdmnesXasInput | None = None
    data: pd.DataFrame | None = None
    metadata: dict | None = None

    @property
    def prefix(self) -> str | None:
        if self.input is not None:
            return self.input.fileout_prefix
        return None

    @property
    def efermi(self) -> float:
        if self.metadata is not None:
            return self.metadata.get("VO_interstitial", 0)
        logger.warning("metadata not available, run `load_data()` first")
        return 0

    @property
    def e0(self) -> float:
        if self.metadata is not None:
            return self.metadata.get("E_edge", 0)
        logger.warning("metadata not available, , run `load_data()` first")
        return 0

    def load_data(self, shift_energy: bool = False) -> pd.DataFrame | None:
        datafile = self.jobdir / f"{self.prefix}.txt"
        if not datafile.exists():
            error_msg = "cannot find the calculation output file"
            logger.error(error_msg)
            error_file = self.jobdir / "error.yaml"
            error_file.write_text(
                yaml.dump({"error": error_msg})
            )
            return

        lines = datafile.read_text().splitlines()
        e_shift: float = 0.0
        skiprows = 1

        for i, line in enumerate(lines):
            stripped = line.strip()
            if stripped.startswith("#"):
                continue
            elif "=" in stripped:
                eq_split = stripped.split(" = ", 1)
                left_vals = [float(v) for v in eq_split[0].split() if v]
                right_names = [n.strip() for n in eq_split[1].split(",")]
                assert len(left_vals) == len(right_names), "Error parsing the metadata line"
                int_names = {"Z", "n_edge", "j_edge", "ninit", "ninit1", "natomsym_f"}
                self.metadata = {
                    name: int(val) if name in int_names else val
                    for name, val in zip(right_names, left_vals)
                }
                e_shift = self.metadata.get("E_edge", 0.0)
                skiprows = i + 1
                print(skiprows)
            else:
                break

        self.data = pd.read_csv(datafile, sep=r"\s+", skiprows=skiprows)
        if shift_energy:
            self.data.iloc[:, 0] += e_shift

def search_jobs(globstr: str, jobsdir: Path | str | None = None) -> list[FdmnesXasSim]:
    """Search for FDMNES job directories matching a glob pattern.

    For each matching subdirectory that contains a ``*_params.yaml`` file,
    a :class:`FdmnesXasSim` instance is built via
    :meth:`FdmnesXasInput.from_yaml`.

    Parameters
    ----------
    globstr : str
        Glob pattern applied to the *name* of each immediate subdirectory
        (e.g. ``"job_*"``, ``"mi_*"``, ``"*.done"``).
    jobsdir : Path | str | None -> ~/.larixite/fdmnes
        Parent directory containing job subdirectories.

    Returns
    -------
    list[FdmnesXasSim]
        Sorted list of matching simulation objects.
    """
    if jobsdir is None:
        jobsdir = Path("~/.larixite/fdmnes").expanduser().resolve()
    if isinstance(jobsdir, str):
        jobsdir = Path(jobsdir).expanduser().resolve()
    if not jobsdir.is_dir():
        return []

    results: list[FdmnesXasSim] = []
    for entry in jobsdir.glob(globstr):
        if not entry.is_dir():
            continue

        for param_file in entry.glob("*_params.yaml"):
            try:
                inp = FdmnesXasInput.from_yaml(param_file)
                results.append(FdmnesXasSim(jobdir=entry, input=inp))
            except Exception:
                continue

    results = sorted(results, key=lambda s: s.jobdir.stat().st_ctime)
    logger.info(f"found {len(results)} job(s) in {jobsdir}")
    return results
