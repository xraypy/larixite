#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Search FDMNES job directories.
"""

from dataclasses import dataclass
from pathlib import Path

import numpy as np
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
        if self.metadata is None:
            self.load()
        return self.metadata.get("VO_interstitial", 0)

    @property
    def e0(self) -> float:
        if self.metadata is None:
            self.load()
        return self.metadata.get("E_edge", 0)

    @property
    def energy(self, erel: bool = False) -> np.ndarray:
        if self.data is None:
            self.load()
        eout = self.data["Energy"].values
        return eout + self.e0 if erel else eout

    @property
    def mu(self) -> np.ndarray:
        if self.data is None:
            self.load()
        return self.data["<xanes>"].values

    def load(self) -> pd.DataFrame | None:
        """Load the FDMNES output data file and parse the metadata line."""
        datafile = self.jobdir / f"{self.prefix}.txt"
        if not datafile.exists():
            error_msg = "cannot find the calculation output file"
            logger.error(error_msg)
            error_file = self.jobdir / "error.yaml"
            error_file.write_text(yaml.dump({"error": error_msg}))
            return

        lines = datafile.read_text().splitlines()
        skiprows = 1

        for i, line in enumerate(lines):
            stripped = line.strip()
            if stripped.startswith("#"):
                continue
            elif "=" in stripped:
                eq_split = stripped.split(" = ", 1)
                left_vals = [float(v) for v in eq_split[0].split() if v]
                right_names = [n.strip() for n in eq_split[1].split(",")]
                assert len(left_vals) == len(right_names), (
                    "Error parsing the metadata line"
                )
                int_names = {"Z", "n_edge", "j_edge", "ninit", "ninit1", "natomsym_f"}
                self.metadata = {
                    name: int(val) if name in int_names else val
                    for name, val in zip(right_names, left_vals)
                }
                skiprows = i + 1
            else:
                break

        self.data = pd.read_csv(datafile, sep=r"\s+", skiprows=skiprows)
        assert len(self.data.columns) >= 2, "Output file must have at least two columns"
        assert "Energy" in self.data.columns, "No 'Energy' column in the output file"
        assert "<xanes>" in self.data.columns, "No '<xanes>' column in the output file"


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
