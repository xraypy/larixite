#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Search FDMNES job directories.
"""

from dataclasses import dataclass
from pathlib import Path

from larixite.fdmnes.input import FdmnesXasInput
from larixite.utils import get_logger

logger = get_logger("larixite.fdmnes.output")
logger.propagate = False


@dataclass
class FdmnesXasSim:
    """A FDMNES XAS simulation associated with a job directory."""

    jobdir: Path
    input: FdmnesXasInput | None = None

    @property
    def prefix(self) -> str | None:
        if self.input is not None:
            return self.input.fileout_prefix
        return None


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
