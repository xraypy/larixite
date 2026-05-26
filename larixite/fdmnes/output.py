#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Search FDMNES job directories.
"""

from pathlib import Path

def search_jobs(globstr: str, jobsdir: Path | str | None = None) -> list[Path]:
    """Search for FDMNES job directories matching a glob pattern.

    Scans the given parent directory for subdirectories whose names match
    the provided glob expression.

    Parameters
    ----------
    globstr : str
        Glob pattern applied to the *name* of each immediate subdirectory
        (e.g. `"job_*"`, `"mi_*"`, `"*.done"`).
    jobsdir : Path | str | None -> ~/.larixite/fdmnes
        Parent directory containing job subdirectories.

    Returns
    -------
    list[Path]
        Sorted list of matching subdirectory paths.
    """
    if jobsdir is None:
        jobsdir = Path("~/.larixite/fdmnes").expanduser().resolve()
    if isinstance(jobsdir, str):
        jobsdir = Path(jobsdir).expanduser().resolve()
    if not jobsdir.is_dir():
        return []

    return sorted(
        (entry for entry in jobsdir.glob(globstr) if entry.is_dir()),
        key=lambda x: x.stat().st_ctime,
    )
