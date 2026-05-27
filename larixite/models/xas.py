#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Experimental XAS data model."""

from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from larixite.utils import get_logger

logger = get_logger("larixite.models.xas")


@dataclass
class XasExpData:
    """Model for an Experimental XAS spectrum."""

    data: pd.DataFrame | None = None

    @classmethod
    def from_csv(cls, csvpath: Path | str) -> "XasExpData":
        """Load experimental XAS data from a CSV file.

        Parameters
        ----------
        csvpath : Path | str
            Path to the CSV file containing the XAS data.

        Returns
        -------
        XasExpData
            Instance with the data loaded from the CSV.
        """
        if isinstance(csvpath, str):
            csvpath = Path(csvpath)
        data = pd.read_csv(csvpath)
        return cls(data=data)
