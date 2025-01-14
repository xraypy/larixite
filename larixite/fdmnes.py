#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Generating FDMNES input files
==============================

[FDMNES](http://fdmnes.neel.cnrs.fr/) is a program for calculating X-ray
spectroscopy (XAS, XES, RIXS) from the atomic structures

"""

from dataclasses import dataclass
from pathlib import Path
from .structure import XasStructureGroup
from .utils import get_logger

logger = get_logger("larixite.fdmnes")

@dataclass
class FdmnesInput:
    xsg: XasStructureGroup  #: container for the atomic structure
    tmplpath: Path  #: path to the FDMNES input template