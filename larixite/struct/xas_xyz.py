#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Specialize XasStructure to handle structures from XYZ files
=============================================================
"""

from dataclasses import dataclass
from larixite.struct.xas import XasStructure
from larixite.utils import get_logger

logger = get_logger("larixite.struct")


@dataclass
class XasStructureXyz(XasStructure):
    @property
    def struct(self):
        return self.molecule

    @property
    def sga(self):
        raise AttributeError("SpacegroupAnalyzer fails for XYZ files")

    @property
    def space_group(self):
        return "P1"

    @property
    def sym_struct(self):
        logger.warning("No symmetrized structure for XYZ, returning the molecule")
        return self.molecule

    @property
    def wyckoff_symbols(self):
        return ["1a" for site in self.molecule.sites]

    @property
    def equivalent_sites(self):
        return [[site] for site in self.molecule.sites]
