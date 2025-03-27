#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Specialize XasStructure to handle structures from CIF files
=============================================================
"""

from dataclasses import dataclass
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from larixite.struct.xas import XasStructure
from larixite.utils import get_logger

logger = get_logger("larixite.struct")


@dataclass
class XasStructureCif(XasStructure):
    @property
    def struct(self):
        return self.structure

    @property
    def sga(self):
        return SpacegroupAnalyzer(self.struct)

    @property
    def space_group(self):
        spg = self.sga.get_symmetry_dataset()
        return f"{spg.international}:{spg.choice}"

    @property
    def sym_struct(self):
        return self.sga.get_symmetrized_structure()

    @property
    def wyckoff_symbols(self):
        return self.sym_struct.wyckoff_symbols

    @property
    def equivalent_sites(self):
        return self.sym_struct.equivalent_sites
