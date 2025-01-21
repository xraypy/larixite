#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Specialize XasStructure to handle structures from XYZ files
=============================================================
"""

from dataclasses import dataclass
import numpy as np
from larixite.struct.xas import XasStructure



@dataclass
class XasStructureXyz(XasStructure):

    def __post_init__(self):
        super().__post_init__()
        if self.absorber_idx is None:
            self.absorber_idx = self.get_absorber_sites()[0][0]
    @property
    def sga(self):
        raise AttributeError("SpacegroupAnalyzer fails for XYZ files")

    @property
    def space_group(self):
        return "P1"

    @property
    def sym_struct(self):
        raise ArithmeticError("Symmetrized structure not available for XYZ files")

    def get_idx_in_struct(self, atom_coords):
        """Get the index corresponding to the given atomic coordinates (cartesian)"""
        for idx, atom in enumerate(self.struct):
            if np.allclose(atom.coords, atom_coords, atol=0.001) is True:
                return idx
        errmsg = f"atomic coordinates {atom_coords} not found in self.struct"
        logger.error(errmsg)
        # raise IndexError(errmsg)
        return None

    def get_absorber_sites(self):
        """Get the indexes of the absorbing atoms in the structure"""
        absorber_sites = []
        for i, site in enumerate(self.struct.sites):
            if self.absorber.name in site.species_string:
                occupancy = self.get_occupancy(site.species_string)
                site_index = self.get_idx_in_struct(site.coords)
                if occupancy != 1:
                    logger.warning(
                        f"Absorber {self.absorber.name} has occupancy of {occupancy} on site {site_index}"
                    )
                absorber_sites.append(site_index)
        if len(absorber_sites) == 0:
            errmsg = f"Absorber {self.absorber.name} not found in structure {self.name}"
            logger.error(errmsg)
            raise AttributeError(errmsg)
        return absorber_sites
