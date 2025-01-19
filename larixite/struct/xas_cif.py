#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Specialize XasStructure to handle structures from CIF files
=============================================================
"""

from dataclasses import dataclass
import numpy as np
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from larixite.struct.xas import XasStructure


@dataclass
class XasStructureCif(XasStructure):
    @property
    def sga(self):
        return SpacegroupAnalyzer(self.struct)

    @property
    def space_group(self):
        return self.sga.get_symmetry_dataset().international

    @property
    def sym_struct(self):
        return self.sga.get_symmetrized_structure()

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
        for i, sites in enumerate(self.sym_struct.equivalent_sites):
            site = sites[0]
            if self.absorber.name in site.species_string:
                occupancy = self.get_occupancy(site.species_string)
                site_index = self.get_idx_in_struct(site.coords)
                if occupancy != 1:
                    logger.warning(
                        f"Absorber {self.absorber.name} has occupancy of {occupancy} on site {site_index}"
                    )
                absorber_sites.append([site_index, self.struct[site_index]])
        if len(absorber_sites) == 0:
            errmsg = f"Absorber {self.absorber.name} not found in structure {self.name}"
            logger.error(errmsg)
            raise AttributeError(errmsg)
        return absorber_sites
