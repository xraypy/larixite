#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Atomic structure with an absorber element
==========================================
"""

import numpy as np
from pathlib import Path
from dataclasses import dataclass
from typing import Union, List, Literal
from pymatgen.core import Molecule, Structure, Element, Site
from larixite.utils import fcompact, get_logger, pprint

TOIMPLEMENT = "To implement as subclass, depending on the structure file format"
logger = get_logger("larixite.struct")


def site_label(site: Site) -> str:
    """
    return a string label for a pymatgen Site object,
    using the species string and fractional coordinates

    Parameters
    ----------
    site : pymatgen Site object

    Returns
    -------
    str
    """
    coords = ",".join([fcompact(s) for s in site.frac_coords])
    return f"{site.species_string}[{coords}]"


@dataclass
class XasStructure:
    """Data model for an atomic structure with an absorber element

    .. note::

        This class is agnostic of the file format for the atomic structure.
        Specific handling for different formats should be implemented in
        subclasses. The goal of this class is to provide a common API
        for manipulating atomic structures with an absorber element.

    """

    name: str  #: unique name, usually the input filename
    label: str  #: a short label, usually the input filename without extension
    filepath: Path  #: Path object to the input file
    struct: Structure  #: pymatgen Structure
    molecule: Molecule  #: pymatgen Molecule
    absorber: Element  #: pymatgen Element for the absorber
    absorber_idx: Union[int, None] = None  #: site index for the absorber
    radius: float = 7  #: radius of the absorption sphere from the absorbing atom
    radius_ext: float = 2.5  #: radius extension of the absorption sphere => cluster_size = radius + radius_ext
    struct_type: str = Literal["crystal", "molecule"]  #: type of the structure

    def __post_init__(self):
        if self.absorber_idx is None:
            self.absorber_idx = self.get_absorber_indexes()[0]
        self.build_sites()

    @property
    def cluster_size(self):
        """The size of the cluster is larger than the radius"""
        return self.radius + self.radius_ext

    @cluster_size.setter
    def cluster_size(self, value: float):
        self.radius = value - self.radius_ext

    @property
    def site_labels(self):
        """List all sites, e.g. `Fe[0.125,0.125,0.125]`"""
        return [site_label(site) for site in self.struct.sites]

    @property
    def formula(self):
        return self.struct.composition.reduced_formula

    @property
    def absorber_site(self):
        if self.absorber_idx is None:
            errmsg = "Absorber site index not set!"
            logger.error(errmsg)
            return IndexError(errmsg)
        return self.struct[self.absorber_idx]

    @property
    def cluster(self):
        return self.struct.get_sites_in_sphere(
            self.absorber_site.coords, self.cluster_size
        )

    @property
    def sga(self):
        raise NotImplementedError(TOIMPLEMENT)

    @property
    def space_group(self):
        raise NotImplementedError(TOIMPLEMENT)

    @property
    def sym_struct(self):
        raise NotImplementedError(TOIMPLEMENT)

    @property
    def wyckoff_symbols(self):
        raise NotImplementedError(TOIMPLEMENT)

    @property
    def equivalent_sites(self):
        raise NotImplementedError(TOIMPLEMENT)

    def get_site(self, site_index: int):
        return self.struct[site_index]

    def get_occupancy(self, species_string: str) -> float:
        """Get the occupancy of the absorbing atom from the species string"""
        try:
            ats_occ = species_string.split(",")
            at_occ = [at for at in ats_occ if self.absorber.name in at][0]
            occupancy = float(at_occ.split(":")[1])
        except Exception:
            occupancy = 1
        return occupancy

    def get_idx_in_struct(self, atom_coords):
        """Get the index corresponding to the given atomic coordinates (cartesian)"""
        for idx, atom in enumerate(self.struct):
            if np.allclose(atom.coords, atom_coords, atol=0.001) is True:
                return idx
        errmsg = f"atomic coordinates {atom_coords} not found in self.struct"
        logger.error(errmsg)
        # raise IndexError(errmsg)
        return None

    def get_absorber_indexes(self) -> List[int]:
        """Get the unique indexes of the absorbing element in the pymatgen structure (self.struct)"""
        absorber_indexes = []
        for sites in self.equivalent_sites:
            site = sites[0]
            if self.absorber.symbol in site.species_string:
                occupancy = self.get_occupancy(site.species_string)
                site_index = self.get_idx_in_struct(site.coords)
                if occupancy != 1:
                    logger.info(
                        f"[{self.label}] absorber {self.absorber.symbol} has occupancy of {occupancy} on site {site_index}"
                    )
                absorber_indexes.append(site_index)
        if len(absorber_indexes) == 0:
            errmsg = (
                f"Absorber {self.absorber.symbol} not found in structure {self.name}"
            )
            logger.error(errmsg)
            raise AttributeError(errmsg)
        return absorber_indexes

    def build_sites(self):
        """Parse the symmetrized structure to set the following attributes

        unique_sites:   list of non equivalent sites
        absorber_sites: sublist of unique_sites with the absorber

        .. note::

            content of the lists: (i, site, site_index_in_struct, site_occupancy, len(sites), self.wyckoff_symbols[i])

        """
        self.unique_sites = []
        self.absorber_sites = []
        absname = self.absorber.symbol
        for i, sites in enumerate(self.equivalent_sites):
            site = sites[0]
            site_index = self.get_idx_in_struct(site.coords)
            site_occupancy = self.get_occupancy(site.species_string)
            site_tuple = (
                i,
                site,
                site_index,
                site_occupancy,
                len(sites),
                self.wyckoff_symbols[i],
            )
            self.unique_sites.append(site_tuple)
            if absname in site.species_string:
                self.absorber_sites.append(site_tuple)

        if len(self.absorber_sites) == 0:
            errmsg = (
                f"Absorber {self.absorber.symbol} not found in structure {self.name}"
            )
            logger.error(errmsg)
            raise AttributeError(errmsg)

    def show_sites(self):
        """Show a tabular print for self.unique_sites"""
        header = [
            "idx",
            "label",
            "frac_coords",
            "idx_in_struct",
            "occupancy",
            "cart_coords",
            "wyckoff_site",
        ]
        infos = []
        for idx, site, site_index, occupancy, len_sites, wyckoff in self.unique_sites:
            if idx == self.absorber_idx:
                idx = f"{idx} (abs)"
            infos.append(
                [
                    idx,
                    site.label,
                    site.frac_coords,
                    site_index,
                    occupancy,
                    site.coords,
                    wyckoff,
                ]
            )
        matrix = [header]
        matrix.extend(infos)
        pprint(matrix)
