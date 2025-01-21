#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Atomic structure with an absorber for XAS
================================================================
"""

from pathlib import Path
from dataclasses import dataclass
from typing import Union, Literal
from pymatgen.core import Molecule, Structure, Element, Site
from larixite.utils import fcompact


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
    """Generic data container for an atomic structure model with an absorber"""

    name: str  #: unique name, usually the input filename
    label: str  #: a short label, usually the input filename without extension
    filepath: Path  #: Path object to the input file
    file_format: Literal["cif", "xyz"]  #: input file format (supported only)
    struct: Structure  #: pymatgen Structure
    molecule: Molecule  #: pymatgen Molecule
    absorber: Element  #: pymatgen Element for the absorber
    absorber_idx: Union[int, None] = None  #: site index for the absorber
    radius: float = 7  #: radius of the absorption sphere from the absorbing atom

    def __post_init__(self):
        pass

    @property
    def cluster_size(self):
        """The size of the cluster is 2.5 Angstrom larger than the radius"""
        return self.radius + 2.5

    @cluster_size.setter
    def cluster_size(self, value: float):
        self.radius = value - 2.5

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
