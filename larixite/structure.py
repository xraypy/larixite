#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Wrapper on top of pymatgen to get structures from structural files
===================================================================

.. note::

    --> IMPORTANT: <--

    A structure here is intended composed of an unit cell with atomic positions in fractional coordinates.

"""

import numpy as np
from pathlib import Path
from dataclasses import dataclass
from typing import Union
from pymatgen.core import Molecule, Structure, Element, Lattice, Site
from pymatgen.io.xyz import XYZ
from pymatgen.io.cif import CifParser
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from .amcsd_utils import PMG_CIF_OPTS
from .utils import get_color_logger, fcompact

logger = get_color_logger("larixite.structure")

if logger.level != 10:
    import warnings

    warnings.filterwarnings("ignore", category=UserWarning, module="pymatgen")


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
class XasStructureGroup:
    """ "Data container for an atomic structure model and absorber"""

    name: str  #: unique name, usually the input filename
    label: str  #: a short label, usually the input filename without extension
    filepath: Path  #: Path object to the input file
    struct: Structure  #: pymatgen Structure
    absorber: Element  #: pymatgen Element  #: pymatgen Element of the absorbing element
    absorber_site: int = 1  #: site index with absorber
    radius: float = 7  #: radius of the calculation from the absorbing atom

    @property
    def cluster_size(self):
        """The size of the cluster is 2.5 Angstrom larger than the radius"""
        return self.radius + 2.5

    @cluster_size.setter
    def cluster_size(self, value: float):
        self.radius = value - 2.5


def mol2struct(molecule: Molecule) -> Structure:
    """Convert pymatgen Molecule to a pymatgen Structure"""
    alat, blat, clat = 1, 1, 1
    # extend the lattice
    alat, blat, clat = np.max(molecule.cart_coords, axis=0)
    lattice = Lattice.from_parameters(
        a=alat, b=blat, c=clat, alpha=90, beta=90, gamma=90
    )
    # Create a list of species
    species = [Element(sym) for sym in molecule.species]
    # Create a list of coordinates
    coords = molecule.cart_coords
    # Create the Structure object
    struct = Structure(lattice, species, coords, coords_are_cartesian=True)
    return struct


def get_structure(
    filepath: Union[str, Path], absorber: str, frame: int = 0
) -> XasStructureGroup:
    """Get a XasStructureGroup from a structural files


    :param filepath: filepath to CIF/XYZ file
    :param absorber: atomic symbol of the absorbing element
    :param frame: index of the structure in the CIF/XYZ file

    """
    if isinstance(filepath, str):
        filepath = Path(filepath)
    if not filepath.exists():
        errmsg = f"{filepath} not found"
        logger.error(errmsg)
        raise FileNotFoundError(errmsg)
    #: CIF
    if filepath.suffix == ".cif":
        try:
            structs = CifParser(filepath, **PMG_CIF_OPTS)
        except Exception:
            raise ValueError(f"could not parse text of CIF from {filepath}")
        try:
            struct = structs.parse_structures()[frame]
        except Exception:
            raise ValueError(f"could not get structure {frame} from text of CIF")
        logger.debug("structure created from a CIF file")
    #: XYZ
    elif filepath.suffix == ".xyz":
        xyz = XYZ.from_file(filepath)
        molecules = xyz.all_molecules
        mol = molecules[frame]
        struct = mol2struct(mol)
        logger.debug("structure created from a XYZ file")
    else:
        #: UNSUPPORTED
        raise ValueError(f"file type {filepath.suffix} not supported")
    return XasStructureGroup(
        name=filepath.name,
        label=filepath.stem,
        filepath=filepath,
        struct=struct,
        absorber=Element(absorber),
    )


def get_structs_from_dir(
    structsdir: Union[str, Path],
    absorbers: Union[list[str], str],
    globstr: str = "*",
    exclude_names: list[str] = None,
    **kwargs,
) -> list[XasStructureGroup]:
    """Get a list of XasStructureGroup from a directory"""
    if isinstance(structsdir, str):
        structsdir = Path(structsdir)
    structs_paths = list(structsdir.glob(globstr))
    if exclude_names is not None:
        structs_paths = [
            struct for struct in structs_paths if struct.name not in exclude_names
        ]
    if isinstance(absorbers, str):
        absorbers = [absorbers] * len(structs_paths)
    assert (
        len(structs_paths) == len(absorbers)
    ), f"number of structures ({len(structs_paths)}) != number of absorbers ({len(absorbers)})"
    structs = []
    for istruct, struct_path in enumerate(structs_paths):
        struct = get_structure(struct_path, absorbers[istruct], **kwargs)
        logger.info(f"{istruct}: {struct.name}")
        structs.append(struct)
    return structs


def build_sites(xsg: XasStructureGroup):
    """parse sites of the structure to get several components:

    struct.sites:   list of all sites as parsed by pymatgen
    site_labels:    list of site labels
    unique_sites:   list of (site[0], wyckoff sym) for unique xtal sites
    unique_map:     mapping of all site_labels to unique_site index
    absorber_sites: list of unique sites with absorber

    """
    # get equivalent sites, mapping of all sites to unique sites,
    # and list of site indexes with absorber

    xsg.formula = xsg.struct.composition.reduced_formula
    sga = SpacegroupAnalyzer(xsg.struct)
    xsg.space_group = sga.get_symmetry_dataset().international

    sym_struct = sga.get_symmetrized_structure()
    wyckoff_symbols = sym_struct.wyckoff_symbols

    xsg.site_labels = []
    for site in xsg.struct.sites:
        xsg.site_labels.append(site_label(site))

    xsg.unique_sites = []
    xsg.unique_map = {}
    xsg.absorber_sites = []
    absorber = xsg.absorber.name
    for i, sites in enumerate(sym_struct.equivalent_sites):
        xsg.unique_sites.append((sites[0], len(sites), wyckoff_symbols[i]))
        for site in sites:
            xsg.unique_map[site_label(site)] = i + 1
        if absorber in site.species_string:
            xsg.absorber_sites.append(i)

    xsg.atom_sites = {}
    xsg.atom_site_labels = {}

    for i, dat in enumerate(xsg.unique_sites):
        site = dat[0]
        label = site_label(site)
        for species in site.species:
            elem = species.name
            if elem in xsg.atom_sites:
                xsg.atom_sites[elem].append(i + 1)
                xsg.atom_site_labels[elem].append(label)
            else:
                xsg.atom_sites[elem] = [i + 1]
                xsg.atom_site_labels[elem] = [label]

    all_sites = {}
    for xat in xsg.atom_site_labels.keys():
        all_sites[xat] = {}
        for i, label in enumerate(xsg.atom_site_labels[xat]):
            all_sites[xat][label] = xsg.atom_sites[xat][i]
    xsg.all_sites = all_sites


def build_cluster(
    xsg: XasStructureGroup,
    absorber_site=None,
    radius=None,
):
    if absorber_site not in xsg.atom_sites[xsg.absorber]:
        raise ValueError(
            f"invalid site for absorber {absorber}: must be in {xsg.atom_sites[xsg.absorber]}"
        )
    if radius is not None:
        xsg.radius = radius
    cluster_size = xsg.cluster_size
    csize2 = cluster_size**2
    site_atoms = {}  # map xtal site with list of atoms occupying that site
    site_tags = {}

    for i, site in enumerate(xsg.struct.sites):
        label = site_label(site)
        s_unique = xsg.unique_map.get(label, 0)
        site_species = [e.symbol for e in site.species]
        if len(site_species) > 1:
            s_els = [s.symbol for s in site.species.keys()]

            s_wts = [s for s in site.species.values()]
            site_atoms[i] = rng.choices(s_els, weights=s_wts, k=1000)
            site_tags[i] = f"({site.species_string:s})_{s_unique:d}"
        else:
            site_atoms[i] = [site_species[0]] * 1000
            site_tags[i] = f"{site.species_string:s}_{s_unique:d}"

    # atom0 = xsg.struct[a_index]
    atom0 = xsg.unique_sites[absorber_site - 1][0]
    sphere = xsg.struct.get_neighbors(atom0, xsg.cluster_size)

    xsg.symbols = [xsg.absorber]
    xsg.coords = [[0, 0, 0]]
    site0_species = [e.symbol for e in atom0.species]
    if len(site0_species) > 1:
        xsg.tags = [f"({atom0.species_string})_{absorber_site:d}"]
    else:
        xsg.tags = [f"{atom0.species_string}_{absorber_site:d}"]

    for i, site_dist in enumerate(sphere):
        s_index = site_dist[0].index
        site_symbol = site_atoms[s_index].pop()

        coords = site_dist[0].coords - atom0.coords
        if (coords[0] ** 2 + coords[1] ** 2 + coords[2] ** 2) < csize2:
            xsg.tags.append(site_tags[s_index])
            xsg.symbols.append(site_symbol)
            xsg.coords.append(coords)

    xsg.molecule = Molecule(xsg.symbols, xsg.coords)
