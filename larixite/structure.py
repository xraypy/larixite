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
from typing import Union
from pymatgen.core import Molecule, Structure, Element, Lattice
from pymatgen.io.xyz import XYZ
from pymatgen.io.cif import CifParser
from pymatgen.ext.matproj import MPRestError
import json
from .amcsd_utils import PMG_CIF_OPTS
from .utils import get_color_logger

logger = get_color_logger()


def mol2struct(molecule: Molecule) -> Structure:
    """Convert pymatgen Molecule to a pymatgen Structure"""
    alat, blat, clat = 1, 1, 1
    #extend the lattice
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


def get_structure(file: Union[str, Path], frame: int = 0) -> Structure:
    """Get a pymatgen Structure from CIF/XYZ/MPJSON file"""
    if isinstance(file, str):
        file = Path(file)
    if not file.exists():
        errmsg = f"{file} not found"
        logger.error(errmsg)
        raise FileNotFoundError(errmsg)

    if file.suffix == ".cif":
        try:
            struct = Structure.from_file(file, **PMG_CIF_OPTS)
        except MPRestError:
            logger.error(f"error reading {file}")
            raise
        logger.debug("structure created from a CIF file")
    elif file.suffix == ".xyz":
        xyz = XYZ.from_file(file)
        molecules = xyz.all_molecules
        mol = molecules[frame]
        struct = mol2struct(mol)
        logger.debug("structure created from a XYZ file")
    elif file.suffix == ".mpjson":
        struct = Structure.from_dict(json.load(open(file, "r")))
        logger.debug("structure created from JSON file")
    else:
        errmsg = "only CIF, XYZ and MPJSON files are currently supported"
        logger.error(errmsg)
        raise NotImplementedError(errmsg)

    return struct
