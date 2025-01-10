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
from pymatgen.core import Molecule, Structure, Element, Lattice
from pymatgen.io.xyz import XYZ
from pymatgen.io.cif import CifParser
from pymatgen.ext.matproj import MPRestError
import json
from .amcsd_utils import PMG_CIF_OPTS
from .utils import get_color_logger

logger = get_color_logger()

if logger.level != 10:
    import warnings

    warnings.filterwarnings("ignore", category=UserWarning, module="pymatgen")


@dataclass
class StructureGroup:
    name: str
    label: str
    filepath: Path
    struct: Structure


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


def get_structure(filepath: Union[str, Path], frame: int = 0) -> Structure:
    """Get a pymatgen Structure from CIF/XYZ file"""
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
            struct = structs.parse_structures()[0]
        except Exception:
            raise ValueError("could not get structure from text of CIF")
        logger.debug("structure created from a CIF file")
        return struct
    #: XYZ
    if filepath.suffix == ".xyz":
        xyz = XYZ.from_file(filepath)
        molecules = xyz.all_molecules
        mol = molecules[frame]
        struct = mol2struct(mol)
        logger.debug("structure created from a XYZ file")
        return struct
    #: UNSUPPORTED
    raise ValueError(f"file type {filepath.suffix} not supported")


def get_structs_from_dir(
    structsdir: Union[str, Path],
    globstr: str = "*",
    exclude_names: list[str] = None,
    **kwargs,
) -> list[StructureGroup]:
    """Get a list of Structures from a directory"""
    if isinstance(structsdir, str):
        structsdir = Path(structsdir)
    structs_paths = list(structsdir.glob(globstr))
    if exclude_names is not None:
        structs_paths = [
            struct for struct in structs_paths if struct.name not in exclude_names
        ]
    structs = []
    for istruct, struct_path in enumerate(structs_paths):
        struct = get_structure(struct_path, **kwargs)
        logger.info(f"{istruct}: {struct_path.name}")
        structs.append(
            StructureGroup(
                name=struct_path.name,
                label=struct_path.stem,
                filepath=struct_path,
                struct=struct,
            )
        )
    return structs
