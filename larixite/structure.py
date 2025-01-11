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
from .amcsd_utils import PMG_CIF_OPTS
from .utils import get_color_logger

logger = get_color_logger()

if logger.level != 10:
    import warnings

    warnings.filterwarnings("ignore", category=UserWarning, module="pymatgen")


@dataclass
class XasStructureGroup:
    """ "Data container for an atomic structure model and absorber"""

    name: str  #: unique name, usually the input filename
    label: str  #: a short label, usually the input filename without extension
    filepath: Path  #: Path object to the input file
    struct: Structure  #: pymatgen Structure
    absorber: Element  #: pymatgen Element  #: pymatgen Element of the absorbing element


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
    assert len(structs_paths) == len(
        absorbers
    ), f"number of structures ({len(structs_paths)}) != number of absorbers ({len(absorbers)})"
    structs = []
    for istruct, struct_path in enumerate(structs_paths):
        struct = get_structure(struct_path, absorbers[istruct], **kwargs)
        logger.info(f"{istruct}: {struct.name}")
        structs.append(struct)
    return structs


