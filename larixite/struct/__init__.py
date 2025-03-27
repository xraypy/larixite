#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Wrapper on top of pymatgen to handle atomic structures for XAS calculations
============================================================================
"""

import numpy as np
from pathlib import Path
from typing import Union
from copy import deepcopy
from pymatgen.io.xyz import XYZ
from pymatgen.io.cif import CifParser
from pymatgen.core import Molecule, Structure, Element, Lattice
from larixite.struct.xas import XasStructure
from larixite.struct.xas_cif import XasStructureCif
from larixite.struct.xas_xyz import XasStructureXyz
from larixite.utils import get_logger
from larixite.amcsd_utils import PMG_CIF_OPTS

logger = get_logger("larixite.struct")

if logger.level != 10:
    import warnings

    warnings.filterwarnings("ignore", category=UserWarning, module="pymatgen")


def mol2struct(molecule: Molecule) -> Structure:
    """Convert a pymatgen Molecule to Structure -> WARNING: not working as expected -> DO NOT USE!"""
    # extend the lattice
    cart_coords = deepcopy(molecule.cart_coords)
    species = deepcopy(molecule.species)
    # build lattice
    xcoords = np.array([coord[0] for coord in cart_coords])
    ycoords = np.array([coord[1] for coord in cart_coords])
    zcoords = np.array([coord[2] for coord in cart_coords])
    alat = xcoords.max() - xcoords.min()
    blat = ycoords.max() - ycoords.min()
    clat = zcoords.max() - zcoords.min()
    offset = np.array([abs(xcoords.min()), abs(ycoords.min()), abs(zcoords.min())])
    cart_coords += offset
    lattice = Lattice.from_parameters(
        a=alat, b=blat, c=clat, alpha=90, beta=90, gamma=90
    )
    struct = Structure(lattice, species, cart_coords, coords_are_cartesian=True)
    return struct


def get_structure(
    filepath: Union[str, Path], absorber: str, frame: int = 0
) -> XasStructure:
    """Get a subclass of XasStructure from a structural file, according to its format

    Parameters
    ----------
    filepath : str or Path
        Filepath to CIF/XYZ file
    absorber : str
        Atomic symbol of the absorbing element
    frame : int, optional
        Index of the structure for multi-frame structures in the CIF/XYZ file

    Returns
    -------
    XasStructure
        The XAS structure group for the specified file and absorber.
    """
    if isinstance(filepath, str):
        filepath = Path(filepath)
    if not filepath.exists():
        errmsg = f"{filepath} not found"
        logger.error(errmsg)
        raise FileNotFoundError(errmsg)
    structout = None
    #: CIF
    if filepath.suffix == ".cif":
        try:
            structs = CifParser(filepath, **PMG_CIF_OPTS)
        except Exception:
            raise ValueError(f"could not parse text of CIF from {filepath}")
        try:
            struct = structs.parse_structures(primitive=False)[frame]
        except Exception:
            raise ValueError(f"could not get structure {frame} from text of CIF")
        molecule = Molecule.from_dict(struct.as_dict())
        logger.debug("structure created from a CIF file")
        structout = XasStructureCif(
            name=filepath.name,
            label=filepath.stem,
            filepath=filepath,
            structure=struct,
            molecule=molecule,
            struct_type="crystal",
            absorber=Element(absorber),
        )
    #: XYZ
    if filepath.suffix == ".xyz":
        xyz = XYZ.from_file(filepath)
        molecules = xyz.all_molecules
        molecule = molecules[frame]
        structure = mol2struct(molecule)
        logger.debug("structure created from a XYZ file")
        structout = XasStructureXyz(
            name=filepath.name,
            label=filepath.stem,
            filepath=filepath,
            structure=structure,
            molecule=molecule,
            struct_type="molecule",
            absorber=Element(absorber),
        )
    #: some checks on the structure
    if not structout.struct.is_ordered:
        logger.warning(
            f"[{structout.name}] contains partially occupied sites that are not fully supported yet"
        )
    if structout is not None:
        return structout
    else:
        #: UNSUPPORTED
        raise ValueError(f"File type {filepath.suffix} not supported yet")


def get_structs_from_dir(
    structsdir: Union[str, Path],
    absorbers: Union[list[str], str],
    globstr: str = "*",
    exclude_names: list[str] = None,
    **kwargs,
) -> list[XasStructure]:
    """Get a list of XasStructure from a directory containing structural files

    Parameters
    ----------
    structsdir : str or Path
        directory containing the structural files
    absorbers : list of str or str
        list of atomic symbols of the absorbing elements for each file or a single symbol for all
    globstr : str, optional
        string to filter the files in the directory
    exclude_names : list of str, optional
        list of filenames to exclude
    **kwargs : dict, optional
        additional keyword arguments to pass to `get_structure`

    Returns
    -------
    list of XasStructure
        list of XasStructure objects

    Examples
    --------

    from pathlib import Path
    curdir = Path().cwd()
    basedir = curdir.parent
    testdir = basedir / "test"
    structsdir = testdir / "structs"
    abs = "Fe"
    structs = get_structs_from_dir(structsdir, abs, globstr=f"*{abs}*", exclude_names=["NAMING.tmpl"])

    """
    if isinstance(structsdir, str):
        structsdir = Path(structsdir)
    structs_paths = list(structsdir.glob(globstr))
    if exclude_names is not None:
        structs_paths = [
            struct for struct in structs_paths if struct.name not in exclude_names
        ]
    if isinstance(absorbers, str):
        absorbers = [absorbers] * len(structs_paths)
    assert len(structs_paths) == len(absorbers), (
        f"number of structures ({len(structs_paths)}) != number of absorbers ({len(absorbers)})"
    )
    structs = []
    for istruct, struct_path in enumerate(structs_paths):
        struct = get_structure(struct_path, absorbers[istruct], **kwargs)
        logger.info(f"{istruct}: {struct.name}")
        structs.append(struct)
    return structs
