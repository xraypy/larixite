#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
3D visualization of XasStructure local environment using py3Dmol
===============================================================
"""

import math

from pymatgen.analysis.local_env import CovalentRadius

from larixite.struct.xas import XasStructure
from larixite.utils import get_logger

logger = get_logger("larixite.struct")

try:
    import py3Dmol

    HAS_PY3DMOL = True
except ImportError:
    HAS_PY3DMOL = False


VESTA_COLORS = {
    "Ac": "#70abfa",
    "Ag": "#c0c0c0",
    "Al": "#81b2d6",
    "Am": "#545cf2",
    "Ar": "#cffec4",
    "As": "#74d057",
    "At": "#754f45",
    "Au": "#ffd123",
    "B": "#1fa20f",
    "Ba": "#00c900",
    "Be": "#5ed77b",
    "Bh": "#e00038",
    "Bi": "#9e4fb5",
    "Bk": "#8a4fe3",
    "Br": "#7e3102",
    "C": "#4c4c4c",
    "Ca": "#5a96bd",
    "Cd": "#ffd98f",
    "Ce": "#ffffc7",
    "Cf": "#a136d4",
    "Cl": "#31fc02",
    "Cm": "#785ce3",
    "Co": "#0000af",
    "Cr": "#00009e",
    "Cs": "#57178f",
    "Cu": "#2247dc",
    "Db": "#d1004f",
    "Dy": "#1fffc7",
    "Er": "#00e675",
    "Es": "#b31fd4",
    "Eu": "#61ffc7",
    "F": "#b0b9e6",
    "Fe": "#b57100",
    "Fm": "#b31fba",
    "Fr": "#420066",
    "Ga": "#9ee373",
    "Gd": "#45ffc7",
    "Ge": "#7e6ea6",
    "H": "#ffcccc",
    "He": "#fce8ce",
    "Hf": "#4dc2ff",
    "Hg": "#b8b8d0",
    "Ho": "#00ff9c",
    "Hs": "#e6002e",
    "I": "#940094",
    "In": "#a67573",
    "Ir": "#175487",
    "K": "#a121f6",
    "Kr": "#fac1f3",
    "La": "#5ac449",
    "Li": "#86df73",
    "Lr": "#c70066",
    "Lu": "#00ab24",
    "Md": "#b30da6",
    "Mg": "#fb7b15",
    "Mn": "#a7089d",
    "Mo": "#54b5b5",
    "Mt": "#eb0026",
    "N": "#b0b9e6",
    "Na": "#f9dc3c",
    "Nb": "#73c2c9",
    "Nd": "#c7ffc7",
    "Ne": "#fe37b5",
    "Ni": "#b7bbbd",
    "No": "#bd0d87",
    "Np": "#0080ff",
    "O": "#fe0300",
    "Os": "#266696",
    "P": "#c09cc2",
    "Pa": "#00a1ff",
    "Pb": "#575961",
    "Pd": "#006985",
    "Pm": "#a3ffc7",
    "Po": "#ab5c00",
    "Pr": "#d9ffc7",
    "Pt": "#d0d0e0",
    "Pu": "#006bff",
    "Ra": "#007d00",
    "Rb": "#702eb0",
    "Re": "#267dab",
    "Rf": "#cc0059",
    "Rh": "#0a7d8c",
    "Rn": "#428296",
    "Ru": "#248f8f",
    "S": "#fffa00",
    "Sb": "#9e63b5",
    "Sc": "#b563ab",
    "Se": "#9aef0f",
    "Sg": "#d90045",
    "Si": "#1b3bfa",
    "Sm": "#8fffc7",
    "Sn": "#9a8eb9",
    "Sr": "#00ff00",
    "Ta": "#4da6ff",
    "Tb": "#30ffc7",
    "Tc": "#3b9e9e",
    "Te": "#d47a00",
    "Th": "#00baff",
    "Ti": "#78caff",
    "Tl": "#a6544d",
    "Tm": "#00d452",
    "U": "#008fff",
    "V": "#e51900",
    "W": "#2194d6",
    "Xe": "#429eb0",
    "Y": "#94ffff",
    "Yb": "#00bf38",
    "Zn": "#8f8f81",
    "Zr": "#00ff00",
}


def _round_up(x: float) -> float:
    return math.ceil(x * 100) / 100


def _cluster_to_xyz(mol) -> tuple[str, list[str]]:
    """Convert a pymatgen Molecule to an XYZ format string and the list of unique elements."""
    coords = []
    elements = []
    for site in mol:
        elem = site.species_string
        cart = site.coords
        coords.append((elem, cart))
        if elem not in elements:
            elements.append(elem)

    elements = sorted(elements)

    output_str = f"{len(coords)}\n\n"
    for elem, coord in coords:
        coord_str = " ".join(f"{c:.6f}" for c in coord)
        output_str += f"{elem} {coord_str}\n"

    return output_str, elements


def visualize(struct: XasStructure, radius: float = 2.5, unitcell: bool = False):
    """
    Display a 3D visualization of the local structure around the absorber atom.

    Args:
        struct: XasStructure object with the absorber defined.
        radius: radius in Angstroms for the visualization sphere around the central atom.
        unitcell: if True, displays the structure unit cell (CIF structures only).

    Returns:
        py3Dmol.view object, or None if py3Dmol is not installed.
    """
    if not HAS_PY3DMOL:
        logger.error("py3Dmol not installed! -> run `pip install py3Dmol`")
        return

    radius = _round_up(radius)

    mol = struct.build_cluster(radius=radius)

    xyz, elems = _cluster_to_xyz(mol)

    xyzview = py3Dmol.view(width=600, height=600)
    xyzview.addModel(xyz, "xyz")

    if unitcell:
        try:
            lat = struct.struct.lattice
            a, b, c = lat.a, lat.b, lat.c
            alpha, beta, gamma = lat.alpha, lat.beta, lat.gamma
            m = xyzview.getModel()
            m.setCrystData(a, b, c, alpha, beta, gamma)
            xyzview.addUnitCell()
        except AttributeError:
            logger.warning("Unit cell visualization not available for this structure type")

    color_elems = {}
    for elem in elems:
        color_elems[elem] = VESTA_COLORS.get(elem, "#888888")

    covalent_r = CovalentRadius().radius

    for elem in elems:
        color = color_elems[elem]
        sphere_r = max(0.2, min(covalent_r.get(elem, 0.5) / 2, 0.75))
        xyzview.setStyle(
            {"elem": elem},
            {
                "stick": {
                    "radius": 0.1,
                    "opacity": 1,
                    "hidden": False,
                    "color": color,
                },
                "sphere": {"color": color, "radius": sphere_r, "opacity": 1},
            },
        )

    absorber_sym = struct.absorber.symbol
    absorber_r = max(0.2, min(covalent_r.get(absorber_sym, 0.5) / 2, 0.75))
    absorber_color = color_elems.get(absorber_sym, "#888888")
    xyzview.setStyle(
        {"index": 0},
        {
            "stick": {
                "radius": 0.1,
                "opacity": 1,
                "hidden": False,
                "color": absorber_color,
            },
            "sphere": {"color": absorber_color, "radius": absorber_r, "opacity": 0.3},
        },
    )

    xyzview.zoomTo()
    xyzview.show()

    logger.info(color_elems)
    return xyzview
