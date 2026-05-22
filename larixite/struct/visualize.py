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

    colors = [
        "red",
        "green",
        "blue",
        "orange",
        "yellow",
        "white",
        "purple",
        "pink",
        "brown",
        "cyan",
        "magenta",
        "olive",
        "navy",
        "teal",
        "maroon",
        "turquoise",
        "indigo",
        "salmon",
    ]
    color_elems = {}
    for idx, elem in enumerate(elems):
        color_elems[elem] = colors[idx]

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
    xyzview.setStyle(
        {"index": 0},
        {
            "stick": {
                "radius": 0.1,
                "opacity": 1,
                "hidden": False,
                "color": "gray",
            },
            "sphere": {"color": "gray", "radius": absorber_r, "opacity": 1},
        },
    )

    xyzview.zoomTo()
    xyzview.show()

    logger.info(color_elems)
    return xyzview
