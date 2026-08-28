#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
3D visualization of XasStructure local environment using py3Dmol
===============================================================

TODO (future improvements):

- Add coordination polyhedra around the absorber atom. py3Dmol does not support
  polyhedron rendering directly. Possible approaches:

  1. **nglview** — the `nglview` Jupyter widget supports `add_polyhedron()` and
     works directly with pymatgen structures. Could be offered as an alternative
     backend (e.g. `visualize(..., backend="nglview")`).

  2. **py3Dmol workaround with triangles** — compute the convex hull of the first
     coordination shell with `scipy.spatial.ConvexHull`, then draw the facets as
     semi-transparent plane surfaces via `py3Dmol`'s `addSurface()` or by adding
     dummy atoms at facet centroids with connecting sticks to form a wireframe.

  3. **pymatgen IPyWidgets** — `pymatgen.vis.ipyvolumetric` or `pymatgen.vis` can
     render structures with coordination environments in interactive HTML.

  4. **Export for VESTA** — write the cluster as a `.vesta` file (VESTA's native
     format supports polyhedra definitions) so the user can open it in VESTA.

- Expose the pymatgen Molecule cluster for downstream rendering.
- Support partial occupancy visualization (transparency, multi-snapshot).
- Allow custom color schemes beyond VESTA.
"""

import sys

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
    "Ac": ("#70abfa", "sky blue"),
    "Ag": ("#c0c0c0", "silver"),
    "Al": ("#81b2d6", "light grey blue"),
    "Am": ("#545cf2", "warm blue"),
    "Ar": ("#cffec4", "very light green"),
    "As": ("#74d057", "fresh green"),
    "At": ("#754f45", "brownish purple"),
    "Au": ("#ffd123", "sun yellow"),
    "B": ("#1fa20f", "green"),
    "Ba": ("#00c900", "vibrant green"),
    "Be": ("#5ed77b", "soft green"),
    "Bh": ("#e00038", "cherry"),
    "Bi": ("#9e4fb5", "purply"),
    "Bk": ("#8a4fe3", "purple"),
    "Br": ("#7e3102", "reddish brown"),
    "C": ("#4c4c4c", "charcoal"),
    "Ca": ("#5a96bd", "faded blue"),
    "Cd": ("#ffd98f", "wheat"),
    "Ce": ("#ffffc7", "ecru"),
    "Cf": ("#a136d4", "dark orchid"),
    "Cl": ("#31fc02", "vivid green"),
    "Cm": ("#785ce3", "purple"),
    "Co": ("#0000af", "royal blue"),
    "Cr": ("#00009e", "royal blue"),
    "Cs": ("#57178f", "indigo"),
    "Cu": ("#2247dc", "true blue"),
    "Db": ("#d1004f", "ruby"),
    "Dy": ("#1fffc7", "greenish cyan"),
    "Er": ("#00e675", "tealish green"),
    "Es": ("#b31fd4", "magenta purple"),
    "Eu": ("#61ffc7", "light greenish blue"),
    "F": ("#b0b9e6", "light steel blue"),
    "Fe": ("#b57100", "brown orange"),
    "Fm": ("#b31fba", "magenta purple"),
    "Fr": ("#420066", "royal purple"),
    "Ga": ("#9ee373", "spring green"),
    "Gd": ("#45ffc7", "greenish cyan"),
    "Ge": ("#7e6ea6", "dark lavender"),
    "H": ("#ffcccc", "light rose"),
    "He": ("#fce8ce", "blanched almond"),
    "Hf": ("#4dc2ff", "sky blue"),
    "Hg": ("#b8b8d0", "cloudy blue"),
    "Ho": ("#00ff9c", "medium spring green"),
    "Hs": ("#e6002e", "cherry red"),
    "I": ("#940094", "dark magenta"),
    "In": ("#a67573", "reddish grey"),
    "Ir": ("#175487", "light navy"),
    "K": ("#a121f6", "electric purple"),
    "Kr": ("#fac1f3", "light lavender"),
    "La": ("#5ac449", "leafy green"),
    "Li": ("#86df73", "light green"),
    "Lr": ("#c70066", "deep pink"),
    "Lu": ("#00ab24", "kelly green"),
    "Md": ("#b30da6", "magenta"),
    "Mg": ("#fb7b15", "pumpkin orange"),
    "Mn": ("#a7089d", "barney purple"),
    "Mo": ("#54b5b5", "medium aquamarine"),
    "Mt": ("#eb0026", "cherry red"),
    "N": ("#b0b9e6", "light steel blue"),
    "Na": ("#f9dc3c", "off yellow"),
    "Nb": ("#73c2c9", "seafoam blue"),
    "Nd": ("#c7ffc7", "very pale green"),
    "Ne": ("#fe37b5", "hot pink"),
    "Ni": ("#b7bbbd", "silver"),
    "No": ("#bd0d87", "medium violet red"),
    "Np": ("#0080ff", "deep sky blue"),
    "O": ("#fe0300", "red"),
    "Os": ("#266696", "ugly blue"),
    "P": ("#c09cc2", "pale purple"),
    "Pa": ("#00a1ff", "azure"),
    "Pb": ("#575961", "gunmetal"),
    "Pd": ("#006985", "peacock blue"),
    "Pm": ("#a3ffc7", "light seafoam"),
    "Po": ("#ab5c00", "orange brown"),
    "Pr": ("#d9ffc7", "very light green"),
    "Pt": ("#d0d0e0", "light grey"),
    "Pu": ("#006bff", "bright blue"),
    "Ra": ("#007d00", "green"),
    "Rb": ("#702eb0", "rebecca purple"),
    "Re": ("#267dab", "bluish"),
    "Rf": ("#cc0059", "deep pink"),
    "Rh": ("#0a7d8c", "ocean"),
    "Rn": ("#428296", "dirty blue"),
    "Ru": ("#248f8f", "blue green"),
    "S": ("#fffa00", "bright yellow"),
    "Sb": ("#9e63b5", "amethyst"),
    "Sc": ("#b563ab", "soft purple"),
    "Se": ("#9aef0f", "acid green"),
    "Sg": ("#d90045", "ruby"),
    "Si": ("#1b3bfa", "vivid blue"),
    "Sm": ("#8fffc7", "light seafoam"),
    "Sn": ("#9a8eb9", "heather"),
    "Sr": ("#00ff00", "lime"),
    "Ta": ("#4da6ff", "cornflower blue"),
    "Tb": ("#30ffc7", "greenish cyan"),
    "Tc": ("#3b9e9e", "sea"),
    "Te": ("#d47a00", "pumpkin"),
    "Th": ("#00baff", "deep sky blue"),
    "Ti": ("#78caff", "light blue"),
    "Tl": ("#a6544d", "light maroon"),
    "Tm": ("#00d452", "shamrock green"),
    "U": ("#008fff", "azure"),
    "V": ("#e51900", "tomato red"),
    "W": ("#2194d6", "water blue"),
    "Xe": ("#429eb0", "cool blue"),
    "Y": ("#94ffff", "robin egg blue"),
    "Yb": ("#00bf38", "shamrock green"),
    "Zn": ("#8f8f81", "warm grey"),
    "Zr": ("#00ff00", "lime"),
}


def _lighten_hex(hex_color: str, factor: float = 0.75) -> str:
    """Return a darker hex color obtained by scaling each RGB channel by *factor*."""
    r = int(int(hex_color[1:3], 16) * factor)
    g = int(int(hex_color[3:5], 16) * factor)
    b = int(int(hex_color[5:7], 16) * factor)
    return "#{:02x}{:02x}{:02x}".format(r, g, b)


def _print_color_legend(color_elems: dict[str, str], natoms: int) -> None:
    """Print a human-readable color legend, adapting to the execution environment.

    In Jupyter: renders HTML color swatches.
    In a terminal: uses ANSI 24-bit color escape codes.
    Fallback: plain text with element names and color descriptions.
    """
    try:
        from IPython.display import display, HTML

        swatches = " ".join(
            f'<span style="display:inline-block;width:14px;height:14px;border-radius:3px;background:{c};border:1px solid #555;vertical-align:middle;"></span>'
            f'<span style="color:{_adjust_font_color(c)};vertical-align:middle;">{e}</span>'
            for e, c in color_elems.items()
        )
        display(HTML(f'<div style="font-family:sans-serif;font-size:14px;padding:4px;">{natoms} atoms — {swatches}</div>'))
        return
    except ImportError:
        pass

    if hasattr(sys.stdout, "isatty") and sys.stdout.isatty() and sys.stdout.encoding == "UTF-8":
        for elem, hex_color in color_elems.items():
            r = int(hex_color[1:3], 16)
            g = int(hex_color[3:5], 16)
            b = int(hex_color[5:7], 16)
            sys.stdout.write("\x1b[38;2;{};{};{};m{}\x1b[0m ".format(r, g, b, elem))
        sys.stdout.write(f" ({natoms} atoms)\n")
        return

    legend = ", ".join(f"{e}: {VESTA_COLORS.get(e, ('#888888', 'grey'))[1]}" for e in color_elems)
    logger.info(f"Colors ({natoms} atoms): {legend}")


def _adjust_font_color(hex_color: str, dark: str = "#222222", light: str = "#ffffff") -> str:
    """Return dark or light text color for contrast against the given background hex color."""
    r = int(hex_color[1:3], 16)
    g = int(hex_color[3:5], 16)
    b = int(hex_color[5:7], 16)
    luminance = (0.299 * r + 0.587 * g + 0.114 * b) / 255
    return dark if luminance > 0.55 else light


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
        color_elems[elem] = VESTA_COLORS.get(elem, ("#888888", "grey"))[0]

    covalent_r = CovalentRadius().radius

    absorber_sym = struct.absorber.symbol
    absorber_color = color_elems.get(absorber_sym, "#888888")
    absorber_r = max(0.2, min(covalent_r.get(absorber_sym, 0.5) / 2, 0.75))
    absorber_light = _lighten_hex(absorber_color, 0.6)

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

    #: reset index 0 style to avoid duplication, then reapply with lighter color
    xyzview.setStyle({"index": 0}, None)  # clear
    xyzview.setStyle(
        {"index": 0},
        {
            "stick": {
                "radius": 0.1,
                "opacity": 1,
                "hidden": False,
                "color": absorber_light,
            },
            "sphere": {"color": absorber_light, "radius": absorber_r, "opacity": 1},
        },
    )

    xyzview.zoomTo()
    xyzview.show()

    _print_color_legend(color_elems, len(mol))
    return xyzview
