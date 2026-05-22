#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Coordination environment analysis for XAS structures
====================================================
"""

import numpy as np
from larixite.struct.xas import XasStructure
from larixite.utils import get_logger, pprint

logger = get_logger("larixite.struct.analyze")

# Mapping of chemenv coordination environment symbols to human-readable names
# Based on pymatgen chemenv ideal polyhedra definitions
CE_SYMBOL_DESCRIPTION = {
    "L:2": "Linear",
    "T:3": "Trigonal planar",
    "B:3": "Bent",
    "T:4": "Tetrahedral",
    "S:4": "Square planar",
    "F:4": "See-saw",
    "T:5": "Trigonal bipyramidal",
    "S:5": "Square pyramidal",
    "B:5": "Bay",
    "R:5": "Rocket",
    "G:5": "Envelope",
    "E:5": "Eclipsed pentagonal",
    "S:6": "Hexagonal planar",
    "A:6": "Hexagonal antiprismatic",
    "T:6": "Trigonal prismatic",
    "O:6": "Octahedral",
    "Y:6": "Pentagonal pyramidal",
    "M:6": "Pentagonal bisphenoidal",
    "H:6": "Triangular cup",
    "S:7": "Pentagonal bipyramidal",
    "D:7": "Heptagonal planar",
    "S:8": "Square antiprismatic",
    "D:8": "Dodecahedral",
    "A:8": "Hexagonal bipyramidal",
    "A:9": "Monocapped square antiprismatic",
    "A:10": "Bicapped square antiprismatic",
    "S:12": "Icosahedral",
}


def _group_into_shells(dists, toler):
    """Group sorted distances into shells.

    Distances within `toler` Angstrom of each other are considered
    belonging to the same coordination shell.

    Parameters
    ----------
    dists : list of float
        Sorted distances from the absorber.
    toler : float
        Tolerance [Angstrom] to group distances into the same shell.

    Returns
    -------
    list of list of float
        Each inner list contains the distances belonging to one shell.
    """
    if len(dists) == 0:
        return []

    shells = [[dists[0]]]
    for d in dists[1:]:
        if d - shells[-1][-1] <= toler:
            shells[-1].append(d)
        else:
            shells.append([d])

    return shells


def get_coord_env(
    struct: XasStructure, shell: int = 1, toler: float = 0.1
) -> list:
    """Determine the coordination environment around the absorber atom.

    For crystal structures (CIF), uses pymatgen's chemenv toolkit
    (LocalGeometryFinder, BVAnalyzer, MultiWeightsChemenvStrategy) to
    identify coordination geometry symbols such as:

        - S:4  - Square Plane
        - T:4  - Tetrahedral
        - T:5  - Trigonal bipyramid
        - S:5  - Square pyramidal
        - O:6  - Octahedral
        - T:6  - Trigonal prism

    For molecular structures (XYZ), uses pymatgen's CrystalNN to
    identify neighboring atoms.

    Neighbors are grouped into coordination shells according to their
    distances from the absorber: distances differing by less than
    `toler` are considered the same shell.

    Parameters
    ----------
    struct : XasStructure
        The XAS structure object.
    shell : int, optional
        Number of coordination shells to return. Default is 1 (first
        coordination shell only).
    toler : float, optional
        Distance tolerance [Angstrom] to group neighbors into the same
        coordination shell. Default is 0.1.

    Returns
    -------
    list
        List of lists, each containing:

            [0]: Info string about which site is being analyzed.
            [1]: Coordination environment dictionary with keys:
                - ce_symbol: coordination geometry symbol
                - ce_fraction: probability for given coordination env (0-1)
                - csm: continuous symmetry measure (lower = more symmetric)
                - permutation: atom arrangement order
            [2]: List of shell dictionaries, one per requested shell, with keys:
                - shell_idx: shell index (1-based)
                - avg_dist: average distance of atoms in the shell [Angstrom]
                - cn: coordination number (number of atoms in the shell)
                - elems: dict mapping element symbol to count in this shell
                - neighbors: list of neighbor site dicts with keys:
                    - site: pymatgen Site with element and coordinates
                    - site_index: structure index of the neighbor atom
                    - dist: distance from absorber [Angstrom]
                    - image: displacement vector from original site (CIF)
                    - weight: significance of the neighbor (CIF)
        """
    idx_abs_site = struct.absorber_idx
    abs_site = struct.struct[idx_abs_site]
    coord_env_ce = {}

    # Determine geometry symbol for the 1st coordination shell
    if struct.struct_type == "crystal":
        try:
            from pymatgen.analysis.bond_valence import BVAnalyzer
            from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import (
                LocalGeometryFinder,
            )
            from pymatgen.analysis.chemenv.coordination_environments.structure_environments import (
                LightStructureEnvironments,
            )
            from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import (
                MultiWeightsChemenvStrategy,
            )

            lgf = LocalGeometryFinder()
            lgf.setup_structure(struct.struct)

            bva = BVAnalyzer()
            try:
                valences = bva.get_valences(structure=struct.struct)
            except ValueError:
                valences = "undefined"

            se = lgf.compute_structure_environments(
                max_cn=6,
                valences=valences,
                only_indices=[idx_abs_site],
                only_symbols=["S:4", "T:4", "T:5", "S:5", "O:6", "T:6"],
            )

            logger.debug(se.voronoi.neighbors_distances[idx_abs_site][0]["max"])

            strategy = MultiWeightsChemenvStrategy.stats_article_weights_parameters()
            lse = LightStructureEnvironments.from_structure_environments(
                strategy=strategy, structure_environments=se
            )
            coord_env_ce = lse.coordination_environments[idx_abs_site]
        except Exception:
            logger.warning(
                "chemenv analysis failed, falling back to basic neighbor search"
            )
            coord_env_ce = {"ce_symbol": "chemenv analysis failed"}

    elif struct.struct_type == "molecule":
        try:
            from pymatgen.analysis.local_env import CrystalNN

            obj = CrystalNN()
            coord_dict = obj.get_cn_dict(struct.struct, idx_abs_site)
            coord_env_ce = {"ce_symbol": f"Elements Dict = {coord_dict}"}
        except Exception:
            logger.warning(
                "CrystalNN analysis failed, falling back to basic neighbor search"
            )
            coord_env_ce = {"ce_symbol": "CrystalNN analysis failed"}

    # Build full neighbor list using struct.get_neighbors() up to cluster_size
    # This works for arbitrary distances, unlike chemenv/CrystalNN which only
    # cover the 1st coordination shell
    search_radius = struct.cluster_size
    sphere = struct.struct.get_neighbors(abs_site, search_radius)

    neighbors = []
    for item in sphere:
        # Structure: (PeriodicNeighbor, distance, site_index, jimage)
        # Molecule:   (PeriodicNeighbor, distance, site_index)
        site = item[0]
        dist = round(float(item[1]), 5)
        si = item[2] if len(item) > 2 else None
        if dist < 1e-6:
            continue
        neighbors.append({"site": site, "dist": dist, "site_index": si})

    neighbors.sort(key=lambda x: x["dist"])

    # Group into shells by distance clustering
    dists = [nb["dist"] for nb in neighbors]
    dist_shells = _group_into_shells(dists, toler)

    # Build shell dicts for requested number of shells
    shell_dicts = []
    idx = 0
    for s in range(min(shell, len(dist_shells))):
        shell_dists = dist_shells[s]
        shell_len = len(shell_dists)
        shell_nbs = neighbors[idx : idx + shell_len]
        idx += shell_len

        avg_dist = round(float(np.mean(shell_dists)), 5)
        elems = {}
        for nb in shell_nbs:
            sym = nb["site"].species_string
            elems[sym] = elems.get(sym, 0) + 1

        shell_dicts.append(
            {
                "shell_idx": s + 1,
                "avg_dist": avg_dist,
                "cn": shell_len,
                "elems": elems,
                "neighbors": shell_nbs,
            }
        )

    coord_env_list = [
        [
            f"Coord. Env. for Site {struct.absorber_sites[0][0]}",
            coord_env_ce,
            shell_dicts,
        ]
    ]
    return coord_env_list


def show_coord_env(struct: XasStructure, **kws):
    """Pretty-print coordination environment information for the absorber atom.

    Calls get_coord_env() and formats the results as a table with shells,
    distances, coordination numbers and element breakdown.

    Parameters
    ----------
    struct : XasStructure
        The XAS structure object.
    **kws
        Additional keyword arguments passed to get_coord_env():
        shell : int, default 1
        toler : float, default 0.1
    """
    coord_env = get_coord_env(struct, **kws)[0]
    shells = coord_env[2]

    if len(shells) == 0:
        logger.warning("Empty coordination environment")
        return None

    # Extract coordination geometry symbol(s)
    ce_info = coord_env[1]
    if isinstance(ce_info, list):
        # chemenv returns a list, extract symbol from the first entry
        coord_sym = ce_info[0]["ce_symbol"] if ce_info else "not determined"
    else:
        # fallback: direct dict with ce_symbol key
        coord_sym = ce_info.get("ce_symbol", "not determined")

    print(
        f"Coord. Env. from absorber atom: {struct.absorber.symbol} "
        f"at site {struct.absorber_sites[0][0]}"
    )
    print(f"Geometry: {coord_sym}")
    print()

    header = ["Shell", "Dist. (Ang)", "CN", "Elements"]
    matrix = [header]

    for sh in shells:
        elems_str = ", ".join(
            f"{sym}: {cnt}" for sym, cnt in sorted(sh["elems"].items())
        )
        matrix.append(
            [
                sh["shell_idx"],
                sh["avg_dist"],
                sh["cn"],
                elems_str,
            ]
        )

    pprint(matrix)
    print()

    # Detailed neighbor table per shell
    for sh in shells:
        print(f"Shell {sh['shell_idx']} (avg dist: {sh['avg_dist']} A, CN: {sh['cn']})")
        det_header = ["Element", "Distance (Ang)"]
        det_matrix = [det_header]
        for nb in sh["neighbors"]:
            det_matrix.append([nb["site"].species_string, nb["dist"]])
        pprint(det_matrix)
        print()
