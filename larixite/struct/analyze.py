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


def get_coord_env(struct: XasStructure, radius: float = 3) -> list:
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

    Parameters
    ----------
    struct : XasStructure
        The XAS structure object.
    radius : float, optional
        Maximum distance [Angstrom] to search for neighboring atoms.
        Default is 3.

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
            [2]: List of neighbor site dictionaries with keys:
                - site: pymatgen Site with element and coordinates
                - site_index: structure index of the neighbor atom
                - image: displacement vector from original site
                - weight: significance of the neighbor to the coordination
        """
    idx_abs_site = struct.absorber_idx
    coord_env_list = []

    if struct.struct_type == "crystal":
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
        coord_env = lse._all_nbs_sites

    elif struct.struct_type == "molecule":
        from pymatgen.analysis.local_env import CrystalNN

        obj = CrystalNN()

        coord_env = obj.get_nn_info(struct.struct, idx_abs_site)
        for site in coord_env:
            site["site"].cart_coords = struct.struct[site["site_index"]].coords
        coord_dict = obj.get_cn_dict(struct.struct, idx_abs_site)
        coord_env_ce = {"ce_symbol": f"Elements Dict = {coord_dict}"}

    else:
        raise ValueError(f"Unknown structure type: {struct.struct_type}")

    coord_env_list.append(
        [
            f"Coord. Env. for Site {struct.absorber_sites[0][0]}",
            coord_env_ce,
            coord_env,
        ]
    )
    return coord_env_list


def show_coord_env(struct: XasStructure, **kws):
    """Pretty-print coordination environment information for the absorber atom.

    Calls get_coord_env() and formats the results as a table with element
    species and distances from the absorber site.

    Parameters
    ----------
    struct : XasStructure
        The XAS structure object.
    **kws
        Additional keyword arguments passed to get_coord_env().
    """
    coord_env = get_coord_env(struct, **kws)[0]
    abs_site_coord = struct.absorber_sites[0][4]

    elems_dist = []
    sites = coord_env[2]

    if len(sites) == 0:
        logger.warning("Empty coordination environment")
        return None

    if struct.struct_type == "crystal":
        coord_sym = [
            coord_env[1][i]["ce_symbol"] for i in range(len(coord_env[1]))
        ]
        for site in sites:
            elems_dist.append(
                (
                    site["site"].species,
                    round(np.linalg.norm(site["site"].coords - abs_site_coord), 5),
                )
            )
    elif struct.struct_type == "molecule":
        coord_sym = coord_env[1]["ce_symbol"]
        for site in sites:
            elems_dist.append(
                (
                    site["site"].species,
                    round(
                        np.linalg.norm(site["site"].cart_coords - abs_site_coord), 5
                    ),
                )
            )

    elems_dist = sorted(elems_dist, key=lambda x: x[1])
    print(
        f"Coord. Env. from absorber atom: {struct.absorber.symbol} "
        f"at site {struct.absorber_sites[0][0]}"
    )
    print(coord_sym)
    header = ["Element", "Distance"]
    matrix = [header]
    matrix.extend(elems_dist)
    pprint(matrix)
