"""Tests for larixite.struct.analyze — get_coord_env and show_coord_env."""

from pathlib import Path
from io import StringIO
import sys
import warnings

import pytest
import numpy as np

from larixite.struct import get_structure
from larixite.struct.analyze import get_coord_env, show_coord_env

STRUCTS = Path(__file__).parent / "structs"


# ── helpers ───────────────────────────────────────────────────────


def cap_show(output, struct, **kws):
    """Capture the stdout produced by show_coord_env and return it."""
    old_out = sys.stdout
    sys.stdout = buf = StringIO()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        output = show_coord_env(struct, **kws)
    sys.stdout = old_out
    return buf.getvalue()


@pytest.fixture(scope="module")
def cif_zno():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        return get_structure(STRUCTS / "ZnO_mp-2133.cif", "Zn")


@pytest.fixture(scope="module")
def xyz_cuo6():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        return get_structure(STRUCTS / "CuO6_D4h.xyz", "Cu")


# ── CIF (crystal) tests ─────────────────────────────────────────


class TestCoordEnvCIF:
    def test_shell1_returned(self, cif_zno):
        ce = get_coord_env(cif_zno, shell=1)
        shells = ce[0][2]
        assert len(shells) == 1
        assert shells[0]["shell_idx"] == 1

    def test_shell1_cn_positive(self, cif_zno):
        ce = get_coord_env(cif_zno, shell=1)
        shells = ce[0][2]
        assert shells[0]["cn"] > 0

    def test_shell1_has_elements(self, cif_zno):
        ce = get_coord_env(cif_zno, shell=1)
        shells = ce[0][2]
        assert len(shells[0]["elems"]) > 0

    def test_shell1_neighbors_sorted(self, cif_zno):
        ce = get_coord_env(cif_zno, shell=1)
        shells = ce[0][2]
        dists = [nb["dist"] for nb in shells[0]["neighbors"]]
        assert dists == sorted(dists)

    def test_multiple_shells_returned(self, cif_zno):
        ce = get_coord_env(cif_zno, shell=3)
        shells = ce[0][2]
        assert len(shells) >= 2

    def test_multiple_shells_increasing_distance(self, cif_zno):
        ce = get_coord_env(cif_zno, shell=3)
        shells = ce[0][2]
        for i in range(1, len(shells)):
            assert shells[i]["avg_dist"] > shells[i - 1]["avg_dist"]

    def test_cn_sum_matches_requested_shells(self, cif_zno):
        shell = 3
        ce = get_coord_env(cif_zno, shell=shell)
        shells = ce[0][2]
        expected_cn = 0
        distances = sorted(nb["dist"] for s in shells for nb in s["neighbors"])
        assert len(distances) == sum(s["cn"] for s in shells)

    def test_show_coord_env_output(self, cif_zno):
        output = cap_show(None, cif_zno, shell=3)
        assert "Coord. Env. from absorber atom: Zn" in output
        assert "Shell" in output
        assert "CN" in output

    def test_toler_groups_closely_spaced_distances(self, cif_zno):
        ce_strict = get_coord_env(cif_zno, shell=10, toler=0.01)
        ce_loose = get_coord_env(cif_zno, shell=10, toler=1.0)
        n_shells_strict = len(ce_strict[0][2])
        n_shells_loose = len(ce_loose[0][2])
        # Larger tolerance should group into fewer (or equal) shells
        assert n_shells_loose <= n_shells_strict


# ── XYZ (molecule) tests ─────────────────────────────────────────


class TestCoordEnvXYZ:
    def test_shell1_returned(self, xyz_cuo6):
        ce = get_coord_env(xyz_cuo6, shell=1)
        shells = ce[0][2]
        assert len(shells) >= 1
        assert shells[0]["shell_idx"] == 1

    def test_shell1_cn_positive(self, xyz_cuo6):
        ce = get_coord_env(xyz_cuo6, shell=1)
        shells = ce[0][2]
        assert shells[0]["cn"] > 0

    def test_shell1_has_elements(self, xyz_cuo6):
        ce = get_coord_env(xyz_cuo6, shell=1)
        shells = ce[0][2]
        assert len(shells[0]["elems"]) > 0

    def test_shell1_neighbors_sorted(self, xyz_cuo6):
        ce = get_coord_env(xyz_cuo6, shell=1)
        shells = ce[0][2]
        dists = [nb["dist"] for nb in shells[0]["neighbors"]]
        assert dists == sorted(dists)

    def test_show_coord_env_output(self, xyz_cuo6):
        output = cap_show(None, xyz_cuo6, shell=2)
        assert "Coord. Env. from absorber atom: Cu" in output
        assert "Shell" in output
        assert "CN" in output
