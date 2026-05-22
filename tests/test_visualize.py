"""Tests for larixite.struct.visualize — _cluster_to_xyz, _lighten_hex, VESTA_COLORS, visualize()."""

import math
from pathlib import Path

import pytest
import numpy as np

from larixite.struct import get_structure
from larixite.struct.visualize import (
    HAS_PY3DMOL,
    VESTA_COLORS,
    _round_up,
    _cluster_to_xyz,
    _lighten_hex,
    visualize,
)

STRUCTS = Path(__file__).parent / "structs"

# ── fixtures ──────────────────────────────────────────────────────


@pytest.fixture(scope="module")
def cif_zno():
    return get_structure(STRUCTS / "ZnO_mp-2133.cif", "Zn")


@pytest.fixture(scope="module")
def xyz_cuo6():
    return get_structure(STRUCTS / "CuO6_D4h.xyz", "Cu")


# ── _round_up ───────────────────────────────────────────────────


class TestRoundUp:
    def test_integer(self):
        assert _round_up(2.0) == 2.0

    def test_rounds_up(self):
        assert _round_up(2.31) == 2.31

    def test_rounds_up_small(self):
        assert _round_up(2.001) == 2.01

    def test_rounds_up_boundary(self):
        assert _round_up(2.501) == 2.51


# ── _cluster_to_xyz ─────────────────────────────────────────────


class TestClusterToXyz:
    def test_returns_strings_and_elems(self, cif_zno):
        mol = cif_zno.build_cluster(radius=2.5)
        xyz_str, elems = _cluster_to_xyz(mol)
        assert isinstance(xyz_str, str)
        assert isinstance(elems, list)
        assert all(isinstance(e, str) for e in elems)

    def test_xyz_header_matches_atom_count(self, cif_zno):
        mol = cif_zno.build_cluster(radius=2.5)
        xyz_str, elems = _cluster_to_xyz(mol)
        header = int(xyz_str.strip().split("\n")[0])
        assert header == len(mol)

    def test_absorber_first(self, cif_zno):
        mol = cif_zno.build_cluster(radius=2.5)
        xyz_str, _ = _cluster_to_xyz(mol)
        lines = xyz_str.strip().split("\n")
        absorber = lines[2].split()[0]
        assert absorber == cif_zno.absorber.symbol

    def test_coords_zero_origin(self, cif_zno):
        mol = cif_zno.build_cluster(radius=2.5)
        xyz_str, _ = _cluster_to_xyz(mol)
        lines = xyz_str.strip().split("\n")
        coords = lines[2].split()[1:]
        for c in coords:
            assert float(c) == pytest.approx(0.0, abs=1e-5)

    def test_elements_sorted(self, cif_zno):
        mol = cif_zno.build_cluster(radius=2.5)
        _, elems = _cluster_to_xyz(mol)
        assert elems == sorted(elems)

    def test_elements_match_structure(self, cif_zno):
        mol = cif_zno.build_cluster(radius=7.0)
        _, elems = _cluster_to_xyz(mol)
        assert "Zn" in elems
        assert "O" in elems

    def test_xyz_single_element(self, xyz_cuo6):
        mol = xyz_cuo6.build_cluster(radius=0.5)
        xyz_str, elems = _cluster_to_xyz(mol)
        assert len(elems) == 1
        assert elems == ["Cu"]


# ── _lighten_hex ────────────────────────────────────────────────


class TestLightenHex:
    def test_returns_hex(self):
        assert _lighten_hex("#ff0000", 0.5).startswith("#")

    def test_darkens_correctly(self):
        result = _lighten_hex("#ff0000", 0.5)
        # int(255 * 0.5) == 127 (truncation), not 128
        assert result == "#7f0000"

    def test_identity(self):
        assert _lighten_hex("#b57100", 1.0) == "#b57100"

    def test_factor_0_8(self):
        result = _lighten_hex("#ff0000", 0.8)
        r = int(int("ff", 16) * 0.8)
        assert "#" in result and len(result) == 7

    def test_fe_color(self):
        fe_hex, fe_name = VESTA_COLORS["Fe"]
        light = _lighten_hex(fe_hex, 0.8)
        assert light != fe_hex
        assert light.startswith("#")
        assert fe_name == "brown orange"

    def test_all_vesta_colors_tuplified(self):
        for sym, entry in VESTA_COLORS.items():
            assert isinstance(entry, tuple) and len(entry) == 2
            hex_color, name = entry
            light = _lighten_hex(hex_color, 0.8)
            assert light.startswith("#") and len(light) == 7
            assert isinstance(name, str) and name

    def test_invalid_hex_unchanged(self):
        with pytest.raises(ValueError):
            _lighten_hex("notahex", 0.5)


# ── VESTA_COLORS ────────────────────────────────────────────────


class TestVestaColors:
    COMMON = ["Fe", "O", "C", "N", "Si", "Al", "Ti", "Cu", "Zn", "S", "P", "Mg", "Ca"]

    def test_common_elements_present(self):
        for sym in self.COMMON:
            assert sym in VESTA_COLORS

    def test_colors_are_hex(self):
        for hex_color, name in VESTA_COLORS.values():
            assert hex_color.startswith("#") and len(hex_color) == 7
            assert isinstance(name, str) and name

    def test_all_vesta_values_valid_int(self):
        for hex_color, name in VESTA_COLORS.values():
            int(hex_color[1:3], 16)
            int(hex_color[3:5], 16)
            int(hex_color[5:7], 16)

    def test_no_black_in_palette(self):
        assert "#000000" not in [h for h, n in VESTA_COLORS.values()]

    def test_color_count(self):
        assert len(VESTA_COLORS) > 100


# ── visualize ───────────────────────────────────────────────────


class TestVisualize:
    @pytest.mark.skipif(not HAS_PY3DMOL, reason="py3Dmol not installed")
    def test_returns_view_cif(self, cif_zno):
        result = visualize(cif_zno, radius=2.5)
        assert result is not None

    @pytest.mark.skipif(not HAS_PY3DMOL, reason="py3Dmol not installed")
    def test_returns_view_xyz(self, xyz_cuo6):
        result = visualize(xyz_cuo6, radius=2.0)
        assert result is not None

    @pytest.mark.skipif(not HAS_PY3DMOL, reason="py3Dmol not installed")
    def test_default_radius(self, cif_zno):
        result = visualize(cif_zno)
        assert result is not None

    @pytest.mark.skipif(not HAS_PY3DMOL, reason="py3Dmol not installed")
    def test_unitcell_flag(self, cif_zno):
        result = visualize(cif_zno, radius=2.5, unitcell=True)
        assert result is not None

    @pytest.mark.skipif(not HAS_PY3DMOL, reason="py3Dmol not installed")
    def test_unitcell_xyz_no_error(self, xyz_cuo6):
        result = visualize(xyz_cuo6, radius=2.0, unitcell=True)
        assert result is not None

    @pytest.mark.skipif(not HAS_PY3DMOL, reason="py3Dmol not installed")
    def test_radius_affects_cluster(self, cif_zno):
        result = visualize(cif_zno, radius=1.5)
        assert result is not None

    @pytest.mark.skipif(HAS_PY3DMOL, reason="py3Dmol is installed")
    def test_no_py3dmol_returns_none(self, cif_zno, monkeypatch):
        import larixite.struct.visualize as viz
        monkeypatch.setattr(viz, "HAS_PY3DMOL", False)
        monkeypatch.setattr(viz, "py3Dmol", None)
        result = viz.visualize(cif_zno)
        assert result is None
