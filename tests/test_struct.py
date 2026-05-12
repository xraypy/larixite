"""Tests for larixite.struct — XasStructure, CIF/XYZ parsing, clusters, error handling."""

import tempfile
from pathlib import Path

import pytest
import numpy as np

from larixite.struct import get_structure, get_structure_from_text, get_structs_from_dir
from larixite.struct.xas import site_label
from larixite.struct.xas_cif import XasStructureCif
from larixite.struct.xas_xyz import XasStructureXyz

STRUCTS = Path(__file__).parent / "structs"

# ── fixtures ──────────────────────────────────────────────────────


@pytest.fixture(scope="module")
def cif_zno():
    with pytest.warns(UserWarning, match="fractional coordinates rounded"):
        return get_structure(STRUCTS / "ZnO_mp-2133.cif", "Zn")


@pytest.fixture(scope="module")
def cif_fe3o4():
    return get_structure(STRUCTS / "Fe3O4_cub_fracOcc_Levy2012_ICSD-183969.cif", "Fe")


@pytest.fixture(scope="module")
def xyz_cuo6():
    return get_structure(STRUCTS / "CuO6_D4h.xyz", "Cu")


# ── parametrized occupancy (legacy) ────────────────────────────

test_structures = (
    ("CuO6_D4h.xyz", "Cu", 0, 1),
    ("VO6_Oh.xyz", "V", 0, 1),
    ("GaBr_multi-frame.xyz", "Ga", 0, 1),
    ("Fe3O4_cub_fracOcc_Levy2012_ICSD-183969.cif", "Fe", 0, 0.5),
    ("Fe2.646Ti.354O4_magnetite_Bosi2009_AMCSD4820.cif", "Fe", 0, 0.823),
    ("Fe2.646Ti.354O4_magnetite_Bosi2009_AMCSD4820.cif", "Fe", 1, 1),
    ("K.87H6.13Fe2.79S2O14_jarosite_Basciano2007_AMCSD4438.cif", "Fe", 3, 1),
)


@pytest.mark.parametrize("filename,abs,abs_index,occupancy", test_structures)
def test_struct(filename, abs, abs_index, occupancy):
    sg = get_structure(STRUCTS / filename, abs)
    sg.absorber_idx = abs_index
    s = sg.get_site(sg.absorber_idx)
    assert sg.get_occupancy(s) == pytest.approx(occupancy, rel=0.02)


# ── CIF ────────────────────────────────────────────────────────


class TestCif:
    def test_properties(self, cif_zno):
        assert isinstance(cif_zno, XasStructureCif)
        assert cif_zno.struct_type == "crystal"
        assert cif_zno.formula == "ZnO"
        assert cif_zno.space_group == 186
        assert cif_zno.absorber_idx == 0

    def test_symmetry(self, cif_zno):
        assert cif_zno.sga is not None and cif_zno.sym_struct is not None
        assert len(cif_zno.unique_sites) == 2 and len(cif_zno.equivalent_sites) == 2
        assert len(cif_zno.wyckoff_symbols) == 2
        assert all(isinstance(w, str) and w for w in cif_zno.wyckoff_symbols)

    def test_site_labels(self, cif_zno):
        labels = cif_zno.site_labels
        assert len(labels) == 4 and all("[" in l and "]" in l for l in labels)

    def test_absorber(self, cif_zno):
        assert len(cif_zno.absorber_sites) >= 1
        idxes = cif_zno.get_absorber_indexes()
        assert len(idxes) >= 1 and all("Zn" in cif_zno.get_site(i).species_string for i in idxes)

    def test_alt_absorber(self, cif_zno):
        with pytest.warns(UserWarning, match="fractional coordinates rounded"):
            sg = get_structure(STRUCTS / "ZnO_mp-2133.cif", "O")
        assert sg.absorber.symbol == "O" and len(sg.absorber_sites) >= 1

    def test_full_occupancy(self, cif_zno):
        assert cif_zno.get_occupancy(cif_zno.get_site(0)) == 1.0

    def test_cluster(self, cif_zno):
        from pymatgen.core import Molecule
        cl = cif_zno.cluster
        assert isinstance(cl, Molecule) and cl[0].species_string == "Zn"
        assert np.allclose(cl[0].coords, [0, 0, 0], atol=0.01)
        assert len(cl) > len(cif_zno.struct)

    def test_cluster_radius(self, cif_zno):
        cif_zno.radius = 3
        small = cif_zno.cluster
        cif_zno.radius = 7
        assert len(small) < len(cif_zno.cluster)


# ── XYZ ────────────────────────────────────────────────────────


class TestXyz:
    def test_properties(self, xyz_cuo6):
        assert isinstance(xyz_cuo6, XasStructureXyz)
        assert xyz_cuo6.struct_type == "molecule" and xyz_cuo6.space_group == "P1"
        assert xyz_cuo6.absorber_idx == 0
        assert "Cu" in xyz_cuo6.formula and "O" in xyz_cuo6.formula

    def test_symmetry(self, xyz_cuo6):
        assert all(w == "1a" for w in xyz_cuo6.wyckoff_symbols)
        sites = xyz_cuo6.equivalent_sites
        assert len(sites) == len(xyz_cuo6.struct) and all(len(g) == 1 for g in sites)

    def test_misc(self, xyz_cuo6):
        assert len(xyz_cuo6.get_absorber_indexes()) >= 1
        assert xyz_cuo6.get_occupancy(xyz_cuo6.get_site(0)) == 1.0

    def test_multi_frame(self):
        for f in (0, 3):
            sg = get_structure(STRUCTS / "GaBr_multi-frame.xyz", "Ga", frame=f)
            assert isinstance(sg, XasStructureXyz) and sg.absorber.symbol == "Ga"

    def test_cluster(self, xyz_cuo6):
        from pymatgen.core import Molecule
        cl = xyz_cuo6.cluster
        assert isinstance(cl, Molecule) and cl[0].species_string == "Cu"
        assert all(hasattr(s, "label") for s in cl)


# ── absorber_idx setter ──────────────────────────────────────


class TestAbsorberIdx:
    def test_valid(self, cif_zno):
        cif_zno.absorber_idx = 0
        assert cif_zno.absorber_idx == 0

    def test_out_of_range(self, cif_zno):
        with pytest.raises(IndexError, match="not valid"):
            cif_zno.absorber_idx = 999

    def test_wrong_element(self, xyz_cuo6):
        with pytest.raises(IndexError, match="not in site"):
            xyz_cuo6.absorber_idx = 4

    def test_absorber_site(self, xyz_cuo6):
        assert "Cu" in xyz_cuo6.absorber_site.species_string

    def test_second_absorber_site(self, cif_fe3o4):
        cif_fe3o4.absorber_idx = 1
        assert cif_fe3o4.absorber_idx == 16


# ── errors ───────────────────────────────────────────────────

class TestErrors:
    def test_nonexistent_file(self):
        with pytest.raises(FileNotFoundError):
            get_structure("/tmp/_no_such_file_12345.cif", "Zn")

    def test_bad_extension(self):
        with tempfile.NamedTemporaryFile(suffix=".bad", delete=False) as f:
            f.write(b"x"); fname = f.name
        try:
            with pytest.raises(ValueError, match="unknown structure format"):
                get_structure(fname, "Zn")
        finally:
            Path(fname).unlink()


# ── get_structure_from_text ──────────────────────────────────


class TestFromText:
    @pytest.fixture(scope="class")
    def cif_text(cls):
        return (STRUCTS / "ZnO_mp-2133.cif").read_text()

    @pytest.fixture(scope="class")
    def xyz_text(cls):
        return (STRUCTS / "CuO6_D4h.xyz").read_text()

    def test_cif(self, cif_text):
        with pytest.warns(UserWarning, match="fractional coordinates rounded"):
            sg = get_structure_from_text(cif_text, "Zn", format="cif", filename="t.cif")
        assert isinstance(sg, XasStructureCif) and sg.formula == "ZnO"

    def test_xyz(self, xyz_text):
        sg = get_structure_from_text(xyz_text, "Cu", format="xyz", filename="t.xyz")
        assert isinstance(sg, XasStructureXyz) and sg.struct_type == "molecule"

    def test_int_absorber(self, cif_text, xyz_text):
        # CIF triggers pymatgen rounding warning
        with pytest.warns(UserWarning, match="fractional coordinates rounded"):
            result_cif = get_structure_from_text(cif_text, 30, format="cif", filename="t.cif")
        assert result_cif.absorber.symbol == "Zn"
        assert get_structure_from_text(xyz_text, 29, format="xyz", filename="t.xyz").absorber.symbol == "Cu"

    def test_bad_format(self, cif_text):
        with pytest.raises(ValueError, match="unknown structure format"):
            get_structure_from_text(cif_text, "Zn", format="pdb")


# ── get_structs_from_dir ────────────────────────────────────


class TestFromDir:
    def test_single(self):
        s = get_structs_from_dir(STRUCTS, "Fe", globstr="Fe3O4*")
        assert len(s) == 1 and s[0].name == "Fe3O4_cub_fracOcc_Levy2012_ICSD-183969.cif"

    def test_cif_glob(self):
        with pytest.warns(UserWarning, match="fractional coordinates rounded"):
            structs = get_structs_from_dir(STRUCTS, "Zn", globstr="ZnO*")
        assert all(isinstance(x, XasStructureCif) for x in structs)

    def test_xyz_glob(self):
        assert all(isinstance(x, XasStructureXyz)
                   for x in get_structs_from_dir(STRUCTS, "Cu", globstr="CuO6*"))

    def test_exclude(self):
        s = get_structs_from_dir(STRUCTS, "Ga", globstr="GaBr*",
                                 exclude_names=["GaBr_single_frame.xyz"])
        assert any(x.name == "GaBr_multi-frame.xyz" for x in s)
        assert not any(x.name == "GaBr_single_frame.xyz" for x in s)


# ── site_label ──────────────────────────────────────────────


class TestSiteLabel:
    def test_format(self, cif_zno, xyz_cuo6):
        lbl = site_label(cif_zno.struct[0])
        assert "[" in lbl and "]" in lbl and "," in lbl
        assert site_label(xyz_cuo6.struct[0]).startswith("Cu")

    def test_species_included(self, cif_zno):
        for site in cif_zno.struct:
            assert site.species_string.split(":")[0] in site_label(site)