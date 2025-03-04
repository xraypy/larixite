import pytest

from pathlib import Path
from larixite.struct import get_structure
from larixite.utils import get_logger

logger = get_logger("larixite.test")
testdir = Path().cwd()
structsdir = testdir / "structs"

test_structures = (#filename_in_structsdir, absorber_str, absorber_index, occupancy
    ("CuO6_D4h.xyz", "Cu", 0, 1),
    ("VO6_Oh.xyz", "V", 0, 1),
    ("GaBr_multi-frame.xyz", "Ga", 0, 1),
    ("Fe3O4_cub_fracOcc_Levy2012_ICSD-183969.cif", "Fe", 0, 0.5),
    ("Fe2.646Ti.354O4_magnetite_Bosi2009_AMCSD4820.cif", "Fe", 0, 0.823),
    ("Fe2.646Ti.354O4_magnetite_Bosi2009_AMCSD4820.cif", "Fe", 1, 1),
    ("K.87H6.13Fe2.79S2O14_jarosite_Basciano2007_AMCSD4438.cif", "Fe", 3, 1),
)

def test_struct():
    for filename, abs, abs_index, occupancy in test_structures:
        sg = get_structure(structsdir / filename, abs)
        sg.absorber_idx = abs_index
        s = sg.get_site(sg.absorber_idx)
        assert sg.get_occupancy(s) == occupancy, "Wrong occupancy"


if __name__ == "__main__":
    test_struct()
