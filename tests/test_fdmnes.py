import pytest

from larixite import get_amcsd
from larixite.struct import get_structure
from larixite.fdmnes import FdmnesXasInput


def test_fdmnes():
    db = get_amcsd()
    cifids = {4438: ("S", "Fe"), 4820: ("Ti", "Fe"), 143: ("Fe", "O")}

    for cifid, atoms in cifids.items():
        cif = db.get_cif(cifid)
        outfile = cif.to_file(verbose=True)
        for abs in atoms:
            sg = get_structure(outfile, abs)
            f = FdmnesXasInput(sg, absorber=abs)
            text = f.get_input()
            assert len(text) > 1000


if __name__ == "__main__":
    test_fdmnes()
