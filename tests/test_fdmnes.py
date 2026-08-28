from larixite import get_amcsd
from larixite.utils import get_logger
from larixite.struct import get_structure
from larixite.fdmnes import FdmnesXasInput

logger = get_logger("larixite.test")


def test_fdmnes():
    db = get_amcsd()
    cifids = {
        4438: ("S", "Fe"),
        4820: ("Ti", "Fe"),
        143: ("Fe", "O"),
        2400: ("Fe", "O"),
        2762: ("Fe", "O"),
    }

    for cifid, atoms in cifids.items():
        cif = db.get_cif(cifid)
        outfile = cif.write_cif(verbose=True)
        for abs in atoms:
            logger.info(f"[{cif.label}] {abs}")
            sg = get_structure(outfile, abs)
            f = FdmnesXasInput(sg, absorber=abs)
            text = f.get_input()
            assert len(text) > 700  # TODO: find a better test
            #: test the inputs writes correctly to disk into a temporary directory
            outdir = f.write_input()
            assert outdir.exists()


def test_green_scf_sync():
    db = get_amcsd()
    cif = db.get_cif(4438)
    outfile = cif.write_cif(verbose=True)
    sg = get_structure(outfile, "Fe")

    f = FdmnesXasInput(sg, absorber="Fe", green=False, scf=True, optimize=False)
    assert f.green is False
    assert f.scf is True
    assert f.params["Green"] is False
    assert f.params["SCF"] is True

    f.green = True
    f.scf = False
    assert f.params["Green"] is True
    assert f.params["SCF"] is False


if __name__ == "__main__":
    test_fdmnes()
    test_green_scf_sync()
