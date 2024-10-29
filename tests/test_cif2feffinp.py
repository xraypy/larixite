import pytest

from larixite import get_amcsd, cif2feffinp

def test_cif2feff_v1():

    db = get_amcsd()
    cifids = {4438: ('S', 'Fe'),
              4820: ('Ti', 'Fe'),
              143: ('Fe', 'O')}

    for cifid, atoms in cifids.items():
        cif = db.get_cif(cifid)
        for absorber in atoms:
            for with_h in (True, False):
                text = cif2feffinp(cif.ciftext,
                                   absorber=absorber,
                                   cifid=cifid, cluster_size=6,
                                   with_h=with_h)
                assert len(text) > 2000

if __name__ == '__main__':
    test_cif2feff_v1()
