import os
from pathlib import Path
from random import Random
from io import StringIO
from typing import Union
import numpy as np
from pymatgen.core import __version__ as pymatgen_version, Structure, Site, Molecule
from pymatgen.io.cif import CifParser
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from xraydb import atomic_symbol, atomic_number, xray_edge

from .utils import strict_ascii, fcompact, isotime, get_logger
from .amcsd_utils import PMG_CIF_OPTS
from .amcsd import get_cif

from .version import __version__ as x_version

rng = Random()

TEMPLATE_FOLDER = Path(Path(__file__).parent, 'templates')

logger = get_logger("larixite.cif_cluster")
if logger.level != 10:
    import warnings
    warnings.filterwarnings("ignore", category=UserWarning)


def read_cif_structure(ciftext: str) -> Structure:
    """read CIF text, return CIF Structure

    Arguments
    ---------
      ciftext (string): text of CIF file

    Returns
    -------
      pymatgen Structure object
    """
    if CifParser is None:
        raise ValueError("CifParser from pymatgen not available. Try 'pip install pymatgen'.")

    if os.path.exists(ciftext):
        ciftext = open(ciftext, 'r').read()
    try:
        cifstructs = CifParser(StringIO(ciftext), **PMG_CIF_OPTS)

    except Exception:
        raise ValueError('could not parse text of CIF')

    try:
        cstruct = cifstructs.parse_structures(primitive=False)[0]
    except Exception:
        raise ValueError('could not get structure from text of CIF')
    return cstruct


def site_label(site: Site) -> str:
    """
    return a string label for a pymatgen Site object,
    using the species string and fractional coordinates

    Parameters
    ----------
    site : pymatgen Site object

    Returns
    -------
    str
    """
    coords = ','.join([fcompact(s) for s in site.frac_coords])
    return f'{site.species_string}[{coords}]'


class CIF_Cluster():
    """
    CIF structure for generating clusters around a specific crystal site,
    as used for XAS calculations

    """
    def __init__(self, ciftext=None, filename=None, absorber=None,
                 absorber_site=1, with_h=False, cluster_size=10.0):
        self.filename = filename
        self.ciftext = ciftext
        self.set_absorber(absorber)
        self.absorber_site = absorber_site
        self.with_h = with_h
        self.cluster_size = cluster_size
        self.struct = None
        if ciftext is None and filename is not None:
            self.ciftext = open(filename, 'r').read()
        if self.ciftext is not None:
            self.parse_ciftext(self.ciftext)

    def set_absorber(self, absorber=None):
        """
        set the absorbing atom element

        Parameters
        ----------
        absorber : None, int, or str
            if None, no change will be made.
            if int, the atomic number of the absorbing element
            if str, the atomic symbol of the absorbing element

        Notes
        -----
        The absorber atom is assumed to be in the CIF structure.
        """
        self.absorber_z = None
        self.absorber = absorber
        if isinstance(self.absorber, int):
            self.absorber = atomic_symbol(self.absorber)
        if isinstance(self.absorber, str):
            self.absorber_z = atomic_number(self.absorber)

    def parse_ciftext(
        self, ciftext: Union[str, None] = None, absorber: Union[str, int, None] = None
    ):
        """
        re-initialize ciftext and absorber if either argument is not None and
        read the structure
        """
        if absorber is not None:
            self.set_absorber(absorber)
        if ciftext is not None:
            self.ciftext = ciftext
        self.struct = read_cif_structure(self.ciftext)
        self.get_cif_sites()

    def get_cif_sites(self):
        """parse sites of CIF structure to get several components:

           struct.sites:   list of all sites as parsed by pymatgen
           site_labels:    list of site labels
           unique_sites:   list of (site[0], wyckoff sym) for unique xtal sites
           unique_map:     mapping of all site_labels to unique_site index
           absorber_sites: list of unique sites with absorber

        """
        # get equivalent sites, mapping of all sites to unique sites,
        # and list of site indexes with absorber

        self.formula = self.struct.composition.reduced_formula
        sga = SpacegroupAnalyzer(self.struct)
        self.space_group = sga.get_symmetry_dataset().international

        sym_struct = sga.get_symmetrized_structure()
        wyckoff_symbols = sym_struct.wyckoff_symbols

        self.site_labels = []
        for site in self.struct.sites:
            self.site_labels.append(site_label(site))

        self.unique_sites = []
        self.unique_map = {}
        self.absorber_sites = []
        absorber = '~'*30 if self.absorber is None else self.absorber
        for i, sites in enumerate(sym_struct.equivalent_sites):
            self.unique_sites.append((sites[0], len(sites), wyckoff_symbols[i]))
            for site in sites:
                self.unique_map[site_label(site)] = (i+1)
            if absorber in site.species_string:
                self.absorber_sites.append(i)

        self.atom_sites = {}
        self.atom_site_labels = {}

        for i, dat in enumerate(self.unique_sites):
            site = dat[0]
            label = site_label(site)
            for species in site.species:
                elem = species.name
                if elem in self.atom_sites:
                    self.atom_sites[elem].append(i+1)
                    self.atom_site_labels[elem].append(label)
                else:
                    self.atom_sites[elem] = [i+1]
                    self.atom_site_labels[elem] = [label]

        all_sites = {}
        for xat in self.atom_site_labels.keys():
            all_sites[xat] = {}
            for i, label in enumerate(self.atom_site_labels[xat]):
                all_sites[xat][label] = self.atom_sites[xat][i]
        self.all_sites = all_sites

    def build_cluster(self, absorber=None, absorber_site=1, cluster_size=None):
        if absorber is not None:
            self.set_absorber(absorber)
        if cluster_size is None:
            cluster_size = self.cluster_size

        if absorber_site not in self.atom_sites[self.absorber]:
            raise ValueError(f"invalid site for absorber {absorber}: must be in {self.atom_sites[self.absorber]}")

        csize2 = cluster_size**2

        site_atoms = {}  # map xtal site with list of atoms occupying that site
        site_tags = {}

        for i, site in enumerate(self.struct.sites):
            label = site_label(site)
            s_unique = self.unique_map.get(label, 0)
            site_species = [e.symbol for e in site.species]
            if len(site_species) > 1:
                s_els = [s.symbol for s in site.species.keys()]

                s_wts = [s for s in site.species.values()]
                site_atoms[i] = rng.choices(s_els, weights=s_wts, k=1000)
                site_tags[i] = f'({site.species_string:s})_{s_unique:d}'
            else:
                site_atoms[i] = [site_species[0]] * 1000
                site_tags[i] = f'{site.species_string:s}_{s_unique:d}'

        # atom0 = self.struct[a_index]
        atom0 = self.unique_sites[absorber_site-1][0]
        sphere = self.struct.get_neighbors(atom0, self.cluster_size)

        self.symbols = [self.absorber]
        self.coords = [[0, 0, 0]]
        site0_species = [e.symbol for e in atom0.species]
        if len(site0_species) > 1:
            self.tags = [f'({atom0.species_string})_{absorber_site:d}']
        else:
            self.tags = [f'{atom0.species_string}_{absorber_site:d}']

        for i, site_dist in enumerate(sphere):
            s_index = site_dist[0].index
            site_symbol = site_atoms[s_index].pop()

            coords = site_dist[0].coords - atom0.coords
            if (coords[0]**2 + coords[1]**2  + coords[2]**2) < csize2:
                self.tags.append(site_tags[s_index])
                self.symbols.append(site_symbol)
                self.coords.append(coords)

        self.molecule = Molecule(self.symbols, self.coords)


def cif_cluster(ciftext=None, filename=None, absorber=None, cluster_size=10.0):
    "return list of sites for the structure"
    return CIF_Cluster(ciftext=ciftext, filename=filename,
                       absorber=absorber, cluster_size=cluster_size)


def cif_extra_titles(cifid):
    """get 'extra titles' from AMCSD cif id"""
    try:
        cif = get_cif(cifid)
    except Exception:
        return []

    # titles from CIF
    out = [f'Mineral Name: {cif.mineral.name.lower()}',
           f'CIF Source: AmMin Crystal Structure DB, id={cif.ams_id}']

    pub = getattr(cif, 'publication', None)
    if pub is not None:
        authors = ', '.join(pub.authors)
        ptext = f'{pub.journalname} {pub.volume} [{pub.page_first}:{pub.page_last}] ({pub.year})'
        out.append(f'Publication: {ptext}')
        out.append(f'Authors: {authors}')

    celld = [f'a={fcompact(cif.a)}', f'b={fcompact(cif.b)}', f'c={fcompact(cif.c)}']
    cella = [f'alpha={fcompact(cif.alpha)}', f'beta={fcompact(cif.beta)}',
             f'gamma={fcompact(cif.gamma)}']
    out.append(f'Cell Parameter lengths (Ang): {", ".join(celld)}')
    out.append(f'Cell Parameter angles  (deg): {", ".join(cella)}')
    out.append(f'Cell Volume (Ang^3): {fcompact(cif.cell_volume)}')
    out.append(f'Crystal Density (gr/cm^3): {fcompact(cif.crystal_density)}')
    if cif.compound != '<missing>':
        out.append(f'Compound: {cif.compound}')
    return out


def cif2feffinp(ciftext, absorber, template=None, edge=None, cluster_size=8.0,
                absorber_site=None, extra_titles=None, with_h=False,
                version8=True, rng_seed=None, cifid=None):

    """convert CIF text to Feff8 or Feff6l input file

    Arguments
    ---------
      ciftext (string):         text of CIF file or name of the CIF file.
      absorber (string or int): atomic symbol or atomic number of absorbing element
                                (see Note 1)
      edge (string or None):    edge for calculation (see Note 2)     [None]
      cluster_size (float):     size of cluster, in Angstroms         [8.0]
      absorber_site (None or int):  index of site for absorber (see Note 3) [None]
      extra_titles (list of str or None): extra title lines to include [None]
      with_h (bool):            whether to include H atoms [False]
      version8 (bool):          whether to write Feff8l input (see Note 5)[True]
      rng_seed (int or None):   seed for RNG to get reproducible occupancy selections [None]
    Returns
    -------
      text of Feff input file

    Notes
    -----
      1. absorber is the atomic symbol or number of the absorbing element, and
         must be an element in the CIF structure. If absorber_site is None (default),
         the first site for that element will be used, as found from cluster.atom_sites.

      2. If edge is a string, it must be one of 'K', 'L', 'M', or 'N' edges (note
         Feff6 supports only 'K', 'L3', 'L2', and 'L1' edges). If edge is None,
         it will be assigned to be 'K' for absorbers with Z < 58 (Ce, with an
         edge energy < 40 keV), and 'L3' for absorbers with Z >= 58.
      3. for structures with multiple sites for the absorbing atom, the site
         can be selected by the order in which they are listed in the sites
         list. This depends on the details of the CIF structure, which can be
         found with `cif_sites(ciftext)`, starting counting by 1.
      5. if version8 is False, outputs will be written for Feff6l

    """
    global rng
    if rng_seed is not None:
        rng.seed(rng_seed)

    if template is None:
        template = open(Path(TEMPLATE_FOLDER, 'feff_exafs.tmpl'), 'r').read()

    cluster = CIF_Cluster(ciftext=ciftext, absorber=absorber, cluster_size=cluster_size+0.5)

    if absorber_site is None:
        absorber_site = cluster.atom_sites[absorber][0]
    cluster.build_cluster(absorber_site=absorber_site, cluster_size=cluster_size)

    mol = cluster.molecule

    absorber = cluster.absorber
    absorber_z = cluster.absorber_z
    if edge is None:
        edge = 'K' if absorber_z < 58 else 'L3'

    edge_energy = xray_edge(absorber, edge).energy
    edge_comment = f'{absorber:s} {edge:s} edge, around {edge_energy:.0f} eV'

    unique_pot_atoms = []
    for site in cluster.struct:
        for elem in site.species.elements:
            if elem.symbol not in unique_pot_atoms:
                unique_pot_atoms.append(elem.symbol)

    atoms_map = {}
    for i, atom in enumerate(unique_pot_atoms):
        atoms_map[atom] = i + 1

    if absorber not in atoms_map:
        atlist = ', '.join(atoms_map.keys())
        raise ValueError(f'atomic symbol {absorber:s} not listed in CIF data: ({atlist})')

    # titles
    titles = [f'Formula: {cluster.formula:s}',
              f'SpaceGroup: {cluster.space_group:s}']

    if extra_titles is not None:
        titles.extend(extra_titles)
    if cifid is not None:
        titles.extend(cif_extra_titles(int(cifid)))


    # comments
    comments = ['*', '* crystallographic sites:',
                '*    To change the absorber site, re-run using `absorber_site`',
                '*    with the corresponding site index (counting from 1)',
                '* site    X       Y       Z     Wyckoff  species']

    for i, dat in enumerate(cluster.unique_sites):
        site, n, wsym = dat
        fc = site.frac_coords
        species_string = site.species_string.strip()
        marker = '  <- absorber' if  ((i+1) == absorber_site) else ''
        s1 = f'{i+1:3d}   {fc[0]:.5f} {fc[1]:.5f} {fc[2]:.5f}'
        s2 = f'{wsym:>5s}   {species_string:s} {marker:s}'
        comments.append(f'* {s1}  {s2}')
    comments.append('*')


    # loop to find atoms actually in cluster, in case some atom
    # (maybe fractional occupation) is not included
    ipot_lines = []
    ipot_map = {}
    next_ipot = 0
    at_lines = [(0, mol[0].x, mol[0].y, mol[0].z, 0, absorber, cluster.tags[0])]
    for i, site in enumerate(mol[1:]):
        sym = site.species_string
        if sym == 'H' and not with_h:
            continue
        if sym in ipot_map:
            ipot = ipot_map[sym]
        else:
            next_ipot += 1
            ipot_map[sym] = ipot = next_ipot

        dist = mol.get_distance(0, i+1)
        at_lines.append((dist, site.x, site.y, site.z, ipot, sym, cluster.tags[i+1]))

    if len(ipot_map) > 10:
        comments.append('*** WARNING: Feff 8l is limited to 11 unique potentials***')
        comments.append('*** WARNING: This input file may need editing')

    # ipots
    ipot, z = 0, absorber_z
    ipot_lines = [f'  {ipot:4d}  {z:>4d}   {absorber:>3s}']
    for sym, ipot in ipot_map.items():
        z = atomic_number(sym)
        ipot_lines.append(f'  {ipot:4d}  {z:>4d}   {sym:>3s}')

    # ordered atoms list
    acount = 0
    atoms = []
    for dist, x, y, z, ipot, sym, tag in sorted(at_lines, key=lambda x: x[0]):
        acount += 1
        sym = (sym + ' ')[:2]
        xyzi = f'  {x:+.5f}  {y:+.5f}  {z:+.5f} {ipot:2d}'.replace(' +', '  ')
        atoms.append(f'{xyzi}  {sym:>3s}  {dist:.5f}  * {tag:s}')


    # now ready to write with template
    vers = x_version[:]
    if '.post' in vers:
        vers = vers.split('.post')[0]
    conf = {'version': vers, 'timestamp': isotime(),
            'pymatgen_version': pymatgen_version, 'edge': edge,
            'radius': f'{cluster_size:.2f}' }

    conf['titles'] = '\n'.join([f'TITLE {x}' for x in titles])
    conf['comments'] = '\n'.join(comments)
    conf['potentials'] = '\n'.join(ipot_lines)
    conf['atoms'] = '\n'.join(atoms)

    return strict_ascii(template.format(**conf))
