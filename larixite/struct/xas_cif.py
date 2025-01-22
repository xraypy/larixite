#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Specialize XasStructure to handle structures from CIF files
=============================================================
"""

from dataclasses import dataclass
from random import Random
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core import Molecule
from larixite.struct.xas import XasStructure, site_label
from larixite.utils import get_logger

rng = Random()
logger = get_logger("larixite.struct")

@dataclass
class XasStructureCif(XasStructure):

    @property
    def sga(self):
        return SpacegroupAnalyzer(self.struct)

    @property
    def space_group(self):
        return self.sga.get_symmetry_dataset().international

    @property
    def sym_struct(self):
        return self.sga.get_symmetrized_structure()

    @property
    def equivalent_sites(self):
        return self.sym_struct.equivalent_sites


    def build_sites(self):
        """Parse sites of CIF structure to add the following attributes:

        unique_sites:   list of (site[0], wyckoff sym) for unique xtal sites
        unique_map:     mapping of all site_labels to unique_site index
        absorber_sites: list of unique sites with absorber

        """
        # get equivalent sites, mapping of all sites to unique sites,
        # and list of site indexes with absorber
        wyckoff_symbols = self.sym_struct.wyckoff_symbols
        self.unique_sites = []
        self.unique_map = {}
        self.absorber_sites = []
        absname = self.absorber.symbol
        for i, sites in enumerate(self.sym_struct.equivalent_sites):
            self.unique_sites.append((sites[0], len(sites), wyckoff_symbols[i]))
            for site in sites:
                self.unique_map[site_label(site)] = i + 1
            if absname in site.species_string:
                self.absorber_sites.append(i)

        self.atom_sites = {}
        self.atom_site_labels = {}

        for i, dat in enumerate(self.unique_sites):
            site = dat[0]
            label = site_label(site)
            for species in site.species:
                elem = species.name
                if elem in self.atom_sites:
                    self.atom_sites[elem].append(i + 1)
                    self.atom_site_labels[elem].append(label)
                else:
                    self.atom_sites[elem] = [i + 1]
                    self.atom_site_labels[elem] = [label]

        all_sites = {}
        for xat in self.atom_site_labels.keys():
            all_sites[xat] = {}
            for i, label in enumerate(self.atom_site_labels[xat]):
                all_sites[xat][label] = self.atom_sites[xat][i]
        self.all_sites = all_sites

    def build_cluster(self):
        """Build a cluster around the absorber as pymatgen Molecule"""
        csize2 = self.cluster_size**2

        site_atoms = {}  # map xtal site with list of atoms occupying that site
        site_tags = {}

        for i, site in enumerate(self.struct.sites):
            label = site_label(site)
            s_unique = self.unique_map.get(label, 0)
            site_species = [e.symbol for e in site.species]
            #: handle partial occupancy
            if len(site_species) > 1:
                s_els = [s.symbol for s in site.species.keys()]

                s_wts = [s for s in site.species.values()]
                site_atoms[i] = rng.choices(s_els, weights=s_wts, k=1000)
                site_tags[i] = f"({site.species_string:s})_{s_unique:d}"
            else:
                site_atoms[i] = [site_species[0]] * 1000
                site_tags[i] = f"{site.species_string:s}_{s_unique:d}"

        # atom0 = self.struct[a_index]
        atom0 = self.get_site(self.absorber_idx)
        sphere = self.struct.get_neighbors(atom0, self.cluster_size)

        self.symbols = [self.absorber.symbol]
        self.coords = [[0, 0, 0]]
        site0_species = [e.symbol for e in atom0.species]
        if len(site0_species) > 1:
            self.tags = [f"({atom0.species_string})_{self.absorber_idx:d}"]
        else:
            self.tags = [f"{atom0.species_string}_{self.absorber_idx:d}"]

        for i, site_dist in enumerate(sphere):
            s_index = site_dist[0].index
            site_symbol = site_atoms[s_index].pop()

            coords = site_dist[0].coords - atom0.coords
            if (coords[0] ** 2 + coords[1] ** 2 + coords[2] ** 2) < csize2:
                self.tags.append(site_tags[s_index])
                self.symbols.append(site_symbol)
                self.coords.append(coords)

        self.molecule = Molecule(self.symbols, self.coords)
