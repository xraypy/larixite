#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Generating FDMNES input files
==============================

[FDMNES](http://fdmnes.neel.cnrs.fr/) is a program for calculating X-ray
spectroscopy (XAS, XES, RIXS) from the atomic structures

"""

from dataclasses import dataclass
from typing import Union
from pathlib import Path
from pymatgen.core import __version__ as pymatgen_version, Element
from larixite.struct import get_structure
from larixite.struct.xas import XasStructure
from larixite.utils import get_logger, strict_ascii, isotime
from larixite.version import __version__ as larixite_version

logger = get_logger("larixite.fdmnes")

TEMPLATE_FOLDER = Path(Path(__file__).parent, "templates")

FDMNES_DEFAULT_PARAMS = {
    "Energpho": False,
    "Quadrupole": False,
    "Density": False,
    "Density_all": False,
    "SCF": True,
    "Green": True,
    "Memory_save": True,
    "Relativism": False,
    "Spinorbit": None,
    "SCFexc": False,
    "SCFexcv": False,
    "Screening": False,
    "Vmax": False,
    "Full_atom": False,
    "TDDFT": False,
    "PBE96": False,
    "Atom": False,
    "COOP": False,
    "Convolution": True,
}


@dataclass
class FdmnesInput:
    """ "Input parameters for FDMNES"""

    structpath: Union[
        str, Path, XasStructure
    ]  #: path to the structural file or XasStructure
    absorber: Union[
        str, int, Element
    ]  #: atomic symbol or number of the absorbing element
    absorber_idx: Union[int, None] = (
        None  #: index of the absorbing atom in the pymatgen structure
    )
    frame: int = 0  #: index of the frame inside the structure
    edge: Union[str, None] = None  #: edge for calculation
    radius: float = 7
    tmplpath: Union[str, Path, None] = None  #: path to the FDMNES input template
    params: Union[dict, None] = None  #: parameters for FDMNES

    def __post_init__(self):
        """Validate and adjust attributes"""
        if isinstance(self.absorber, str):
            self.absorber = Element(self.absorber)
        elif isinstance(self.absorber, int):
            self.absorber = Element.from_Z(self.absorber)

        if isinstance(self.structpath, str):
            self.structpath = Path(self.structpath)
        if isinstance(self.structpath, XasStructure):
            self.xs = self.structpath
            self.structpath = self.xs.filepath
        else:
            self.xs = get_structure(self.structpath, absorber=self.absorber)

        if self.tmplpath is None:
            self.tmplpath = Path(TEMPLATE_FOLDER, "fdmnes_new.tmpl")
        if isinstance(self.tmplpath, str):
            self.tmplpath = Path(self.tmplpath)

        self.validate_edge()

        if self.params is None:
            self.params = FDMNES_DEFAULT_PARAMS
            self.params = self.get_opt_params()

    def validate_edge(self):
        """Validates and adjusts the edge attribute"""

        valid_edges = [
            "K",
            "L1",
            "L2",
            "L3",
            "L23",
            "M1",
            "M2",
            "M3",
            "M23",
            "M4",
            "M5",
            "M45",
            "N1",
            "N2",
            "N3",
            "N23",
            "N4",
            "N5",
            "N45",
        ]
        if self.edge is None:
            self.edge = "K" if self.absorber.Z < 58 else "L"
        if self.edge == "L":
            self.edge = "L23"
            logger.warning("Edge 'L' changed to 'L23'")
        if self.edge not in valid_edges:
            bad_edge = self.edge
            self.edge = "K" if self.absorber.Z < 58 else "L"
            logger.error(f"Edge {bad_edge} not valid -> changed to {self.edge}")

    def get_opt_params(self) -> dict:
        """Optimize the input parameters"""
        params = self.params.copy()
        atoms_z = [species.Z for species in self.xs.struct.species]
        abs_z = self.absorber.Z
        transition_metals = [range(21, 31), range(39, 49), range(57, 81)]

        if any(abs_z in r for r in transition_metals):
            params["Quadrupole"] = True

        if self.edge == "L23" and abs_z in range(21, 26):
            params["TDDFT"] = True

        if any(z > 36 for z in atoms_z):
            params["Relativism"] = True

        if any(z > 50 for z in atoms_z):
            params["Spinorbit"] = True

        if 8 in atoms_z:
            params["Full_atom"] = True

        return params

    def get_input(self, comment: str = "") -> str:
        params = self.params.copy()
        template = open(self.tmplpath, "r").read()

        comment = (
            f"   {self.xs.name}: {self.absorber.symbol} ({self.absorber.Z}) {self.edge} edge"
            + comment
        )
        #: fill the template
        vers = larixite_version[:]
        if ".post" in vers:
            vers = vers.split(".post")[0]
        conf = {
            "timestamp": isotime(),
            "version": vers,
            "pymatgen_version": pymatgen_version,
            "comment": comment,
            "edge": self.edge,
            "radius": f"{self.radius:.2f}",
        }

        for parkey, parval in params.items():
            conf[parkey] = str(parkey) if parval is True else f"! {parkey}"

        conf["absorber"] = self.absorber.symbol
        conf["absorber_idx"] = self.absorber_idx

        conf["Structure"] = "! TODO"

        return strict_ascii(template.format(**conf))
