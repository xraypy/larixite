#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Generating FDMNES input files
==============================

[FDMNES](http://fdmnes.neel.cnrs.fr/) is a program for calculating X-ray
spectroscopy (XAS, XES, RIXS) from the atomic structures

"""

import time
import yaml
from dataclasses import dataclass
from pathlib import Path
import os
import shutil
import subprocess
from pymatgen.core import __version__ as pymatgen_version, Element, Molecule
from larixite.struct import get_structure, get_structure_from_text
from larixite.struct.xas import XasStructure
from larixite.utils import get_logger, strict_ascii, isotime, read_textfile
from larixite.version import __version__ as larixite_version


logger = get_logger("larixite.fdmnes")

TEMPLATE_FOLDER = Path(Path(__file__).parent.parent, "templates")

FDMNES_DEFAULT_PARAMS = {  #: "FDMNES key name": True/False
    "Energpho": False,
    "Quadrupole": False,
    "Polarize": False,  #: -> TODO: implement giving polarization vectors
    "Density": False,
    "Density_all": False,
    "Green": True,
    "Memory_save": True,
    "Relativism": False,
    "Spinorbit": False,
    "SCF": False,
    "SCF_exc": False,
    "R_self": True,
    "N_self": True,
    "P_self": True,
    "Vmax": False,
    "Screening": False,  #: -> TODO
    "Dilatorb": False,  #: -> TODO
    "TDDFT": False,
    "SpGr_atom": False,  #: `Full_atom` is by default since March 2026
    "Hedin": False,  #: `PBE69` is by default since March 2026
    "Atom_conf": False,  #: preferred over `Atom` (permits to keep atomic number in the list of atoms) -> TODO
    "COOP": False,  #: -> TODO
    "Convolution": True,
    "E_cut": False,  #: taken from E_fermi
    "Dec": False,  #: apply energy shift by E_cut
    "Ecent": True,  #: arctan center
    "Elarg": True,  #: arctan width
    "Gamma_hole": False,  #: use core-hole broadening
    "Gamma_max": True,
    "Gaussian": False,
}


def _clean_env():
    """Return environment without SLURM_ variables."""
    env = dict(os.environ)
    for key in list(env):
        if key.startswith("SLURM_"):
            del env[key]
    return env


@dataclass
class FdmnesXasInput:
    """Input generator for a XAS calculation with FDMNES"""

    structpath: (
        str | Path | XasStructure
    )  #: str or path to the structural file, or directly the XasStructure object
    absorber: str | int | Element  #: atomic symbol or number of the absorbing element
    edge: str | None = None  #: edge for calculation
    struct_type: str | None = None  #: type of the structure
    frame: int = (
        0  #: index of the frame inside the structure (e.g. for multi-frame XYZ files)
    )
    radius: float = 7  #: radius of the calculation
    green: bool = True  #: use muffin-tin approximation
    scf: bool = False  #: enable SCF
    erange: str = "-10.0 0.25 60.0 1.0 100.0"  #: energy range
    rself: float | None = None  #: radius for the SCF calculation
    nself: int = 100  #: number of maximum iterations for the SCF calculation
    pself: float = 0.025  #: initial weight for the SCF calculation
    vmax: float | str | None = None  #: maximum potential value for molecules
    fileout_prefix: str | None = None  #: prefix for the FDMNES job (deafault: "job" / extension: .inp)
    tmplpath: str | Path | None = None  #: path to the templates directory
    params: dict[str, bool] | None = None  #: enable/disable parameters for FDMNES
    optimize: bool = False  #: optimize the input parameters
    spacer: str = "   "  #: spacer string for the FDMNES input text
    outdir: str | Path | None = None  #: path to the output directory of the FDMNES jobs
    ecut: float | str | None = None  #: energy cutoff for the convolution
    ecent: float = 30.  #: energy center for the convolution
    elarg: float = 30.  #: energy width for the convolution
    gamma_hole: float | str | None = None  #: start width of the energy convolution (None -> core-hole broadening)
    gamma_max: float = 8  #: maximum energy for the convolution
    gaussian: float | str | None = None  #: gaussian broadening (= experimental braodening)
    estart: float | str | None = None  #: starting energy for the output spectrum

    def __post_init__(self):
        """Validate and optimize attributes"""
        #: absorber
        if isinstance(self.absorber, str):
            self.absorber = Element(self.absorber)
        elif isinstance(self.absorber, int):
            self.absorber = Element.from_Z(self.absorber)
        #: init structure object and type
        if isinstance(self.structpath, str):
            self.structpath = Path(self.structpath)
        if isinstance(self.structpath, XasStructure):
            self._xs = self.structpath
            self.structpath = self._xs.filepath
        else:
            self._xs = get_structure(self.structpath, absorber=self.absorber)
        if self.fileout_prefix is None:
            self.fileout_prefix = "job"
        #: parameters
        if self.params is None:
            self.params = FDMNES_DEFAULT_PARAMS
        self.params["Green"] = self.green
        self.params["SCF"] = self.scf
        #: radius
        self.set_radius(self.radius)
        #: R_self
        if self.rself is None:
            self.set_rself()
        else:
            self.params["R_self"] = True
        #: atictivate Vmax if given by the user
        if self.vmax is None:
            self.params["Vmax"] = False
            self.vmax = "!"
        else:
            self.params["Vmax"] = True
        #: structure type
        if self.struct_type is None:
            self.struct_type = self._xs.struct_type
        else:
            self._xs.struct_type = self.struct_type
        #: template
        if self.tmplpath is None:
            self.tmplpath = Path(TEMPLATE_FOLDER)
        if isinstance(self.tmplpath, str):
            self.tmplpath = Path(self.tmplpath)
        #: absorption edge
        self.validate_edge()
        #: optimize params
        if self.optimize:
            self.params = self.optimize_params()
        #: set outdir to None before writing the input
        if isinstance(self.outdir, str):
            self.outdir = Path(self.outdir)
        if self.outdir is None:
            self.outdir = Path.home() / ".larixite" / "fdmnes"
        #: CONVOLUTION
        if self.ecut is None:
            self.params["E_cut"] = False
            self.ecut = "!"
        else:
            self.params["E_cut"] = True
        if self.gamma_hole is None:
            self.params["Gamma_hole"] = False
            self.gamma_hole = "!"
        else:
            self.params["Gamma_hole"] = True
        if self.gaussian is None:
            self.params["Gaussian"] = False
            self.gaussian = "!"
        else:
            self.params["Gaussian"] = True
        if self.estart is None:
            self.params["Estart"] = False
            self.estart = "!"
        else:
            self.params["Estart"] = True
        #: store a list of jobs
        self._jobs = []

    def __setattr__(self, name: str, value):
        """Sync green/scf with params dict on assignment"""
        super().__setattr__(name, value)
        if name == "green" and hasattr(self, "params") and self.params is not None:
            self.params["Green"] = value
        elif name == "scf" and hasattr(self, "params") and self.params is not None:
            self.params["SCF"] = value

    @property
    def xs(self):
        return self._xs

    def set_radius(self, value: float):
        self.radius = value
        self._xs.radius = value

    def set_rself(self):
        if self.radius < 3.5:
            self.rself = self.radius
        else:
            self.rself = 3.5
        logger.debug(f"R_self set to {self.rself} for radius {self.radius}")

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
        if self.edge == "M":
            self.edge = "M45"
            logger.warning("Edge 'M' changed to 'M45'")
        if self.edge not in valid_edges:
            bad_edge = self.edge
            self.edge = "K" if self.absorber.Z < 58 else "L"
            logger.error(f"Edge {bad_edge} not valid -> changed to {self.edge}")

    def optimize_params(self) -> dict:
        """Optimize the given input parameters"""
        params = self.params.copy()
        atoms_z = [species.Z for species in self._xs.struct.types_of_species]
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
            logger.info(
                "Spinorbit enabled. **NOTE**: the simulations are typically 4 to 8 times longer and need 2 times more memory space"
            )

        if "mol" in self.struct_type.lower():
            params["Vmax"] = True
            self.vmax = -6
            logger.info(f"Molecule: Vmax enablend and set to {self.vmax}")

        #: enable SCF
        self.scf = True
        logger.info("SCF enabled")

        return params

    def get_structure(self, struct_type: str | None = None) -> str:
        """Get the structure section of the input

        Parameters
        ------------

        struct_type: None | str [None -> self.struct_type]
            type of the structure -> see Notes

        Notes
        -----

        FDMNES supports various structure types:
            - Crystal  -> Implemented (default)
            - Molecule -> TODO
            - Cif_file -> Implemented
            - Film  -> Not implemented yet
            - Surface  -> Not implemented yet
            - Interface  -> Not implemented yet
            - Pdb_file  -> Not implemented yet
            - Film_Pdb_file  -> Not implemented yet
            - Film_Cif_file  -> Not implemented yet

        """
        if struct_type is not None:
            self.struct_type = struct_type
            self.optimize_params()
        else:
            struct_type = self.struct_type
        if "crys" in struct_type.lower() and isinstance(self._xs.struct, Molecule):
            errmsg = "cannot generate a crystal input from a molecule -> use `struct_type='molecule'`"
            logger.error(errmsg)
            raise AttributeError(errmsg)
        logger.debug(f"Generating structure section for {struct_type}")
        structout = []
        if "crys" in struct_type.lower():
            #: absorber
            structout.append("Z_absorber")
            structout.append(f"{self.spacer}{self.absorber.Z}")
            #: space group
            spgrp = self._xs.space_group
            if spgrp == 1:  #: FDMNES doesn't recognize 1 as a space group -> `P1`
                spgrp = "P1"
            if spgrp == 2:  #: FDMNES doesn't recognize 2 as a space group -> `P-1`
                spgrp = "P-1"
            structout.append("Spgroup")
            structout.append(f"{self.spacer}{spgrp}")
            #: occupancy
            structout.append("Occupancy")
            #: crystal
            structout.append("Crystal")
            lattice = self._xs.sym_struct.lattice
            structout.append(
                f"{self.spacer}{lattice.a} {lattice.b} {lattice.c} {lattice.alpha} {lattice.beta} {lattice.gamma}"
            )
            for (
                idx,
                site,
                site_index,
                occupancy,
                len_sites,
                wyckoff,
            ) in self._xs.unique_sites:
                zelems = [elem.Z for elem in site.species.elements]
                if not len(set(zelems)) == 1:
                    logger.warning(
                        f"[{self._xs.label}] site {idx} has species with different Z -> {site.species_string}"
                    )
                for elem, elstr in zip(
                    site.species.elements, site.species_string.split(", ")
                ):
                    sitestr = f"{elem.Z:>3d} {site.a:15.10f} {site.b:15.10f} {site.c:15.10f} {occupancy:>5.3f} !{site.label:>4s} {wyckoff:>4s} {elstr:>4s}"
                    structout.append(sitestr)
        elif "cif" in struct_type.lower():
            structout.append("Z_absorber")
            structout.append(f"{self.spacer}{self.absorber.Z}")
            structout.append("Cif_file")
            structout.append(f"{self.spacer}{Path(self._xs.filepath).name}")
        elif "mol" in struct_type.lower():
            #: build the cluster and map by distance from absorber at (0,0,0)
            mol = self._xs.cluster
            map_mol_by_dist = [(0, 0)]
            for i, site in enumerate(mol[1:]):
                isite = i + 1
                dist = mol.get_distance(0, isite)
                map_mol_by_dist.append((isite, dist))
            #: absorber
            structout.append("Absorber")
            structout.append("   1")
            #: molecule
            structout.append("Molecule")
            structout.append(f"{self.spacer}1.0 1.0 1.0 90.0 90.0 90.0")
            for i, dist in sorted(map_mol_by_dist, key=lambda x: x[1]):
                site = mol[i]
                if not len(site.species) == 1:
                    raise NotImplementedError("TODO: partial site occupancy")
                el = site.species.elements[0]
                structout.append(
                    f"{el.Z:>3d} {site.x:15.10f} {site.y:15.10f} {site.z:15.10f} ! {site.label:>6s} {dist:10.5f}"
                )
        else:
            errmsg = f"Structure type `{struct_type}` not supported"
            logger.error(errmsg)
            raise AttributeError(errmsg)
        return "\n".join(structout)

    def get_input(
        self,
        comment: str = "",
        struct_type: str = None,
        template: str | Path | None = None,
    ) -> str:
        """Get the FDMNES input text"""
        params = self.params.copy()
        if template is None:
            template = self.tmplpath / "fdmnes_xas.tmpl"
        if isinstance(template, str):
            template = Path(template)
        template_text = open(template, "r").read()

        comment = (
            f"{self.spacer}{self._xs.name}: {self.absorber.symbol} ({self.absorber.Z}) {self.edge} edge"
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
            "filout": self.spacer + self.fileout_prefix,
            "edge": f"{self.spacer}{self.edge}",
            "radius": f"{self.spacer}{self.radius:.2f}",
            "erange": self.spacer + self.erange,
            "rself": f"{self.spacer}{self.rself}",
            "nself": f"{self.spacer}{self.nself}",
            "pself": f"{self.spacer}{self.pself}",
            "vmax": f"{self.spacer}{self.vmax}",
            "struct_type": self.struct_type,
            "structure": self.get_structure(struct_type=struct_type),
            "ecut": f"{self.spacer}{self.ecut}",
            "ecent": f"{self.spacer}{self.ecent}",
            "elarg": f"{self.spacer}{self.elarg}",
            "gamma_hole": f"{self.spacer}{self.gamma_hole}",
            "gamma_max": f"{self.spacer}{self.gamma_max}",
            "gaussian": f"{self.spacer}{self.gaussian}",
            "estart": f"{self.spacer}{self.estart}",
        }
        for parkey, parval in params.items():
            conf[parkey] = str(parkey) if parval is True else f"! {parkey}"

        return strict_ascii(template_text.format(**conf))

    def write_input(
        self, inputtext: str | None = None, outdir: str | Path | None = None
    ) -> Path:
        """Write the FDMNES input text to disk and return the job output directory"""
        if inputtext is None:
            inputtext = self.get_input()
        if isinstance(outdir, str):
            outdir = Path(outdir)
        if outdir is None:
            outdir = self.outdir
        tstamp = time.strftime("%y%m%d_%H%M%S")  #: 260518_144910
        jobdir = outdir / f"{self.fileout_prefix}_{tstamp}"
        jobdir.mkdir(parents=True, exist_ok=True)
        fileout_name = f"{self.fileout_prefix}.inp"
        fnout = jobdir / fileout_name
        with open(fnout, "w") as fp:
            fp.write(inputtext)
        with open(jobdir / "fdmfile.txt", "w") as fp:
            fp.write(f"1\n{fileout_name}\n")
            logger.info(f"written {fp.name}")
        if self.struct_type is not None and "cif" in self.struct_type.lower():
            src = Path(self._xs.filepath)
            shutil.copy2(src, jobdir)
            logger.info(f"copied {src.name} to {jobdir}")
        self._jobs.append(jobdir)
        _ = self.dump_params(jobdir / f"{self.fileout_prefix}_params.yaml")
        return jobdir

    def write_sbatch(
        self,
        jobdir: str | Path | None = None,
        template: str | Path | None = None,
        ncpus: int = 12,
        nnodes: int = 1,
        mem_per_cpu: str = "16GB",
        constraint: str = "",
        **kwargs,
    ) -> Path:
        """Generates a SBATCH file (SLURM workload manager) using a template

        Arguments
        ---------
        jobdir: str, Path or None
            path to the job directory to write the SBATCH file
            if None: uses last job directory
        template: str, Path or None
            path to the SBATCH template file
            if None, the default template will be used (`fdmnes_sbatch_esrf.tmpl`)
        **kwargs
            keyword arguments to be replaced in the template file

        Returns
        -------
        None: writes `{self.fileout_prefix}.sbatch`

        """
        if template is None:
            template = self.tmplpath / "fdmnes_sbatch_esrf.tmpl"
        if jobdir is None:
            try:
                jobdir = self._jobs[-1]
            except IndexError:
                logger.error("execute `write_input` first")
                return
        if isinstance(jobdir, str):
            jobdir = Path(jobdir)
        sbatchout = jobdir / f"{self.fileout_prefix}.sbatch"
        kwargs = {
            "jobname": self.fileout_prefix,
            "nnodes": nnodes,
            "ncpus": ncpus,
            "mem_per_cpu": mem_per_cpu,
            "constraint": constraint,
            **kwargs,
        }
        with open(sbatchout, "w") as fp, open(template) as tp:
            fp.write(tp.read().format(**kwargs))
            logger.info(f"written {fp.name}")
        return sbatchout

    def run_sbatch(self, jobdir: Path) -> str | None:
        """Submit the sbatch script via `sbatch --parsable`.

        Writes `status.yaml` with status=``submitted`` and the SLURM job ID.
        Returns the job ID on success, None on failure.
        """
        sbatch_script = jobdir / f"{self.fileout_prefix}.sbatch"
        if not sbatch_script.exists():
            logger.warning(f"No sbatch script found for job {self.fileout_prefix}")
            return None

        try:
            result = subprocess.run(
                ["sbatch", "--parsable", str(sbatch_script)],
                capture_output=True,
                text=True,
                check=True,
                cwd=jobdir,
                env=_clean_env(),
            )
            slurm_job_id = result.stdout.strip()
            status_file = jobdir / "status.yaml"
            status_file.write_text(
                yaml.dump(
                    {
                        "status": "submitted",
                        "slurm_job_id": slurm_job_id,
                    }
                )
            )
            logger.info(
                f"SUBMITTED job {self.fileout_prefix} (slurm job id: {slurm_job_id})"
            )
            return slurm_job_id
        except subprocess.CalledProcessError as e:
            logger.error(f"FAILED to submit job {self.fileout_prefix}: {e}")
            return None

    def dump_params(self, yamlpath: str | Path | None = None) -> Path:
        """Dump input parameters to a YAML file

        Parameters
        ----------
        yamlpath : str | Path | None
            path to the output YAML file
            if None: writes `{jobdir}/{fileout_prefix}_params.yaml`
            if jobdir doesn't exist yet, writes next to the structural file

        Returns
        -------
        Path
            path to the written YAML file
        """
        params_dict = {
            "structpath": str(self.structpath),
            "absorber": self.absorber.symbol,
            "edge": self.edge,
            "struct_type": self.struct_type,
            "frame": self.frame,
            "radius": self.radius,
            "erange": self.erange,
            "rself": self.rself,
            "nself": self.nself,
            "pself": self.pself,
            "vmax": self.vmax,
            "fileout_prefix": self.fileout_prefix,
            "params": self.params,
            "optimize": self.optimize,
            "spacer": self.spacer,
            "outdir": str(self.outdir) if self.outdir else None,
        }

        if yamlpath is None:
            if self._jobs:
                yamlpath = self._jobs[-1] / f"{self.fileout_prefix}_params.yaml"
            else:
                yamlpath = (
                    Path(self.structpath).parent / f"{self.fileout_prefix}_params.yaml"
                )

        if isinstance(yamlpath, str):
            yamlpath = Path(yamlpath)

        with open(yamlpath, "w") as fp:
            yaml.dump(params_dict, fp, default_flow_style=False)
        logger.info(f"written {yamlpath.name}")
        return yamlpath

    @classmethod
    def from_yaml(cls, yamlpath: str | Path) -> "FdmnesXasInput":
        """Restore a FdmnesXasInput object from a YAML parameters file

        Parameters
        ----------
        yamlpath : str | Path
            path to the YAML file previously written by `dump_params()`

        Returns
        -------
        FdmnesXasInput
            reconstructed input object
        """
        if isinstance(yamlpath, str):
            yamlpath = Path(yamlpath)

        with open(yamlpath) as fp:
            params_dict = yaml.safe_load(fp)

        kwargs = {k: v for k, v in params_dict.items() if k != "structpath"}
        return cls(params_dict["structpath"], **kwargs)


def struct2fdmnes(
    inp: str | Path,
    absorber: str | int | Element,
    frame: int = 0,
    format: str = "cif",
    filename: None | str = None,
    radius: float = 7,
    edge: str | None = None,
    optimize: bool = True,
    green: bool = True,
) -> dict:
    """convert CIF/XYZ  into a dictionary of {name: text} for FDMNES output files

    Parameters
    ----------
    inp : str or Path
        text of CIF file, name of CIF file, or Path to CIF file
    absorber : str, int, or Element
        Atomic symbo or number of the absorbing element
    frame : int, optional
        Index of the structure for multi-frame structures in the CIF file [0]
    format : str
        format of text : 'cif' or 'xyz' ['cif']
    filename : str
        full path to filename  ['unknown.{format}']
    optimize : bool
        optimize input parameters based on absorber and structure [True]
    green : bool
        use muffin-tin approximation [True]

    Returns
    -------
    dict
        Dictionary with the FDMNES input

    """
    if len(inp) < 512 and Path(inp).exists():
        if filename is None:
            filename = Path(inp).absolute().as_posix()
        inp = read_textfile(inp)
    if filename is None:
        filename = f"unknown.{format}"

    if isinstance(absorber, str):
        absorber = Element(absorber)
    elif isinstance(absorber, int):
        absorber = Element.from_Z(absorber)

    structs = get_structure_from_text(
        inp, absorber, frame=frame, format=format, filename=filename
    )
    fout_prefix = f"{filename.replace('.', '_')}_{absorber.symbol}"
    fout_name = f"{fout_prefix}.inp"
    fdm = FdmnesXasInput(
        structs,
        absorber=absorber,
        frame=frame,
        edge=edge,
        radius=radius,
        fileout_prefix=fout_prefix,
        optimize=optimize,
    )
    fdm.green = green
    return {"fdmfile.txt": f"1\n{fout_name}\n", fout_name: fdm.get_input()}
