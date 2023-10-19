from __future__ import annotations

import os
import tempfile
from collections import namedtuple
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from lxml import etree

from pygopher import utils

__all__ = ["TransitionMoment", "PGopher", "Simulation", "Molecule"]


"""
Plan for organization

Compartmentalize the problem into three classes:
    - Simulation
        - Temperature, Max J, etc.
    - Molecule
        - Rotational parameters
        - Transition moments
    - PGopher
        - ElementTree builder that puts everything together

"""

TransitionMoment = namedtuple("TransitionMoment", "a b c", defaults=[1.0, 0.0, 0.0])


class PGopher:
    def __init__(self, simulation: Simulation, molecule: Molecule) -> None:
        self.simulation = simulation
        self.molecule = molecule
        # Add molecule as a child of simulation
        self.simulation[-1].append(molecule)

    def to_xml(self) -> str:
        return etree.tostring(self.simulation, encoding="unicode", pretty_print=True)

    def save_xml(self, filepath="simulation.pgo") -> None:
        xml = self.to_xml()
        with open(filepath, "w+") as write_file:
            write_file.write('<?xml version="1.0"?>\n')
            write_file.write(xml)

    def __str__(self) -> str:
        return self.to_xml()

    def simulate(
        self, filepath: str | Path | None = None
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """
        Given a PGopher input file, run the simulation with ``pgo``.

        Parameters
        ----------
        filepath : str | Path | None, optional
            Points to a PGopher input file, by default None, which
            will write the current simulation parameters to a temp
            file, and use that as input.

        Returns
        -------
        tuple[pd.DataFrame, pd.DataFrame]
            _description_
        """
        temp = False
        if filepath is None:
            temp = True
            tempfile.tempdir = os.getcwd()
            with tempfile.NamedTemporaryFile(
                mode="w+", delete=False, suffix=".pgo"
            ) as output:
                xml = self.to_xml()
                output.write('<?xml version="1.0"?>\n')
                output.write(xml)
                filepath = output.name
        # First run for the line list
        process = utils.run_pgopher(filepath)
        linelist_df = utils.parse_linelist(process.stdout)
        # now run for the partition function
        process = utils.run_pgopher(filepath, "--qtable")
        q_df = utils.parse_partition_func(process.stdout)
        if temp is True:
            os.remove(filepath)
        return linelist_df, q_df

    @classmethod
    def from_yml(cls, filepath: str | Path) -> PGopher:
        if isinstance(filepath, str):
            filepath = Path(filepath)
        assert filepath.exists(), f"Specified {filepath} for input, but not found."
        parameters = utils.read_yaml(filepath)
        # First generate the Simulation object
        mixture = parameters.get("mixture", None)
        species = parameters.get("species", None)
        sim_obj = Simulation(mixture, species)
        # Create the Molecule object
        mol_type = parameters.get("mol_type", "Asymmetric")
        mol_param = parameters.get("parameters", None)
        settings = parameters.get("settings", None)
        trans_mom = parameters.get("trans_mom", None)
        mol_obj = Molecule(mol_type, settings, mol_param, trans_mom)
        pgo_obj = cls(sim_obj.to_element(), mol_obj.to_element())
        return pgo_obj

    @classmethod
    def from_rng(
        cls,
        seed: int | None = None,
        constant_max: float = 30000.0,
        constant_min: float = 1000.0,
        distortion: bool = False,
        mixture_kwargs: dict[str, Any] = {},
        species_kwargs: dict[str, Any] = {},
        settings: dict[str, Any] = {},
        parameters: dict[str, Any] = {},
        trans_mom: dict[str, Any] = {},
    ) -> PGopher:
        """
        Generate an asymmetric top from randomly generated constants.

        Parameters
        ----------
        seed : int | None, optional
            Random seed, by default None
        constant_max : float, optional
            Maximum value for rotational constants, by default 30000.0
        constant_min : float, optional
            Minimum value for rotational constants, by default 1000.0
        distortion : bool, optional
            Whether to use centrifugal distortion, by default False
        mixture_kwargs : dict[str, Any], optional
            Kwargs to pass into the ``mixture`` specification, by default {}
        species_kwargs : dict[str, Any], optional
            Kwargs to pass into the ``species`` specification, by default {}
        settings : dict[str, Any], optional
            Kwargs to pass into the ``settings`` of ``Molecule``, by default {}
        parameters : dict[str, Any], optional
            Kwargs to override for spectroscopic parameters in ``Molecule``, by default {}
        trans_mom : dict[str, Any], optional
            Dict mapping for principal axes and transition moments in Debye, by default {}

        Returns
        -------
        _type_
            _description_
        """
        np.random.seed(seed)
        constants = np.sort(np.random.uniform(constant_min, constant_max, 3))
        constants = {key: value for key, value in zip(["C", "B", "A"], constants)}
        # Do a single component spectrum for now
        axis = int(np.random.randint(low=0, high=3))
        dipoles = [0.0] * 3
        dipoles[axis] = 1.0
        dipoles = {key: value for key, value in zip(["a", "b", "c"], dipoles)}
        if distortion:
            cd_terms = np.random.uniform(
                low=0.0, high=[1e-2, 1e-2, 2.0, 1e-2, 1e-2], size=5
            )
            cd_terms = {
                key: value
                for key, value in zip(["DJ", "DJK", "DK", "deltaJ", "deltaK"], cd_terms)
            }
            constants.update(**cd_terms)
        if trans_mom:
            dipoles.update(**trans_mom)
        if parameters:
            constants.update(**parameters)
        mol_obj = Molecule(
            parameters=constants,
            settings=settings,
            trans_mom=dipoles,
            mol_type="Asymmetric",
        )
        sim_obj = Simulation(
            mixture_kwargs=mixture_kwargs, species_kwargs=species_kwargs
        )
        pgo_obj = cls(sim_obj.to_element(), mol_obj.to_element())
        pgo_obj.dipoles = dipoles
        pgo_obj.constants = constants
        return pgo_obj


class Simulation:
    """

    Hierarchy is:
    Mixture -> Species -> Molecule -> mol_type/Transition moments
    """

    def __init__(
        self, mixture_kwargs: dict[str, Any] = {}, species_kwargs: dict[str, Any] = {}
    ):
        # configure default values for mixture
        default_mixture = {
            "Temperature": 300.0,
            "Fmin": 0.001,
            "Fmax": 250000.0,
            "OThreshold": 1e-4,
            "SmallE": 2e-18,
        }
        for key, value in default_mixture.items():
            mixture_kwargs.setdefault(key, value)
        self.mixture = mixture_kwargs
        # configure default values for species
        species_kwargs = {"Name": "Species", "Jmax": 30}
        for key, value in species_kwargs.items():
            species_kwargs.setdefault(key, value)
        self.species = species_kwargs

    def to_element(self) -> etree.Element:
        """
        Convert the Simulation object into an XML Element object.
        This is the main avenue to convert the Python into the PGopher
        readable format.

        Returns
        -------
        Element object

        """
        mixture = etree.Element("Mixture")
        settings = {
            "Units": "MHz",
            "PlotUnits": "MHz",
            "IntensityUnits": "nm2MHzperMolecule",
            "PrintLevel": "CSV",
        }
        for key, value in settings.items():
            mixture.set(key, value)

        for key, value in self.mixture.items():
            parameter = etree.SubElement(mixture, "Parameter")
            parameter.set("Name", key)
            parameter.set("Value", str(value))
        species = etree.Element("Species")
        for key, value in self.species.items():
            species.set(key, str(value))
        mixture.append(species)
        return mixture

    def __str__(self) -> str:
        return etree.dump(self.to_element())


class Molecule:
    def __init__(
        self,
        mol_type: str = "Asymmetric",
        settings: dict[str, Any] = {},
        parameters: dict[str, Any] = {},
        trans_mom: dict[str, Any] = {},
    ):
        """
        Class for handling all the molecule specific parameters, such
        as spectroscopic constants.

        Parameters
        ----------
        mol_type : str, optional
            _description_, by default "Asymmetric"
        settings : dict[str, Any], optional
            _description_, by default {}
        parameters : dict[str, Any], optional
            _description_, by default {}
        trans_mom : dict[str, Any], optional
            _description_, by default {}
        """
        # settings corresponds to the name of the molecule, and corresponds
        # to the AsymmetricMolecule Element
        mol_type = mol_type.capitalize()
        if mol_type not in ["Asymmetric", "Linear", "Symmetric"]:
            raise ValueError(f"mol_type name not recognized - {mol_type}")
        self.mol_type = mol_type
        settings.setdefault("Name", "Molecule")
        self.settings = settings
        default_parameters = {"A": 24023.0, "B": 2102.0, "C": 1962.0}
        for key, value in default_parameters.items():
            parameters.setdefault(key, value)
        self.parameters = parameters
        for key in list(trans_mom.keys()):
            if key not in ["a", "b", "c"]:
                del trans_mom[key]
        # use namedtuple for regular structure
        self.trans_mom = TransitionMoment(**trans_mom)._asdict()

    def to_element(self) -> etree.Element:
        """
        Convert the Python object into an XML representation.
        The resulting Element can be dumped as a string with the
        etree.dump method.

        Returns
        -------
        etree.Element
            Element representation of the Molecule object
        """
        molecule = etree.Element(f"{self.mol_type}Molecule")
        for key, value in self.settings.items():
            molecule.set(key, value)
        manifold = etree.Element(f"{self.mol_type}Manifold")
        for key, value in {"Name": "Ground", "Initial": "True"}.items():
            manifold.set(key, value)
        if self.mol_type != "Linear":
            mol_type = f"{self.mol_type}Top"
        else:
            mol_type = self.mol_type
        # Set up the Hamiltonian parameters
        hamiltonian = etree.Element(f"{mol_type}")
        hamiltonian.set("Name", "v=0")
        if self.mol_type == "Asymmetric":
            hamiltonian.set("Symmetry", "A")
        for key, value in self.parameters.items():
            param = etree.SubElement(hamiltonian, "Parameter")
            param.set("Name", key)
            if type(value) != str:
                value = str(value)
            param.set("Value", value)
        # Build up the transition moments and the corresponding
        # dipole moments
        transition = etree.Element("TransitionMoments")
        for state in ["Bra", "Ket"]:
            transition.set(state, "Ground")
        for axis, value in self.trans_mom.items():
            if value > 0.0:
                moment = etree.SubElement(transition, "CartesianTransitionMoment")
                for state in ["Bra", "Ket"]:
                    moment.set(state, "v=0")
                moment.set("Axis", axis)
                dipole = etree.SubElement(moment, "Parameter")
                dipole.set("Name", "Strength")
                dipole.set("Value", str(value))
        # Set up the hierarchy altogether
        manifold.append(hamiltonian)
        for obj in [manifold, transition]:
            molecule.append(obj)
        return molecule

    def __str__(self) -> str:
        return etree.dump(self.to_element())
