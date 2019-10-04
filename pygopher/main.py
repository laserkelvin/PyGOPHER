
from lxml import etree
from pathlib import Path
from warnings import warn
import tempfile
import os
from . import utils

import numpy as np


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

class PGopher:
    def __init__(self, simulation: "Element", molecule: "Element"):
        self.simulation = simulation
        self.molecule = molecule
        # Add molecule as a child of simulation
        self.simulation[-1].append(molecule)
        
    def to_xml(self):
        return etree.tostring(
            self.simulation, encoding='unicode', pretty_print=True
            )
    
    def save_xml(self, filepath="simulation.pgo"):
        xml = self.to_xml()
        with open(filepath, "w+") as write_file:
            write_file.write('<?xml version="1.0"?>\n')
            write_file.write(xml)
    
    def __str__(self):
        return self.to_xml()
    
    def simulate(self, filepath=None):
        temp = False
        if filepath is None:
            temp = True
            tempfile.tempdir = os.getcwd()
            with tempfile.NamedTemporaryFile(mode="w+", delete=False, suffix=".pgo") as output:
                xml = self.to_xml()
                output.write('<?xml version="1.0"?>\n')
                output.write(xml)
                filepath = output.name
        # First run for the line list
        process = utils.run_pgopher(filepath)
        linelist_df = utils.parse_linelist(process.stdout)
        process = utils.run_pgopher(filepath, "--qtable")
        q_df = utils.parse_partition_func(process.stdout)
        if temp is True:
            os.remove(filepath)
        return linelist_df, q_df
    
    @classmethod
    def from_yml(cls, filepath):
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
    def from_rng(cls, seed=None, constant_max=30000., constant_min=1000., distortion=False,
                mixture=None, species=None, settings=None, parameters=None, trans_mom=None):
        """
        Create a PGopher simulation based on a random set of constants.
        This is a specific case of an Asymmetric top.
        
        Parameters
        ----------
        seed : [type], optional
            [description], by default None
        constant_max : [type], optional
            [description], by default 30000.
        constant_min : [type], optional
            [description], by default 1000.
        mixture : [type], optional
            [description], by default None
        species : [type], optional
            [description], by default None
        settings : [type], optional
            [description], by default None
        parameters : [type], optional
            [description], by default None
        trans_mom : [type], optional
            [description], by default None  
        """
        np.random.seed(seed)
        constants = np.sort(np.random.uniform(constant_min, constant_max, 3))
        constants = {key: value for key, value in zip(["C", "B", "A"], constants)}
        # Do a single component spectrum for now
        axis = int(np.random.randint(low=0, high=3))
        dipoles = [0.] * 3
        dipoles[axis] = 1.
        dipoles = {key: value for key, value in zip(["a", "b", "c"], dipoles)}
        if distortion:
            cd_terms = np.random.uniform(
                low=0.,
                high=[1e-2, 1e-2, 2., 1e-2, 1e-2],
                size=5
            )
            cd_terms = {
                key: value for key, value in zip(["DJ", "DJK", "DK", "deltaJ", "deltaK"], cd_terms)
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
            mol_type="Asymmetric"
            )
        sim_obj = Simulation(
            mixture=mixture,
            species=species
        )
        pgo_obj = cls(
            sim_obj.to_element(),
            mol_obj.to_element()
        )
        pgo_obj.dipoles = dipoles
        return pgo_obj
        

class Simulation:
    """
    
    Hierarchy is:
    Mixture -> Species -> Molecule -> mol_type/Transition moments
    """
    def __init__(self, mixture=None, species=None):
        self.mixture = {
            "Temperature": 300.,
            "Fmin": 0.001,
            "Fmax": 250000.,
            "OThreshold": 1e-4,
            "SmallE": 2e-18
        }
        if mixture:
            self.mixture.update(**mixture)
        self.species = {
            "Name": "Species",
            "Jmax": 30
        }
        if species:
            self.species.update(**species)
            
    def to_element(self):
        """
        Convert the Simulation object into an XML Element object.
        This is the main avenue to convert the Python into the PGopher
        readable format.
        
        Returns
        -------
        Element object
            
        """
        mixture = etree.Element(
            "Mixture"
        )
        settings = {
            "Units": "MHz", 
            "PlotUnits": "MHz",
            "IntensityUnits": "nm2MHzperMolecule",
            "PrintLevel": "CSV"
            }
        for key, value in settings.items():
            mixture.set(key, value)
        
        for key, value in self.mixture.items():
            parameter = etree.SubElement(
                mixture,
                "Parameter"
            )
            parameter.set(key, str(value))
        species = etree.Element(
            "Species"
        )
        for key, value in self.species.items():
            species.set(key, str(value))
        mixture.append(species)
        return mixture
    
    def __str__(self):
        return etree.dump(self.to_element())
        
    
class Molecule:
    """
    Class for handling all of the details of the molecule
    """
    def __init__(self, mol_type="Asymmetric",
                settings=None, parameters=None, trans_mom=None):
        # settings corresponds to the name of the molecule, and corresponds
        # to the AsymmetricMolecule Element
        if mol_type not in ["Asymmetric", "Linear", "Symmetric"]:
            warn(f"mol_type name not recognized - {mol_type}")
        self.mol_type = mol_type
        self.settings = {
            "name": "Molecule"
        }
        if settings:
            self.settings.update(**settings)
        self.parameters = {
            "A": 24023.,
            "B": 2102.,
            "C": 1962.
        }
        if parameters:
            self.parameters.update(**parameters)
        self.trans_mom = {
            "a": 1.,
            "b": 0.,
            "c": 0.
        }
        if trans_mom:
            self.trans_mom.update(**trans_mom)
    
    def to_element(self):
        """
        
        Convert the Python object into an XML representation.
        The resulting Element can be dumped as a string with the
        etree.dump method.
        
        Returns
        -------
        [type]
            [description]
        """
        molecule = etree.Element(
            f"{self.mol_type}Molecule"
            )
        for key, value in self.settings.items():
            molecule.set(key, value)
        manifold = etree.Element(
            f"{self.mol_type}Manifold"
        )
        for key, value in {"Name": "Ground", "Initial": "True"}.items():
            manifold.set(key, value)
        if self.mol_type != "Linear":
            mol_type = f"{self.mol_type}Top"
        else:
            mol_type = self.mol_type
        # Set up the Hamiltonian parameters
        hamiltonian = etree.Element(
            f"{mol_type}"
        )
        hamiltonian.set("Name", "v=0")
        for key, value in self.parameters.items():
            param = etree.SubElement(
                hamiltonian,
                "Parameter"
            )
            param.set("Name", key)
            if type(value) != str:
                value = str(value)
            param.set("Value", value)
        # Build up the transition moments and the corresponding
        # dipole moments
        transition = etree.Element(
            "TransitionMoments"
        )
        for state in ["Bra", "Ket"]:
            transition.set(state, "Ground")
        for axis, value in self.trans_mom.items():
            if value > 0.:
                moment = etree.SubElement(
                    transition,
                    "CartesianTransitionMoment"
                )
                for state in ["Bra", "Ket"]:
                    moment.set(state, "v=0")
                moment.set("Axis", axis)
                dipole = etree.SubElement(
                    moment,
                    "Parameter"
                )
                dipole.set("Name", "Strength")
                dipole.set("Value", str(value))
        # Set up the hierarchy altogether
        manifold.append(hamiltonian)
        for obj in [manifold, transition]:
            molecule.append(obj)
        return molecule

    def __str__(self):
        return etree.dump(self.to_element())
