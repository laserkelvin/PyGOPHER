
from lxml import etree, builder
from pathlib import Path


class PgopherSimulation:
    """
    
    Hierarchy is:
    Mixture -> Species -> Molecule -> Hamiltonian/Transition moments
    """
    def __init__(self, T=300., **kwargs):
        self.__dict__.update(**kwargs)
        self.mixture = PgopherObject(
            **{"Units": "MHz", "PlotUnits": "MHz"}
        )
        self.species = PgopherObject(
            **{"Name": "Species", "Jmax": "20"}
        )
        self.molecule = AsymmetricTop(
            name="AsymmetricTop",
        )
        


class PgopherObject:
    def __init__(self, name=None, value=None, **kwargs):
        self.name = name
        self.value = value
        self.__dict__.update(**kwargs)
    
    def __add__(self, other):
        return float(self.value) + float(other.value)
    
    def __sub__(self, other):
        return float(self.value) - float(other.value)
    
    def __div__(self, other):
        return float(self.value) / float(other.value)
    
    def __mul__(self, other):
        return float(self.value) * float(other.value)

    def __repr__(self):
        return f"{self.name}, {self.value}" 
    
    def to_element(self):
        element = etree.Element(self.__class__.__name__)
        for key, value in self.__dict__.items():
            if value is not None:
                element.set(key, str(value))
        return element
    
    @classmethod
    def from_element(cls, element):
        obj = cls(name=element.tag)
        obj.__dict__.update(**element)
        return obj
        
    
class Parameter(PgopherObject):
    """
    Core building block of other objects; used to define parameter
    and its associated value.
    
    Parameters
    ----------
    name : str
        Name of the parameter
    value : float
        Value of the parameter
    comment : str
        A comment line for the parameter, e.g. source/origin
    """
    def __init__(self, name=None, value=None, comment=None, **kwargs):
        super().__init__(name=name, value=value, **kwargs)
        self.comment = comment
        

class TransitionMoments(PgopherObject):
    def __init__(self, bra="Ground", ket="Ground", **kwargs):
        super().__init__(**kwargs)
        self.bra = bra
        self.ket = ket


class CartesianTransitionMoment(TransitionMoments):
    """
    Object counterpart to the CartesianTransitionMoment element, which
    contains information about the dipole moment and its projection
    
    Parameters
    ----------
    TransitionMoments : [type]
        [description]
    """
    def __init__(self, axis="a", bra="v=0", ket="v=0", dipole=1.):
        super().__init__(name="TransitionMoment")
        assert axis in ["a", "b", "c"]
        self.axis = axis
        self.bra = bra
        self.ket = ket
        self.dipole = dipole
        
    def to_element(self):
        element = etree.Element(self.__class__.__name__)
        for key, value in self.__dict__.items():
            if key != "dipole":
                element.set(key, str(value))
        element.append(
            Parameter(name="Strength", value=dipole)
        )
        return element
    
    @classmethod
    def from_element(cls, element):
        obj = cls(name=element.tag)
        for key, value in element.items():
            obj.__dict__.set(
                key.lower(), value
            )
        for child in element.iter():
            if child.tag == "Parameter":
                obj.dipole = float(child.get("Strength"))
        return obj
        

class Hamiltonian(PgopherObject):
    def __init__(self, name=None, value=None, comment=None, parameters=None,
                 nuclei=None, **kwargs):
        super().__init__(name=name, value=value, comment=comment)
        self.nuclei = nuclei
        self.symmetry = "A"
        self.__dict__.update(**kwargs)
   
   def to_element(self):
        element = etree.Element(self.__class__.__name__)
        try:
            for key, value in self.parameters.items():
                if value is not None:
                    element.set(key, str(value))
            return element
        except AttributeError:
            raise Exception(f"{self.__class__.__name__} contains no parameters.")
        
    
class AsymmetricTop(Hamiltonian):
    def __init__(self, name=None, value=None, comment=None, parameters=None, **kwargs):
        super().__init__(name=name, value=value, comment=comment,)
        self.parameters = {
            "A": "27000",
            "B": "3600",
            "C": "2560"
        }
        if parameters:
            self.parameters.update(**parameters)
        
        
class SymmetricTop(Hamiltonian):
    def __init__(self, name=None, value=None, comment=None, parameters=None):
        super().__init__(name=name, value=value, comment=comment, parameters=parameters)
        

class LinearTop(Hamiltonian):
    def __init__(self, name=None, value=None, comment=None, parameters=None):
        super().__init__(name=name, value=value, comment=comment, parameters=parameters)
