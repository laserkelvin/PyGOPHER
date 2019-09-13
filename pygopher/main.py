
import xml.etree.ElementTree as ET
from pathlib import Path


class PgoAttr:
    def __init__(self, name=None, value=None, comment=None):
        self.name = name
        self.value = value
        self.comment = comment
        self._attr_name = self.__class__.__name__

    def __str__(self):
        return f"""<{self._attr_name} Name="{self.name}" Value="{self.value}"/>"""
    
    def __add__(self, other):
        return self.value + other.value
    
    def __sub__(self, other):
        return self.value - other.value
    
    def __div__(self, other):
        return self.value / other.value
    
    def __mul__(self, other):
        return self.value * other.value
    
    def __repr__(self):
        return f"{self.name}, {self.value}" 
    
    
class Parameter(PgoAttr):
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
    def __init__(self, name=None, value=None, comment=None):
        super().__init__(name=name, value=value, comment=comment)


class CartesianTransitionMoment(PgoAttr):
    def __init__(self, axis="a", bra=None, ket=None, parameters=None):
        super().__init__(name="TransitionMoment")
        self.axis = axis
        self.bra = bra
        self.ket = ket
        self.parameters = parameters
        

class Hamiltonian(PgoAttr):
    def __init__(self, name=None, value=None, comment=None, parameters=None,
                 nuclei=None):
        super().__init__(name=name, value=value, comment=comment)
        self.parameters = parameters
        self.nuclei = nuclei
        self.symmetry = "A"
        
    @classmethod
    def from_dict(cls, param_dict, name="v=0", comment=None):
        parameters = [
            Parameter(key, value) for key, value in param_dict.items()
        ]
        object = cls(name=name, comment=comment, parameters=parameters)
        return object
    
    
class AsymmetricTop(Hamiltonian):
    def __init__(self, name=None, value=None, comment=None, parameters=None):
        super().__init__(name=name, value=value, comment=comment, parameters=parameters)
        
        
class SymmetricTop(Hamiltonian):
    def __init__(self, name=None, value=None, comment=None, parameters=None):
        super().__init__(name=name, value=value, comment=comment, parameters=parameters)
        

class LinearTop(Hamiltonian):
    def __init__(self, name=None, value=None, comment=None, parameters=None):
        super().__init__(name=name, value=value, comment=comment, parameters=parameters)
