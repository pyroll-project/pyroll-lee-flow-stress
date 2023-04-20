from pyroll.core import Profile, DeformationUnit, Hook

from .lee_flow_stress import ChemicalComposition, flow_stress

VERSION = "2.0"

Profile.chemical_composition = Hook[ChemicalComposition]()


@DeformationUnit.Profile.flow_stress
def lee_flow_stress(self: DeformationUnit.Profile):
    if hasattr(self, "chemical_composition"):
        return flow_stress(
            self.chemical_composition,
            self.strain,
            self.unit.strain_rate,
            self.temperature
        )


@DeformationUnit.Profile.flow_stress_function
def lee_flow_stress_function(self: DeformationUnit.Profile):
    if hasattr(self, "chemical_composition"):
        def f(strain: float, strain_rate: float, temperature: float) -> float:
            return flow_stress(self.chemical_composition, strain, strain_rate, temperature)

        return f


from . import materials
