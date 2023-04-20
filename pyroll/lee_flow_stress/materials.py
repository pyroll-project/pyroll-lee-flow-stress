from pyroll.core import Profile

from .lee_flow_stress import ChemicalComposition


def is_material(profile: Profile, materials: set[str]):
    if isinstance(profile.material, str):
        return profile.material.lower() in materials
    return materials.intersection([m.lower() for m in profile.material])


@Profile.chemical_composition
def c20(self: Profile):
    if is_material(self, {"c20", "c22"}):
        return ChemicalComposition(
            weight_percent_carbon=0.2,
            weight_percent_silicium=0.2,
            weight_percent_manganese=0.5,
            weight_percent_phosphor=0.02,
            weight_percent_sulfur=0.017,
            weight_percent_copper=0.2
        )

    @Profile.chemical_composition
    def c45(self: Profile):
        if is_material(self, {"c45"}):
            return ChemicalComposition(
                weight_percent_carbon=0.45,
                weight_percent_silicium=0.27,
                weight_percent_manganese=0.61,
                weight_percent_phosphor=0.02,
                weight_percent_chromium=0.32,
                weight_percent_molybdenum=0.05
            )

    @Profile.chemical_composition
    def s355j2(self: Profile):
        if is_material(self, {"s355j2"}):
            return ChemicalComposition(
                weight_percent_carbon=0.21,
                weight_percent_silicium=0.47,
                weight_percent_manganese=1.23,
                weight_percent_phosphor=0.035,
                weight_percent_sulfur=0.02
            )
