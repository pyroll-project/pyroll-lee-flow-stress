from dataclasses import dataclass
from typing import Optional
import numpy as np


@dataclass
class ChemicalComposition:
    """
    Class representing the materials chemical composition in weight percent.
    """
    weight_percent_carbon: Optional[float]
    weight_percent_silicium: Optional[float] = 0
    weight_percent_manganese: Optional[float] = 0
    weight_percent_phosphor: Optional[float] = 0
    weight_percent_vanadium: Optional[float] = 0
    weight_percent_sulfur: Optional[float] = 0
    weight_percent_chromium: Optional[float] = 0
    weight_percent_molybdenum: Optional[float] = 0
    weight_percent_nickel: Optional[float] = 0
    weight_percent_niob: Optional[float] = 0
    weight_percent_copper: Optional[float] = 0
    weight_percent_nitrogen: Optional[float] = 0
    weight_percent_tin: Optional[float] = 0
    weight_percent_lead: Optional[float] = 0

    def __post_init__(self):
        self.iww_equivalent_carbon_content = self.weight_percent_carbon + self.weight_percent_manganese / 6 + (
                self.weight_percent_chromium + self.weight_percent_molybdenum + self.weight_percent_vanadium) / 5 + (
                                                     self.weight_percent_copper + self.weight_percent_nickel) / 15
        self.equivalent_carbon_content = self.iww_equivalent_carbon_content


def flow_stress(chemical_composition: ChemicalComposition, strain: float, strain_rate: float, temperature: float):
    """
    Calculates the flow stress according to the constitutive equation from Lee et al. the provided material composition, strain, strain rate and temperature.

    :param chemical_composition: the chemical composition of the material
    :param strain: the equivalent strain experienced
    :param strain_rate: the equivalent strain rate experienced
    :param temperature: the absolute temperature of the material (K)
    """

    transformation_temperature = 0.95 * (chemical_composition.equivalent_carbon_content + 0.41) / (
            chemical_composition.equivalent_carbon_content + 0.32)

    normalized_temperature = temperature / 1000

    g_fun = 30 * (chemical_composition.equivalent_carbon_content + 0.9) * (
            normalized_temperature - 0.95 * (chemical_composition.equivalent_carbon_content + 0.49) / (
            chemical_composition.equivalent_carbon_content + 0.42)) ** 2 + (
                    chemical_composition.equivalent_carbon_content + 0.06) / (
                    chemical_composition.equivalent_carbon_content + 0.09)

    if normalized_temperature <= transformation_temperature:
        deformation_resistance_contribution = 0.28 * g_fun * np.exp(
            (chemical_composition.equivalent_carbon_content + 0.32) / (
                    0.19 * (chemical_composition.equivalent_carbon_content + 0.41)) - 0.01 / (
                    chemical_composition.equivalent_carbon_content + 0.05))
        strain_rate_sensitivity = (
                                          0.081 * chemical_composition.equivalent_carbon_content - 0.154) * normalized_temperature + (
                                          -0.019 * chemical_composition.equivalent_carbon_content + 0.207) + 0.027 / (
                                          chemical_composition.equivalent_carbon_content + 0.32)

    else:
        deformation_resistance_contribution = 0.28 * np.exp(
            5 / normalized_temperature - 0.01 / (chemical_composition.equivalent_carbon_content + 0.05))
        strain_rate_sensitivity = (
                                          -0.019 * chemical_composition.equivalent_carbon_content + 0.126) * normalized_temperature + (
                                          0.076 * chemical_composition.equivalent_carbon_content - 0.05)

    strain_hardening_factor = 0.41 - 0.07 * chemical_composition.equivalent_carbon_content
    strain_hardening_contribution = 1.3 * (strain / 0.2) ** strain_hardening_factor - 0.3 * (strain / 0.2)

    strain_rate_hardening_contribution = (strain_rate / 10) ** strain_rate_sensitivity * (
            strain_rate / 100) ** (strain_rate_sensitivity / 2.4) * (strain_rate / 1000) ** (
                                                 strain_rate_sensitivity / 15)

    return 9.81 * deformation_resistance_contribution * strain_hardening_contribution * strain_rate_hardening_contribution
