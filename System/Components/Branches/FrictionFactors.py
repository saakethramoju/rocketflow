from __future__ import annotations

import numpy as np
from scipy.special import wrightomega
from typing import TYPE_CHECKING

from System import Component, State

if TYPE_CHECKING:
    from System import Network



class Colebrook(Component):
    """
    Poiseuille number input is only for incompressible flow.
    """

    def __init__(
        self,
        name: str,
        network: Network,
        mass_flow: State,
        friction_factor: State,
        hydraulic_diameter: State | float,
        dynamic_viscosity: State,
        cross_sectional_area: State | float,
        poiseuille_number: float = 16,
        roughness: State | float = 0.0,
        reynolds_number: State | float | None = None,
        reynolds_number_threshold: State | float = 2300.0,
    ):
        self.setup()

        #self.log_friction_factor = State(np.log(self.friction_factor.value))
        self.Deff = 16.0 * self.hydraulic_diameter.value / self.poiseuille_number.value

    def evaluate_states(self):
        mdot = abs(self.mass_flow.value)
        mu = self.dynamic_viscosity.value
        A = self.cross_sectional_area.value
        Dh = self.hydraulic_diameter.value
        Po = self.poiseuille_number.value

        Re_Dh = mdot * Dh / (mu * A)
        Deff = 16.0 * Dh / Po
        Re_eff = mdot * Deff / (mu * A)

        self.Deff = Deff

        if Re_Dh <= self.reynolds_number_threshold.value:
            self.reynolds_number.value = Re_Dh
        else:
            self.reynolds_number.value = Re_eff

        #self.friction_factor.value = np.exp(self.log_friction_factor.value)

    @property
    def iteration_variables(self):
        #return [self.log_friction_factor]
        return [self.friction_factor]

    @property
    def residuals(self):
        f = max(self.friction_factor.value, 1e-12)

        mdot = abs(self.mass_flow.value)
        mu = self.dynamic_viscosity.value
        A = self.cross_sectional_area.value
        Dh = self.hydraulic_diameter.value
        Po = self.poiseuille_number.value
        eps = self.roughness.value

        Re_Dh = max(mdot * Dh / (mu * A), 1e-12)

        if Re_Dh <= self.reynolds_number_threshold.value:
            f_predicted = 4.0 * Po / Re_Dh
            return [self.friction_factor.value - f_predicted]

        Deff = 16.0 * Dh / Po
        Re_eff = max(mdot * Deff / (mu * A), 1e-12)

        return [
            1.0 / np.sqrt(f)
            + 2.0 * np.log10(
                eps / (3.7 * Deff)
                + 2.51 / (Re_eff * np.sqrt(f))
            )
        ]





class Churchill(Component):
    """
    Poiseuille number input is only
    for incompresible flow.
    """
    def __init__(
        self,
        name: str,
        network: Network,
        mass_flow: State,
        friction_factor: State,
        hydraulic_diameter: State | float,
        dynamic_viscosity: State,
        cross_sectional_area: State | float,
        roughness: State | float = 0.0,
        poiseuille_number: float = 16,
        reynolds_number: State | float | None = None,
    ):
        self.setup()

        #self.log_friction_factor = State(np.log(self.friction_factor.value))

        self.Deff = 16*self.hydraulic_diameter.value / self.poiseuille_number.value

    def evaluate_states(self):
        self.Deff = 16 * self.hydraulic_diameter.value / self.poiseuille_number.value
        self.reynolds_number.value = (
            abs(self.mass_flow.value)
            * self.Deff
            / (self.dynamic_viscosity.value * self.cross_sectional_area.value)
        )
        #self.friction_factor.value = np.exp(self.log_friction_factor.value)

    @property
    def iteration_variables(self):
        #return [self.log_friction_factor]
        return [self.friction_factor]

    @property
    def residuals(self):
        Re = max(self.reynolds_number.value, 1e-12)
        relative_roughness = self.roughness.value / self.Deff

        A = (
            2.457
            * np.log(
                1.0 / ((7.0 / Re) ** 0.9 + 0.27 * relative_roughness)
            )
        ) ** 16

        B = (37530.0 / Re) ** 16

        return [self.friction_factor.value -8.0 * (
            (8.0 / Re) ** 12
            + (A + B) ** (-1.5)
        ) ** (1.0 / 12.0)]