from __future__ import annotations

import numpy as np
from scipy.special import wrightomega
from typing import TYPE_CHECKING

from System import Component, State

if TYPE_CHECKING:
    from System import Network




class Colebrook(Component):

    def __init__(
        self,
        name: str,
        network: Network,
        mass_flow: State,
        hydraulic_diameter: State | float,
        dynamic_viscosity: State,
        cross_sectional_area: State | float,
        roughness: State | float = 0.0,
        friction_factor: State | float | None = None,
        reynolds_number: State | float | None = None,
        reynolds_number_threshold: State | float = 2300.0,
    ):
        self.setup()

        self.update_reynolds_number()

        if not self.friction_factor.is_assigned:
            self.friction_factor.value = self.initial_friction_factor()

    def update_reynolds_number(self) -> None:
        if not self.mass_flow.is_assigned:
            self.reynolds_number.value = self.reynolds_number_threshold.value
        else:
            self.reynolds_number.value = (
                abs(self.mass_flow.value)
                * self.hydraulic_diameter.value
                / (self.dynamic_viscosity.value * self.cross_sectional_area.value)
            )

    def laminar_friction_factor(self) -> float:
        Re = max(self.reynolds_number.value, 1e-12)
        return 64.0 / Re

    def turbulent_friction_factor(self) -> float:
        Re = max(self.reynolds_number.value, 1e-12)
        eps = self.roughness.value
        Dh = self.hydraulic_diameter.value

        a = 2.51 / Re
        b = eps / (3.7 * Dh)
        c = np.log(10.0) / 2.0

        z = np.log(c / a) + (c * b) / a
        x = (1.0 / c) * wrightomega(z).real - (b / a)

        return 1.0 / x**2

    def initial_friction_factor(self) -> float:
        if self.reynolds_number.value <= self.reynolds_number_threshold.value:
            return self.laminar_friction_factor()

        return self.turbulent_friction_factor()

    def colebrook_residual(self) -> float:
        f = self.friction_factor.value
        Re = max(self.reynolds_number.value, 1e-12)

        if Re <= self.reynolds_number_threshold.value:
            return f - self.laminar_friction_factor()

        eps = self.roughness.value
        Dh = self.hydraulic_diameter.value

        return (
            1.0 / np.sqrt(f)
            + 2.0 * np.log10(
                eps / (3.7 * Dh)
                + 2.51 / (Re * np.sqrt(f))
            )
        )

    def evaluate_states(self):
        self.update_reynolds_number()

    @property
    def iteration_variables(self):
        return [self.friction_factor]

    @property
    def residuals(self):
        return [self.colebrook_residual()]






class Churchill(Component):

    def __init__(
        self,
        name: str,
        network: Network,
        mass_flow: State,
        hydraulic_diameter: State | float,
        dynamic_viscosity: State,
        cross_sectional_area: State | float,
        roughness: State | float = 0.0,
        friction_factor: State | float | None = None,
        reynolds_number: State | float | None = None,
    ):
        self.setup()

        self.update_reynolds_number()

        if not self.friction_factor.is_assigned:
            self.friction_factor.value = self.churchill_friction_factor()

    def update_reynolds_number(self) -> None:
        self.reynolds_number.value = (
            abs(self.mass_flow.value)
            * self.hydraulic_diameter.value
            / (self.dynamic_viscosity.value * self.cross_sectional_area.value)
        )

    def churchill_friction_factor(self) -> float:
        Re = max(self.reynolds_number.value, 1e-12)
        relative_roughness = self.roughness.value / self.hydraulic_diameter.value

        A = (
            2.457
            * np.log(
                1.0 / ((7.0 / Re) ** 0.9 + 0.27 * relative_roughness)
            )
        ) ** 16

        B = (37530.0 / Re) ** 16

        return 8.0 * (
            (8.0 / Re) ** 12
            + 1.0 / (A + B) ** 1.5
        ) ** (1.0 / 12.0)

    def evaluate_states(self):
        self.update_reynolds_number()

    @property
    def iteration_variables(self):
        return [self.friction_factor]

    @property
    def residuals(self):
        return [
            self.friction_factor.value
            - self.churchill_friction_factor()
        ]