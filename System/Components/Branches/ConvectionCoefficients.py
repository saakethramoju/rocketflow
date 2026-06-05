from __future__ import annotations

from typing import TYPE_CHECKING

from System import Component

if TYPE_CHECKING:
    from System import Network, State



class DittusBoelter(Component):
    """
    Dittus-Boelter turbulent forced-convection heat transfer coefficient.

    Uses the Colburn form:

        Nu = 0.023 Re^0.8 Pr^(1/3)
        h  = Nu k / Dh

    Intended for fully developed, single-phase, turbulent internal flow.
    Fluid properties should usually be evaluated at the fluid node temperature.
    """

    def __init__(
        self,
        name: str,
        network: Network,
        hydraulic_diameter: State | float,
        fluid_conductivity: State,
        fluid_specific_heat: State,
        fluid_dynamic_viscosity: State,
        cross_sectional_area: State | float,
        mass_flow: State,
        convection_coefficient: State | None = None,
    ):
        self.setup()

    def evaluate_states(self):
        Dh = self.hydraulic_diameter.value
        k = self.fluid_conductivity.value
        Cp = self.fluid_specific_heat.value
        mu = self.fluid_dynamic_viscosity.value
        A = self.cross_sectional_area.value
        mdot = abs(self.mass_flow.value)

        if Dh <= 0.0:
            raise ValueError(
                f"{self.name}: hydraulic_diameter must be greater than zero. Got {Dh}."
            )

        if A <= 0.0:
            raise ValueError(
                f"{self.name}: cross_sectional_area must be greater than zero. Got {A}."
            )

        Re = mdot * Dh / (mu * A)
        Pr = mu * Cp / k
        Nu = 0.023 * Re**0.8 * Pr**(1.0 / 3.0)

        self.convection_coefficient.value = Nu * k / Dh