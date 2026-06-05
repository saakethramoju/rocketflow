from __future__ import annotations

from typing import TYPE_CHECKING

from System import Component

if TYPE_CHECKING:
    from System import Network, State


class SiederTate(Component):
    """
    Sieder-Tate turbulent forced-convection heat transfer coefficient.

    Correlation
    -----------
        Nu = 0.027 Re^0.8 Pr^(1/3) (mu / mu_w)^0.14

        Re = mdot * Dh / (mu * A)

        Pr = cp * mu / k

        h = Nu * k / Dh

    Parameters
    ----------
    hydraulic_diameter : State | float
        Hydraulic diameter of the flow passage [m].

    fluid_conductivity : State
        Bulk fluid thermal conductivity [W/m-K].

    fluid_specific_heat : State
        Bulk fluid specific heat capacity [J/kg-K].

    bulk_fluid_dynamic_viscosity : State
        Dynamic viscosity evaluated at the bulk fluid temperature [Pa-s].

    wall_fluid_dynamic_viscosity : State
        Dynamic viscosity evaluated at the wall temperature [Pa-s].

    cross_sectional_area : State | float
        Flow cross-sectional area [m²].

    mass_flow : State
        Fluid mass flow rate [kg/s]. The absolute value is used.

    convection_coefficient : State, optional
        Output convection coefficient h [W/m²-K].
        If omitted, a new State is created.

    Outputs
    -------
    convection_coefficient : State
        Convective heat transfer coefficient [W/m²-K].

    Assumptions
    -----------
    * Single-phase flow.
    * Fully developed turbulent internal flow.
    * Bulk fluid properties are evaluated at the bulk fluid temperature.
    * Wall viscosity is evaluated at the wall temperature.
    * Preferred over Dittus-Boelter when wall and fluid temperatures differ
      enough for viscosity variation to matter.

    Recommended Validity Range
    --------------------------
    * Re >= 10,000
    * 0.7 <= Pr <= 16,700
    * L / Dh > 10

    Notes
    -----
    Sieder-Tate corrects Dittus-Boelter-like turbulent convection for
    temperature-dependent viscosity using (mu / mu_w)^0.14.
    """

    def __init__(
        self,
        name: str,
        network: Network,
        hydraulic_diameter: State | float,
        fluid_conductivity: State,
        fluid_specific_heat: State,
        bulk_fluid_dynamic_viscosity: State,
        wall_fluid_dynamic_viscosity: State,
        cross_sectional_area: State | float,
        mass_flow: State,
        convection_coefficient: State | None = None,
    ):
        self.setup()

    def evaluate_states(self):
        Dh = self.hydraulic_diameter.value
        k = self.fluid_conductivity.value
        Cp = self.fluid_specific_heat.value
        mu = self.bulk_fluid_dynamic_viscosity.value
        mu_w = self.wall_fluid_dynamic_viscosity.value
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
        Nu = 0.027 * Re**0.8 * Pr**(1.0 / 3.0) * (mu / mu_w)**0.14

        self.convection_coefficient.value = Nu * k / Dh



class DittusBoelter(Component):
    """
    Dittus-Boelter (Colburn form) turbulent forced-convection heat transfer
    coefficient for fully developed internal flow.

    Correlation
    -----------
        Nu = 0.023 Re^0.8 Pr^(1/3)

        Re = mdot * Dh / (mu * A)

        Pr = cp * mu / k

        h = Nu * k / Dh

    Parameters
    ----------
    hydraulic_diameter : State | float
        Hydraulic diameter of the flow passage [m].

    fluid_conductivity : State
        Fluid thermal conductivity [W/m-K].

    fluid_specific_heat : State
        Fluid specific heat capacity [J/kg-K].

    fluid_dynamic_viscosity : State
        Fluid dynamic viscosity [Pa-s].

    cross_sectional_area : State | float
        Flow cross-sectional area [m²].

    mass_flow : State
        Fluid mass flow rate [kg/s].

    convection_coefficient : State, optional
        Output convection coefficient h [W/m²-K].
        If omitted, a new State is created.

    Outputs
    -------
    convection_coefficient : State
        Convective heat transfer coefficient [W/m²-K].

    Derived Quantities
    ------------------
    Reynolds number:

        Re = mdot * Dh / (mu * A)

    Prandtl number:

        Pr = cp * mu / k

    Nusselt number:

        Nu = 0.023 Re^0.8 Pr^(1/3)

    Assumptions
    -----------
    * Single-phase flow.
    * Fully developed turbulent internal flow.
    * Fluid properties evaluated at the bulk fluid temperature.
    * Uses the Colburn form with a fixed Prandtl exponent of 1/3.

    Recommended Validity Range
    --------------------------
    * Re >= 10,000
    * 0.7 <= Pr <= 160

    Recommended Temperature Difference Limits
    -----------------------------------------
    For best accuracy, the wall-to-fluid temperature difference should be:

    * |Twall - Tfluid| < 260.928 K (10 °F) for liquids.
    * |Twall - Tfluid| < 310.928 K (100 °F) for gases.

    Notes
    -----
    For large temperature differences, developing flow, rough tubes,
    transitional flow, or strong property variation, more advanced
    correlations such as Sieder-Tate or Gnielinski are generally preferred.
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