from __future__ import annotations

import numpy as np
from scipy.special import lambertw
from typing import TYPE_CHECKING

from System import Component, State

if TYPE_CHECKING:
    from System import Network




class IsentropicCompressibleOrifice(Component):
    """
    Assumes ideal gas
    """
    def __init__(self,
                 name: str,
                 network: Network,
                 upstream_total_pressure: State,
                 upstream_total_temperature: State,
                 downstream_pressure: State,
                 discharge_coefficient: float,
                 cross_sectional_area: float,
                 specific_gas_constant: float,
                 specific_heat_ratio: State,
                 mass_flow: State | None = None,
                 total_enthalpy: State | None = None):
        
        self.setup()
    

    def evaluate_states(self):

        P1 = self.upstream_total_pressure.value
        T1 = self.upstream_total_temperature.value
        P2 = self.downstream_pressure.value

        CdA = self.discharge_coefficient.value * self.cross_sectional_area.value
        R = self.specific_gas_constant.value
        g = self.specific_heat_ratio.value

        if np.isclose(P1, P2):
            self.mass_flow.value = 0.0
            return

        sign = np.sign(P1 - P2)

        Po = max(P1, P2)
        Pb = min(P1, P2)
        To = T1

        pressure_ratio = Pb / Po
        critical_pressure_ratio = (2 / (g + 1)) ** (g / (g - 1))

        if pressure_ratio <= critical_pressure_ratio:
            flow_function = np.sqrt((g / (R * To)) * (2 / (g + 1)) ** ((g + 1) / (g - 1)))

        else:
            flow_function = np.sqrt((2 * g / (R * To * (g - 1))) * (pressure_ratio ** (2 / g) - pressure_ratio ** ((g + 1) / g)))

        self.mass_flow.value = sign * CdA * Po * flow_function

        cp = g * R / (g - 1.0)
        self.total_enthalpy.value = cp * To



class FannoFlow(Component):
    """
    Assumptions
    -----------
    1) Forward flow only
    2) Constant friction factor
    3) Ideal gas
    4) Circular duct
    """

    def __init__(
        self,
        name: str,
        network: Network,
        mass_flow: State,
        upstream_density: State,
        upstream_speed_of_sound: State,
        upstream_specific_heat_ratio: State,
        downstream_density: State,
        downstream_speed_of_sound: State,
        length: float,
        inner_diameter: float,
        friction_factor: State | None = None,
        mass_flux: State | None = None,
        upstream_mach_number: State | None = None,
        downstream_mach_number: State | None = None,
        upstream_static_enthalpy: State | None = None,
        total_enthalpy: State | None = None
    ):
        self.setup()

    def evaluate_states(self):
        mdot = self.mass_flow.value
        f = self.friction_factor.value
        L = self.length.value
        D = self.inner_diameter.value
        A = (np.pi / 4.0) * D**2

        G = mdot / A
        self.mass_flux.value = G

        if mdot <= 0.0:
            self.upstream_mach_number.value = 0.0
            self.downstream_mach_number.value = 0.0
            self._fL_D_residual = mdot
            return

        rho_in = self.upstream_density.value
        a_in = self.upstream_speed_of_sound.value
        k = self.upstream_specific_heat_ratio.value
        rho_out = self.downstream_density.value
        a_out = self.downstream_speed_of_sound.value

        M1 = G / (rho_in * a_in)
        M2_raw = G / (rho_out * a_out)

        self.upstream_mach_number.value = M1
        self.downstream_mach_number.value = min(M2_raw, 1.0)

        M2 = min(M2_raw, 1.0 - 1e-10)

        fL_D1 = (1.0 - M1**2) / (k * M1**2) + (k + 1.0) / (2.0 * k) * np.log(((k + 1.0) * M1**2) / (2.0 + (k - 1.0) * M1**2))
        fL_D2 = (1.0 - M2**2) / (k * M2**2) + (k + 1.0) / (2.0 * k) * np.log(((k + 1.0) * M2**2) / (2.0 + (k - 1.0) * M2**2))

        self._fL_D_residual = (fL_D1 - fL_D2) - f * L / D

        if self.upstream_static_enthalpy.is_assigned:
            h1 = self.upstream_static_enthalpy.value
            v1 = M1 * a_in
            self.total_enthalpy.value = h1 + 0.5*(v1**2)


    @property
    def iteration_variables(self) -> list[State]:
        return [self.mass_flow]

    @property
    def residuals(self) -> list[float]:
        fL_D_actual = self.friction_factor.value * self.length.value / self.inner_diameter.value
        scale = max(abs(fL_D_actual), 1.0)
        return [self._fL_D_residual / scale]
    

'''

class FannoFlow(Component):
    """
    Assumptions
    -----------
    1) Constant friction factor
    2) Ideal gas
    3) Circular duct
    """

    def __init__(
        self,
        name: str,
        network: Network,
        mass_flow: State,
        upstream_pressure: State,
        upstream_density: State,
        upstream_speed_of_sound: State,
        upstream_specific_heat_ratio: State,
        downstream_pressure: State,
        downstream_density: State,
        downstream_speed_of_sound: State,
        downstream_specific_heat_ratio: State,
        length: float,
        inner_diameter: float,
        friction_factor: State | None = None,
        mass_flux: State | None = None,
        upstream_mach_number: State | None = None,
        downstream_mach_number: State | None = None,
    ):
        self.setup()

    def evaluate_states(self):
        mdot = self.mass_flow.value
        f = self.friction_factor.value
        L = self.length.value
        D = self.inner_diameter.value
        A = (np.pi / 4.0) * D**2

        G = mdot / A
        self.mass_flux.value = G

        if mdot >= 0.0:
            Pin = self.upstream_pressure.value
            rho_in = self.upstream_density.value
            a_in = self.upstream_speed_of_sound.value
            k = self.upstream_specific_heat_ratio.value

            Pout = self.downstream_pressure.value
            rho_out = self.downstream_density.value
            a_out = self.downstream_speed_of_sound.value
        else:
            Pin = self.downstream_pressure.value
            rho_in = self.downstream_density.value
            a_in = self.downstream_speed_of_sound.value
            k = self.downstream_specific_heat_ratio.value

            Pout = self.upstream_pressure.value
            rho_out = self.upstream_density.value
            a_out = self.upstream_speed_of_sound.value



        M1 = np.abs(G) / (rho_in * a_in)
        self.upstream_mach_number.value = M1
        M2_raw = np.abs(G) / (rho_out * a_out)
        M2 = min(M2_raw, 1.0)
        self.downstream_mach_number.value = M2
        M2 = min(M2_raw, 1.0 - 1e-10)

        fL_D1 = (1.0 - M1**2) / (k * M1**2) + (k + 1.0) / (2.0 * k) * np.log(((k + 1.0) * M1**2) / (2.0 + (k - 1.0) * M1**2))
        fL_D2 = (1.0 - M2**2) / (k * M2**2) + (k + 1.0) / (2.0 * k) * np.log(((k + 1.0) * M2**2) / (2.0 + (k - 1.0) * M2**2))

        fL_D_pred = fL_D1 - fL_D2
        fL_D_actual = f*L/D

        self._fL_D_residual = fL_D_pred - fL_D_actual
        
    @property
    def iteration_variables(self) -> list[State]:
        return [self.mass_flow]
        
    @property
    def residuals(self) -> list[float]:
        fL_D_actual = self.friction_factor.value * self.length.value / self.inner_diameter.value
        scale = max(abs(fL_D_actual), 1.0)
        return [self._fL_D_residual / scale]

'''