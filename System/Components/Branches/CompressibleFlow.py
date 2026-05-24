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

        if not self.total_enthalpy.is_assigned:
            cp = g * R / (g - 1.0)
            self.total_enthalpy.value = cp * To




class ChokedFannoFlow(Component):
    """
    Assumptions
    -----------
    1) Forward flow only
    3) Constant friction factor
    4) Ideal gas
    5) Circular duct

    Supersonic flow tends to be biased towards shorter tubes,
    larger diameters, or smaller friction factors.
    """

    def __init__(
        self,
        name: str,
        network: Network,
        upstream_density: State,
        upstream_speed_of_sound: State,
        specific_heat_ratio: State,
        friction_factor: State,
        length: float,
        inner_diameter: float,

        mass_flux: State | None = None,
        mass_flow: State | None = None,
        upstream_static_enthalpy: State | None = None,
        total_enthalpy: State | None = None,
        upstream_mach_number: State | None = None,
        downstream_mach_number = 1.0,
        regime = "subsonic"
    ):
        self.setup()

    def evaluate_states(self):
        k = self.specific_heat_ratio.value
        rho1 = self.upstream_density.value
        a1 = self.upstream_speed_of_sound.value
        L = self.length.value
        D = self.inner_diameter.value
        f = self.friction_factor.value
        A = (np.pi / 4.0) * D**2

        fL_D = f * L / D

        M1 = self._inverse_fanno_function(fL_D, k, self.regime)

        G = rho1 * M1 * a1
        mdot = G * A

        self.mass_flux.value = G
        self.mass_flow.value = mdot
        self.upstream_mach_number.value = M1

        if self.upstream_static_enthalpy.is_assigned:
            h1 = self.upstream_static_enthalpy.value
            v1 = M1 * a1
            self.total_enthalpy.value = h1 + 0.5*(v1**2)



    def _fanno_function(self, M: float, k: float) -> float:
        return (
            (1.0 - M**2) / (k * M**2)
            + (k + 1.0) / (2.0 * k)
            * np.log(((k + 1.0) * M**2) / (2.0 + (k - 1.0) * M**2))
        )

    def _inverse_fanno_function(
        self,
        target: float,
        k: float,
        branch: str = "subsonic",
    ) -> float:
        if target <= 0.0:
            return 1.0

        c = (k + 1.0) / 2.0
        B = 1.0 + (k / c) * target
        arg = -np.exp(-B)

        if branch == "subsonic":
            W = lambertw(arg, k=-1)
        elif branch == "supersonic":
            W = lambertw(arg, k=0)
        else:
            raise ValueError("branch must be 'subsonic' or 'supersonic'")

        u = -np.real(W)
        x = c * u - (k - 1.0) / 2.0
        M = 1.0 / np.sqrt(x)

        return float(M)




class SubsonicFannoFlow(Component):

    def __init__(self, 
                 name: str, 
                 network: Network):
        self.setup()
