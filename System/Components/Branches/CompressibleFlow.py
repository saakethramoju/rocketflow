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
                 upstream_pressure: State,
                 upstream_temperature: State,
                 downstream_pressure: State,
                 discharge_coefficient: float,
                 cross_sectional_area: float,
                 specific_gas_constant: float,
                 specific_heat_ratio: State,
                 mass_flow: State | None = None):
        
        self.setup()
    

    def evaluate_states(self):

        P1 = self.upstream_pressure.value
        T1 = self.upstream_temperature.value
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
        cross_sectional_area: float,
        length: float,
        inner_diameter: float,
        friction_factor: State | None = None,
        mass_flux: State | None = None,
        upstream_mach_number: State | None = None,
        downstream_mach_number: State | None = None,
        predicted_downstream_pressure: State | None = None,
        flow_regime: str = "None",
    ):
        self.setup()

    def evaluate_states(self):
        mdot = self.mass_flow.value
        P1 = self.upstream_pressure.value
        rho1 = self.upstream_density.value
        a1 = self.upstream_speed_of_sound.value
        k = self.upstream_specific_heat_ratio.value
        A = self.cross_sectional_area.value
        f = self.friction_factor.value
        L = self.length.value
        D = self.inner_diameter.value

        G = mdot / A
        self.mass_flux.value = G

        v1 = G / rho1
        M1 = v1 / a1

        if M1 <= 0.0:
            raise ValueError("FannoFlow currently requires positive inlet Mach number.")

        self.upstream_mach_number.value = M1

        F1 = (1.0 - M1**2) / (k * M1**2) + (k + 1.0) / (2.0 * k) * np.log(((k + 1.0) * M1**2) / (2.0 + (k - 1.0) * M1**2))

        F2 = F1 - f * L / D

        if F2 <= 0.0:
            M2 = 1.0
            self.flow_regime = "choked"
        else:
            branch = -1 if M1 < 1.0 else 0
            W = lambertw(-np.exp(-1.0 - (2.0 * k / (k + 1.0)) * F2), branch).real
            M2 = (-((k + 1.0) / 2.0) * W - ((k - 1.0) / 2.0)) ** -0.5

            if M2 < 1.0:
                self.flow_regime = "subsonic"
            elif M2 > 1.0:
                self.flow_regime = "supersonic"
            else:
                self.flow_regime = "sonic"

        self.downstream_mach_number.value = M2

        pressure_ratio = M1 / M2 * np.sqrt((1.0 + 0.5 * (k - 1.0) * M1**2) / (1.0 + 0.5 * (k - 1.0) * M2**2))
        self.predicted_downstream_pressure.value = P1 * pressure_ratio

    @property
    def iteration_variables(self) -> list[State]:
        return [self.mass_flow]

    @property
    def residuals(self) -> list[float]:
        return [self.predicted_downstream_pressure.value - self.downstream_pressure.value]