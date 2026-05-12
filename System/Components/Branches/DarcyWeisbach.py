from __future__ import annotations

import math
from scipy.special import lambertw
from typing import TYPE_CHECKING

from System import Component, State

if TYPE_CHECKING:
    from System import Network


class DarcyWeisbach(Component):
    
    def __init__(self, 
                 name: str, 
                 network: Network,
                 upstream_pressure: State,
                 downstream_pressure: State,
                 length: float,
                 inner_diameter: float,
                 roughness: float | None = 0.0,
                 density: State | None = None,
                 dynamic_viscosity: State | None = None,
                 Reynolds_number_threshold: float = 2300,
                 mass_flow: State | None = None,
                 friction_factor: State | None = None,
                 flow_regime: str | None = None):
        
        self._friction_factor_provided = friction_factor is not None
        self.setup()

    def evaluate_states(self):

        rho = self.density.value
        Dh = self.inner_diameter.value
        l = self.length.value
        P1 = self.upstream_pressure.value
        P2 = self.downstream_pressure.value

        dP = P1 - P2
        A = math.pi * Dh**2 / 4.0

        if abs(dP) < 1e-12:
            self.mass_flow.value = 0.0
            self.flow_regime = "zero flow"
            return

        if self._friction_factor_provided:
            f = self.friction_factor.value
            self.flow_regime = "<fixed friction factor>"

        else:
            mu = self.dynamic_viscosity.value
            Re_thresh = self.Reynolds_number_threshold.value
            eps = self.roughness.value

            if self.mass_flow.is_assigned:
                mdot_for_Re = self.mass_flow.value
            else:
                mdot_for_Re = math.copysign(
                    A * math.sqrt(2.0 * rho * abs(dP)),
                    dP,
                )

            Re = abs(mdot_for_Re) * Dh / (mu * A)
            Re = max(Re, 1e-12)

            if Re <= Re_thresh:
                f = 64.0 / Re
                self.flow_regime = "laminar"
            else:
                a = 2.51 / Re
                b = eps / (3.7 * Dh)
                c = math.log(10.0) / 2.0

                x = (
                    (1.0 / c)
                    * lambertw((c / a) * math.exp((c * b) / a)).real
                    - (b / a)
                )

                f = 1.0 / x**2
                self.flow_regime = "turbulent"

            self.friction_factor.value = f

        print(f)

        Kf = 8.0 * f * l / (rho * math.pi**2 * Dh**5)

        self.mass_flow.value = math.copysign(
            math.sqrt(abs(dP) / Kf),
            dP,
        )