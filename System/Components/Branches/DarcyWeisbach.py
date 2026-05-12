from __future__ import annotations

import math
from scipy.special import wrightomega
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
                 flow_regime: str | None = None,
                 Reynolds_number: float | None = None):
        
        self._friction_factor_provided = friction_factor is not None
        self.setup()


    def evaluate_states(self):

        rho = self.density.value
        Dh = self.inner_diameter.value
        L = self.length.value
        P1 = self.upstream_pressure.value
        P2 = self.downstream_pressure.value

        dP = P1 - P2
        A = math.pi * Dh**2 / 4.0

        # avoid divide-by-zero / undefined flow direction
        if abs(dP) < 1e-12:
            self.mass_flow.value = 0.0
            self.flow_regime = "zero flow"
            return

        def friction_factor_from_Re(Re: float) -> tuple[float, str]:
            Re = max(Re, 1e-12)

            if Re <= Re_thresh:
                return 64.0 / Re, "laminar"

            a = 2.51 / Re
            b = eps / (3.7 * Dh)
            c = math.log(10.0) / 2.0

            # Equivalent to the Lambert-W form, but avoids exp(...) overflow.
            z = math.log(c / a) + (c * b) / a
            x = (1.0 / c) * wrightomega(z).real - (b / a)

            return 1.0 / x**2, "turbulent"

        # user explicitly provided friction factor
        if self._friction_factor_provided:
            f = self.friction_factor.value

            Kf = 8.0 * f * L / (rho * math.pi**2 * Dh**5)

            self.mass_flow.value = math.copysign(
                math.sqrt(abs(dP) / Kf),
                dP,
            )

            if self.dynamic_viscosity.is_assigned:
                mu = self.dynamic_viscosity.value
                Re = abs(self.mass_flow.value) * Dh / (mu * A)
                self.Reynolds_number.value = max(Re, 1e-12)

                if self.Reynolds_number.value <= self.Reynolds_number_threshold.value:
                    self.flow_regime = "laminar"
                else:
                    self.flow_regime = "turbulent"
            else:
                self.flow_regime = "<fixed friction factor>"

            return

        mu = self.dynamic_viscosity.value
        Re_thresh = self.Reynolds_number_threshold.value
        eps = self.roughness.value

        # use previous converged mdot if available;
        # otherwise use inviscid/orifice-like initial guess
        if self.mass_flow.is_assigned:
            mdot = self.mass_flow.value
        else:
            mdot = math.copysign(A * math.sqrt(2.0 * rho * abs(dP)), dP)

        # Fixed-point iteration is used because Darcy-Weisbach flow is implicit:
        #
        #     mdot -> Re -> friction factor -> mdot
        #
        # A single pass would use an inconsistent friction factor because
        # Re depends on mdot. This loop only resolves the internal pipe
        # constitutive relation; the network solver still handles global
        # mass/pressure balance residuals.

        max_iter = 50
        rel_tol = 1e-10
        abs_tol = 1e-12

        for _ in range(max_iter):
            Re = abs(mdot) * Dh / (mu * A)
            f, self.flow_regime = friction_factor_from_Re(Re)

            Kf = 8.0 * f * L / (rho * math.pi**2 * Dh**5)

            mdot_new = math.copysign(
                math.sqrt(abs(dP) / Kf),
                dP,
            )

            if abs(mdot_new - mdot) <= max(
                abs_tol,
                rel_tol * max(abs(mdot_new), 1.0),
            ):
                mdot = mdot_new
                break

            mdot = mdot_new

        # recompute final Re and f so stored values match final mdot
        Re = abs(mdot) * Dh / (mu * A)
        f, self.flow_regime = friction_factor_from_Re(Re)

        self.Reynolds_number.value = max(Re, 1e-12)
        self.friction_factor.value = f

        # common final Darcy-Weisbach solve
        Kf = 8.0 * f * L / (rho * math.pi**2 * Dh**5)

        self.mass_flow.value = math.copysign(
            math.sqrt(abs(dP) / Kf),
            dP,
        )