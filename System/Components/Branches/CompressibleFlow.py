from __future__ import annotations

import numpy as np
from scipy.special import lambertw
from scipy.optimize import brentq
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










class IsentropicAreaChange(Component):
    """
    Ratios are downstream to upstream
    """
    def __init__(
        self,
        name: str,
        network: Network,
        upstream_mach_number: State,
        upstream_static_pressure: State,
        upstream_static_temperature: State,
        gas_constant: float,
        specific_heat_ratio: float,
        upstream_area: float,
        downstream_static_pressure: State | None = None,
        downstream_area: State | float | None = None,
        exit_mach_regime: str = "subsonic",
        downstream_mach_number: State | None = None,
        mass_flow: State | None = None,
        total_enthalpy: State | None = None,
        static_temperature_ratio: State | None = None,
        velocity_ratio: State | None = None,
        density_ratio: State | None = None
    ):
        
        self._use_downstream_pressure = downstream_static_pressure is not None
        self._use_downstream_area = downstream_area is not None

        self.setup()
        temp = self.exit_mach_regime
        self.exit_mach_regime = self.exit_mach_regime.lower()

        if self.exit_mach_regime not in ("subsonic", "supersonic"):
            raise ValueError(
                f"Regime must be 'subsonic' or 'supersonic', got {temp}"
            )

        if self._use_downstream_pressure and self._use_downstream_area:
            self._mode = "pressure"
        elif self._use_downstream_pressure:
            self._mode = "pressure"
        elif self._use_downstream_area:
            self._mode = "area"
        else:
            raise ValueError(
                "IsentropicAreaChange requires either downstream_static_pressure "
                "or downstream_area."
            )

    def evaluate_states(self):
        M1 = self.upstream_mach_number.value
        p1 = self.upstream_static_pressure.value
        T1 = self.upstream_static_temperature.value
        A1 = self.upstream_area.value
        R = self.gas_constant.value
        k = self.specific_heat_ratio.value

        T0 = T1 * (1.0 + (k - 1.0) / 2.0 * M1**2)
        p0 = p1 * (1.0 + (k - 1.0) / 2.0 * M1**2) ** (k / (k - 1.0))

        if self._mode == "pressure":
            p2 = self.downstream_static_pressure.value

            p2_p1 = p2 / p1
            C = (1.0 + (k - 1.0) / 2.0 * M1**2) ** (k / (k - 1.0))

            M2_squared = (2.0 / (k - 1.0)) * (
                (C / p2_p1) ** ((k - 1.0) / k) - 1.0
            )

            M2 = np.sqrt(max(M2_squared, 0.0))

            A1_Astar = self._area_mach_function(M1, k)
            Astar = A1 / A1_Astar
            A2_consistent = Astar * self._area_mach_function(M2, k)

            self.downstream_area.value = A2_consistent

        elif self._mode == "area":
            A2 = self.downstream_area.value

            A1_Astar = self._area_mach_function(M1, k)
            Astar = A1 / A1_Astar

            A2_Astar = A2 / Astar
            M2 = self._inverse_area_mach_function(A2_Astar, k, self.exit_mach_regime)

            p2 = p0 / (
                1.0 + (k - 1.0) / 2.0 * M2**2
            ) ** (k / (k - 1.0))

            self.downstream_static_pressure.value = p2

        else:
            raise ValueError(
                "IsentropicAreaChange requires either downstream_static_pressure "
                "or downstream_area to be assigned."
            )

        mdot = p1 * np.sqrt(k / (R * T1)) * A1 * M1
        h0 = k * R * T0 / (k - 1.0)

        T2_T1 = (1 + (k-1)/2 * M1**2) / (1 + (k-1)/2 * M2**2)
        v2_v1 = (M2/M1) * np.sqrt(T2_T1)
        rho2_rho1 = (p2/p1)**(1/k)

        self.downstream_mach_number.value = M2
        self.mass_flow.value = mdot
        self.total_enthalpy.value = h0
        self.static_temperature_ratio.value = T2_T1
        self.velocity_ratio.value = v2_v1
        self.density_ratio.value = rho2_rho1

    def _area_mach_function(self, M: float, gamma: float) -> float:
        M = float(M)
        gamma = float(gamma)

        if M <= 0.0:
            raise ValueError("Mach number must be positive.")

        return (
            1.0 / M
            * (
                (2.0 / (gamma + 1.0))
                * (1.0 + (gamma - 1.0) / 2.0 * M**2)
            )
            ** ((gamma + 1.0) / (2.0 * (gamma - 1.0)))
        )

    def _inverse_area_mach_function(
        self,
        area_ratio: float,
        gamma: float,
        branch: str = "subsonic",
    ) -> float:
        area_ratio = float(area_ratio)
        gamma = float(gamma)

        if area_ratio < 1.0:
            raise ValueError("A/A* must be >= 1.")

        if np.isclose(area_ratio, 1.0):
            return 1.0

        def residual(M):
            return self._area_mach_function(M, gamma) - area_ratio

        if branch == "subsonic":
            return float(brentq(residual, 1e-12, 1.0 - 1e-12))

        elif branch == "supersonic":
            lo = 1.0 + 1e-12
            hi = 2.0

            while residual(hi) < 0.0:
                hi *= 2.0
                if hi > 1e6:
                    raise RuntimeError("Could not bracket supersonic Mach solution.")

            return float(brentq(residual, lo, hi))

        else:
            raise ValueError("branch must be 'subsonic' or 'supersonic'")
        









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
        regime: str = "subsonic"
    ):
        self.setup()

        self.regime = self.regime.lower()

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






class CompressibleFlowTube(Component):
    """
    Explicit unchoked subsonic Fanno branch.

    Given upstream state, downstream static pressure, friction factor,
    length, and diameter, this component internally solves for the
    upstream Mach number and then computes mass flow and downstream Mach.

    This component does not add global solver residuals.
    """

    def __init__(
        self,
        name: str,
        network: Network,
        mass_flow: State,
        upstream_static_pressure: State,
        upstream_static_temperature: State,
        upstream_density: State,
        downstream_static_pressure: State,
        downstream_static_temperature: State,
        downstream_density: State,
        friction_factor: float,
        length: float,
        inner_diameter: float,
        upstream_static_enthalpy: State | None = None,
        upstream_speed_of_sound: State | None = None,
        downstream_speed_of_sound: State | None = None,
        specific_heat_ratio: State | None = None,

        total_enthalpy: State | None = None,
        upstream_mach_number: State | None = None,
        downstream_mach_number: State | None = None,
        upstream_total_pressure: State | None = None,
        upstream_total_temperature: State | None = None,
        downstream_total_pressure: State | None = None,
        downstream_total_temperature: State | None = None
    ):
        self.setup()

    def evaluate_states(self):
        mdot = self.mass_flow.value
        p1 = self.upstream_static_pressure.value
        rho1 = self.upstream_density.value
        h1 = self.upstream_static_enthalpy.value
        p2 = self.downstream_static_pressure.value
        rho2 = self.downstream_density.value
        f = self.friction_factor.value
        L = self.length.value
        D = self.inner_diameter.value
        A = (np.pi/4) * (D**2)

        u1 = mdot / (rho1 * A)
        u2 = mdot / (rho2 * A)

        if self.upstream_static_enthalpy.is_assigned:
            h1 = self.upstream_static_enthalpy.value
            self.total_enthalpy.value = h1 + 0.5 * u1**2

        if self.upstream_speed_of_sound.is_assigned:
            a1 = self.upstream_speed_of_sound.value
            M1 = u1 / a1
            self.upstream_mach_number.value = M1

        if self.downstream_speed_of_sound.is_assigned:
            a2 = self.downstream_speed_of_sound.value
            M2 = u2 / a2
            self.downstream_mach_number.value = M2

        if (
            self.specific_heat_ratio.is_assigned
            and self.upstream_speed_of_sound.is_assigned
        ):
            k = self.specific_heat_ratio.value
            self.upstream_total_temperature.value = self.upstream_static_temperature.value * (1 + 0.5 * (k - 1) * M1**2)
            self.upstream_total_pressure.value = p1 * (1 + 0.5 * (k - 1) * M1**2)**(k / (k - 1))

        if (
            self.specific_heat_ratio.is_assigned
            and self.downstream_speed_of_sound.is_assigned
        ):
            k = self.specific_heat_ratio.value
            self.downstream_total_temperature.value = self.downstream_static_temperature.value * (1 + 0.5 * (k - 1) * M2**2)
            self.downstream_total_pressure.value = p2 * (1 + 0.5 * (k - 1) * M2**2)**(k / (k - 1))      

        Kf = 8 * f * L / (rho1 * np.pi**2 * D**5)

        inertia = max(mdot, 0.0) * (u2 - u1) - max(-mdot, 0.0) * (u1 - u2)
        pressure = (p1 - p2) * A
        friction = Kf * mdot * np.abs(mdot) * A


        self._residual = pressure - friction - inertia


    @property
    def iteration_variables(self) -> list[State]:
        return [self.mass_flow]
    
    @property
    def residuals(self) -> list[float]:
        return [self._residual]