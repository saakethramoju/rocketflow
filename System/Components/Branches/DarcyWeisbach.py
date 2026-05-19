from __future__ import annotations

import math
from scipy.special import wrightomega
from typing import TYPE_CHECKING

from System import Component, State

if TYPE_CHECKING:
    from System import Network

class GravityPressureChange(Component):

    def __init__(
        self,
        name: str,
        network: Network,
        upstream_pressure: State | float,
        density: State | float,
        elevation_change: State | float,
        gravitional_acceleration: State | float = 9.80665,
        downstream_pressure: State | None = None,
    ):
        """
        Elevation change is positive upwards
        """
        self.setup()

    def evaluate_states(self) -> None:
        self.downstream_pressure.value = (
            self.upstream_pressure.value
            - self.density.value
            * self.gravitional_acceleration.value
            * self.elevation_change.value
        )




class GenericDarcyWeisbach(Component):

    def __init__(self, 
                 name: str, 
                 network: Network,
                 upstream_pressure: State,
                 downstream_pressure: State,
                 length: float,
                 cross_sectional_area: float,
                 hydraulic_diameter: float,
                 density: State,
                 dynamic_viscosity: State,
                 poiseuille_number: float = 16.0,
                 roughness: float | None = 0.0,
                 Reynolds_number_threshold: float = 2300,
                 mass_flow: State | None = None,
                 friction_factor: State | None = None,
                 flow_regime: str | None = None,
                 Reynolds_number: float | None = None):

        self._friction_factor_provided = friction_factor is not None
        self.setup()

    def evaluate_states(self):
        rho = self.density.value
        L = self.length.value
        P1 = self.upstream_pressure.value
        P2 = self.downstream_pressure.value
        dP = P1 - P2

        A = self.cross_sectional_area.value
        Dh = self.hydraulic_diameter.value
        Po = self.poiseuille_number.value

        friction_factor = (
            self.friction_factor.value
            if self._friction_factor_provided
            else None
        )

        initial_mass_flow = (
            self.mass_flow.value
            if self.mass_flow.is_assigned
            else None
        )

        mdot, f, Re, regime = darcy_weisbach_mass_flow(
            pressure_drop=dP,
            density=rho,
            dynamic_viscosity=self.dynamic_viscosity.value,
            length=L,
            area=A,
            hydraulic_diameter=Dh,
            roughness=self.roughness.value,
            poiseuille_number=Po,
            Reynolds_number_threshold=self.Reynolds_number_threshold.value,
            friction_factor=friction_factor,
            initial_mass_flow=initial_mass_flow,
        )

        self.mass_flow.value = mdot
        self.friction_factor.value = f
        self.Reynolds_number.value = Re
        self.flow_regime = regime


class CircularPipeDarcyWeisbach(Component):
    
    def __init__(self, 
                 name: str, 
                 network: Network,
                 upstream_pressure: State,
                 downstream_pressure: State,
                 length: float,
                 inner_diameter: float,
                 density: State,
                 dynamic_viscosity: State,
                 roughness: float | None = 0.0,
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

        friction_factor = (
            self.friction_factor.value
            if self._friction_factor_provided
            else None
        )

        initial_mass_flow = (
            self.mass_flow.value
            if self.mass_flow.is_assigned
            else None
        )

        mdot, f, Re, regime = darcy_weisbach_mass_flow(
            pressure_drop=dP,
            density=rho,
            dynamic_viscosity=self.dynamic_viscosity.value,
            length=L,
            area=A,
            hydraulic_diameter=Dh,
            roughness=self.roughness.value,
            poiseuille_number=16.0,
            Reynolds_number_threshold=self.Reynolds_number_threshold.value,
            friction_factor=friction_factor,
            initial_mass_flow=initial_mass_flow,
        )

        self.mass_flow.value = mdot
        self.friction_factor.value = f
        self.Reynolds_number.value = Re
        self.flow_regime = regime



class RectangularDuctDarcyWeisbach(Component):

    def __init__(self, 
                 name: str, 
                 network: Network,
                 upstream_pressure: State,
                 downstream_pressure: State,
                 length: float,
                 height: float,
                 width: float,
                 density: State,
                 dynamic_viscosity: State,
                 roughness: float | None = 0.0,
                 Reynolds_number_threshold: float = 2300,
                 mass_flow: State | None = None,
                 friction_factor: State | None = None,
                 flow_regime: str | None = None,
                 Reynolds_number: float | None = None):

        self._friction_factor_provided = friction_factor is not None
        self.setup()

    def evaluate_states(self):
        rho = self.density.value
        L = self.length.value
        P1 = self.upstream_pressure.value
        P2 = self.downstream_pressure.value
        dP = P1 - P2
        height = self.height.value
        width = self.width.value

        A = height * width
        perimeter = 2.0 * (height + width)
        Dh = 4.0 * A / perimeter

        # rectangular duct convention: x = smaller/larger side
        a = max(width, height) / 2.0
        b = min(width, height) / 2.0
        x = b / a

        Po = 23.9201 - 29.436 * x + 30.3872 * x**2 - 10.7128 * x**3

        friction_factor = (
            self.friction_factor.value
            if self._friction_factor_provided
            else None
        )

        initial_mass_flow = (
            self.mass_flow.value
            if self.mass_flow.is_assigned
            else None
        )

        mdot, f, Re, regime = darcy_weisbach_mass_flow(
            pressure_drop=dP,
            density=rho,
            dynamic_viscosity=self.dynamic_viscosity.value,
            length=L,
            area=A,
            hydraulic_diameter=Dh,
            roughness=self.roughness.value,
            poiseuille_number=Po,
            Reynolds_number_threshold=self.Reynolds_number_threshold.value,
            friction_factor=friction_factor,
            initial_mass_flow=initial_mass_flow,
        )

        self.mass_flow.value = mdot
        self.friction_factor.value = f
        self.Reynolds_number.value = Re
        self.flow_regime = regime



class EllipticalDuctDarcyWeisbach(Component):

    def __init__(self, 
                 name: str, 
                 network: Network,
                 upstream_pressure: State,
                 downstream_pressure: State,
                 length: float,
                 semi_major_axis: float,
                 semi_minor_axis: float,
                 density: State,
                 dynamic_viscosity: State,
                 roughness: float | None = 0.0,
                 Reynolds_number_threshold: float = 2300,
                 mass_flow: State | None = None,
                 friction_factor: State | None = None,
                 flow_regime: str | None = None,
                 Reynolds_number: float | None = None):

        self._friction_factor_provided = friction_factor is not None
        self.setup()

    def evaluate_states(self):
        rho = self.density.value
        L = self.length.value
        P1 = self.upstream_pressure.value
        P2 = self.downstream_pressure.value
        dP = P1 - P2

        a = self.semi_major_axis.value
        b = self.semi_minor_axis.value

        # enforce convention: a is larger semi-axis, b is smaller semi-axis
        a, b = max(a, b), min(a, b)

        A = math.pi * a * b

        # Ramanujan ellipse perimeter approximation
        perimeter = math.pi * (
            3.0 * (a + b)
            - math.sqrt((3.0 * a + b) * (a + 3.0 * b))
        )

        Dh = 4.0 * A / perimeter

        x = b / a
        Po = (
            19.7669
            - 4.53458 * x
            - 11.5239 * x**2
            + 22.3709 * x**3
            - 10.0874 * x**4
        )

        friction_factor = (
            self.friction_factor.value
            if self._friction_factor_provided
            else None
        )

        initial_mass_flow = (
            self.mass_flow.value
            if self.mass_flow.is_assigned
            else None
        )

        mdot, f, Re, regime = darcy_weisbach_mass_flow(
            pressure_drop=dP,
            density=rho,
            dynamic_viscosity=self.dynamic_viscosity.value,
            length=L,
            area=A,
            hydraulic_diameter=Dh,
            roughness=self.roughness.value,
            poiseuille_number=Po,
            Reynolds_number_threshold=self.Reynolds_number_threshold.value,
            friction_factor=friction_factor,
            initial_mass_flow=initial_mass_flow,
        )

        self.mass_flow.value = mdot
        self.friction_factor.value = f
        self.Reynolds_number.value = Re
        self.flow_regime = regime






def darcy_weisbach_mass_flow(
    pressure_drop: float,
    density: float,
    dynamic_viscosity: float,
    length: float,
    area: float,
    hydraulic_diameter: float,
    roughness: float = 0.0,
    poiseuille_number: float = 16.0,
    Reynolds_number_threshold: float = 2300.0,
    friction_factor: float | None = None,
    initial_mass_flow: float | None = None,
    max_iter: int = 50,
    rel_tol: float = 1e-10,
    abs_tol: float = 1e-12,
) -> tuple[float, float, float, str]:
    """
    Compute Darcy-Weisbach mass flow for circular or non-circular ducts.

    Returns
    -------
    mass_flow
        Signed mass flow rate.

    friction_factor
        Darcy friction factor.

    Reynolds_number
        Final Reynolds number based on the effective diameter used for
        friction factor evaluation.

    flow_regime
        "laminar", "turbulent", "zero flow", or "<fixed friction factor>".
    """

    if length <= 0.0:
        raise ValueError(f"Darcy-Weisbach length must be positive. Got length={length}.")
        
    if hydraulic_diameter <= 0.0:
        raise ValueError("hydraulic_diameter must be positive.")

    if area <= 0.0:
        raise ValueError("area must be positive.")

    # local aliases for readability
    dP = pressure_drop
    rho = density
    mu = dynamic_viscosity
    L = length
    A = area
    Dh = hydraulic_diameter
    eps = roughness
    Po = poiseuille_number
    Re_thresh = Reynolds_number_threshold

    # no pressure drop -> no flow
    if abs(dP) < 1e-12:
        return 0.0, 0.0, 0.0, "zero flow"

    def friction_factor_from_Re(Re: float) -> tuple[float, str]:

        # prevent divide-by-zero
        Re = max(Re, 1e-12)

        # laminar friction factor using Poiseuille number
        #
        # circular pipe:
        #     Po = 16
        #     f = 64 / Re
        if Re <= Re_thresh:
            return 4.0 * Po / Re, "laminar"

        # effective diameter/Reynolds number method
        # for noncircular turbulent ducts
        Deff = 16.0 * Dh / Po
        Re_eff = Re * Deff / Dh

        a = 2.51 / Re_eff
        b = eps / (3.7 * Deff)
        c = math.log(10.0) / 2.0

        # equivalent to Lambert-W Colebrook solution,
        # but avoids exp overflow
        z = math.log(c / a) + (c * b) / a
        x = (1.0 / c) * wrightomega(z).real - (b / a)

        return 1.0 / x**2, "turbulent"

    # directly use user-specified friction factor
    if friction_factor is not None:

        f = friction_factor

        Kf = 8.0 * f * L / (rho * math.pi**2 * Dh**5)

        mdot = math.copysign(
            math.sqrt(abs(dP) / Kf),
            dP,
        )

        Re = max(abs(mdot) * Dh / (mu * A), 1e-12)

        if Re <= Re_thresh:
            regime = "laminar"
        else:
            regime = "turbulent"

        return mdot, f, Re, regime

    # use previous converged mdot if available;
    # otherwise use inviscid/orifice-like initial guess
    if initial_mass_flow is not None:
        mdot = initial_mass_flow
    else:
        mdot = math.copysign(
            A * math.sqrt(2.0 * rho * abs(dP)),
            dP,
        )

    # fixed-point iteration:
    #
    #     mdot -> Re -> f -> mdot
    #
    # needed because friction factor depends on Reynolds number,
    # which itself depends on mass flow rate
    for _ in range(max_iter):
        
        Re = max(abs(mdot) * Dh / (mu * A), 1e-12)

        f, regime = friction_factor_from_Re(Re)

        # Darcy-Weisbach resistance coefficient
        Kf = 8.0 * f * L / (rho * math.pi**2 * Dh**5)

        # updated mass flow from pressure drop
        mdot_new = math.copysign(
            math.sqrt(abs(dP) / Kf),
            dP,
        )

        # combined relative + absolute convergence check
        if abs(mdot_new - mdot) <= max(
            abs_tol,
            rel_tol * max(abs(mdot_new), 1.0),
        ):
            mdot = mdot_new
            break

        mdot = mdot_new

    # recompute final Re/f so returned values are self-consistent
    Re = max(abs(mdot) * Dh / (mu * A), 1e-12)

    f, regime = friction_factor_from_Re(Re)

    Kf = 8.0 * f * L / (rho * math.pi**2 * Dh**5)

    mdot = math.copysign(
        math.sqrt(abs(dP) / Kf),
        dP,
    )

    return mdot, f, Re, regime