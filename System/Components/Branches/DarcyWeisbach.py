from __future__ import annotations

import math
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
        mass_flow: State | None = None,
        effective_area: State | None = None,
    ):
        """
        Computes pressure change from elevation.

        elevation_change is positive upward.

        If mass_flow is assigned, this also computes an equivalent
        effective area using:

            mdot = CdA * sqrt(2 * rho * abs(dP))
        """
        self.setup()

        self.evaluate_states()

    def evaluate_states(self) -> None:
        self.downstream_pressure.value = (
            self.upstream_pressure.value
            - self.density.value
            * self.gravitional_acceleration.value
            * self.elevation_change.value
        )

        if self.mass_flow.is_assigned:
            dP = self.upstream_pressure.value - self.downstream_pressure.value
            self.effective_area.value = effective_area_from_mass_flow(
                self.mass_flow.value,
                dP,
                self.density.value,
            )


class GenericDarcyWeisbach(Component):

    def __init__(
        self,
        name: str,
        network: Network,
        upstream_pressure: State,
        downstream_pressure: State,
        length: float,
        cross_sectional_area: float,
        hydraulic_diameter: float,
        density: State,
        friction_factor: State | float | None = None,
        mass_flow: State | float | None = None,
        effective_area: State | None = None,
    ):
        self.setup()
        self._predicted_mass_flow = None

        if not self.mass_flow.is_assigned:
            self.mass_flow.value = initial_mass_flow_guess(
                self.upstream_pressure.value - self.downstream_pressure.value,
                self.density.value,
                self.cross_sectional_area.value,
            )

    def evaluate_states(self):
        self._predicted_mass_flow = predicted_darcy_mass_flow(
            pressure_drop=self.upstream_pressure.value - self.downstream_pressure.value,
            density=self.density.value,
            length=self.length.value,
            hydraulic_diameter=self.hydraulic_diameter.value,
            friction_factor=self.friction_factor.value,
        )

        self.effective_area.value = effective_area_from_mass_flow(
            self.mass_flow.value,
            self.upstream_pressure.value - self.downstream_pressure.value,
            self.density.value,
        )

    @property
    def iteration_variables(self):
        return [self.mass_flow]

    @property
    def residuals(self):
        return [self.mass_flow.value - self._predicted_mass_flow]


class CircularPipeDarcyWeisbach(Component):

    def __init__(
        self,
        name: str,
        network: Network,
        upstream_pressure: State,
        downstream_pressure: State,
        length: float,
        inner_diameter: float,
        density: State,
        friction_factor: State | float | None = None,
        mass_flow: State | float | None = None,
        effective_area: State | None = None,
    ):
        self.setup()
        self._predicted_mass_flow = None

        self.hydraulic_diameter = self.inner_diameter
        self.cross_sectional_area = State(math.pi * self.inner_diameter.value**2 / 4.0)

        if not self.mass_flow.is_assigned:
            self.mass_flow.value = initial_mass_flow_guess(
                self.upstream_pressure.value - self.downstream_pressure.value,
                self.density.value,
                self.cross_sectional_area.value,
            )

    def evaluate_states(self):
        self.cross_sectional_area.value = math.pi * self.inner_diameter.value**2 / 4.0

        self._predicted_mass_flow = predicted_darcy_mass_flow(
            pressure_drop=self.upstream_pressure.value - self.downstream_pressure.value,
            density=self.density.value,
            length=self.length.value,
            hydraulic_diameter=self.hydraulic_diameter.value,
            friction_factor=self.friction_factor.value,
        )

        self.effective_area.value = effective_area_from_mass_flow(
            self.mass_flow.value,
            self.upstream_pressure.value - self.downstream_pressure.value,
            self.density.value,
        )

    @property
    def iteration_variables(self):
        return [self.mass_flow]

    @property
    def residuals(self):
        return [self.mass_flow.value - self._predicted_mass_flow]


class RectangularDuctDarcyWeisbach(Component):

    def __init__(
        self,
        name: str,
        network: Network,
        upstream_pressure: State,
        downstream_pressure: State,
        length: float,
        height: float,
        width: float,
        density: State,
        friction_factor: State | float | None = None,
        mass_flow: State | float | None = None,
        effective_area: State | None = None,
    ):
        self.setup()
        self._predicted_mass_flow = None

        self.cross_sectional_area = State()
        self.hydraulic_diameter = State()
        self.update_geometry()

        if not self.mass_flow.is_assigned:
            self.mass_flow.value = initial_mass_flow_guess(
                self.upstream_pressure.value - self.downstream_pressure.value,
                self.density.value,
                self.cross_sectional_area.value,
            )

    def update_geometry(self):
        height = self.height.value
        width = self.width.value

        area = height * width
        perimeter = 2.0 * (height + width)

        self.cross_sectional_area.value = area
        self.hydraulic_diameter.value = 4.0 * area / perimeter

    def evaluate_states(self):
        self.update_geometry()

        self._predicted_mass_flow = predicted_darcy_mass_flow(
            pressure_drop=self.upstream_pressure.value - self.downstream_pressure.value,
            density=self.density.value,
            length=self.length.value,
            hydraulic_diameter=self.hydraulic_diameter.value,
            friction_factor=self.friction_factor.value,
        )

        self.effective_area.value = effective_area_from_mass_flow(
            self.mass_flow.value,
            self.upstream_pressure.value - self.downstream_pressure.value,
            self.density.value,
        )

    @property
    def iteration_variables(self):
        return [self.mass_flow]

    @property
    def residuals(self):
        return [self.mass_flow.value - self._predicted_mass_flow]


class EllipticalDuctDarcyWeisbach(Component):

    def __init__(
        self,
        name: str,
        network: Network,
        upstream_pressure: State,
        downstream_pressure: State,
        length: float,
        semi_major_axis: float,
        semi_minor_axis: float,
        density: State,
        friction_factor: State | float | None = None,
        mass_flow: State | float | None = None,
        effective_area: State | None = None,
    ):
        self.setup()
        self._predicted_mass_flow = None

        self.cross_sectional_area = State()
        self.hydraulic_diameter = State()
        self.update_geometry()

        if not self.mass_flow.is_assigned:
            self.mass_flow.value = initial_mass_flow_guess(
                self.upstream_pressure.value - self.downstream_pressure.value,
                self.density.value,
                self.cross_sectional_area.value,
            )

    def update_geometry(self):
        a = self.semi_major_axis.value
        b = self.semi_minor_axis.value
        a, b = max(a, b), min(a, b)

        area = math.pi * a * b
        perimeter = math.pi * (
            3.0 * (a + b)
            - math.sqrt((3.0 * a + b) * (a + 3.0 * b))
        )

        self.cross_sectional_area.value = area
        self.hydraulic_diameter.value = 4.0 * area / perimeter

    def evaluate_states(self):
        self.update_geometry()

        self._predicted_mass_flow = predicted_darcy_mass_flow(
            pressure_drop=self.upstream_pressure.value - self.downstream_pressure.value,
            density=self.density.value,
            length=self.length.value,
            hydraulic_diameter=self.hydraulic_diameter.value,
            friction_factor=self.friction_factor.value,
        )

        self.effective_area.value = effective_area_from_mass_flow(
            self.mass_flow.value,
            self.upstream_pressure.value - self.downstream_pressure.value,
            self.density.value,
        )

    @property
    def iteration_variables(self):
        return [self.mass_flow]

    @property
    def residuals(self):
        return [self.mass_flow.value - self._predicted_mass_flow]


def predicted_darcy_mass_flow(
    pressure_drop: float,
    density: float,
    length: float,
    hydraulic_diameter: float,
    friction_factor: float,
) -> float:
    if length <= 0.0:
        raise ValueError(f"Darcy-Weisbach length must be positive. Got length={length}.")

    if hydraulic_diameter <= 0.0:
        raise ValueError("hydraulic_diameter must be positive.")

    if density <= 0.0:
        raise ValueError("density must be positive.")

    if friction_factor <= 0.0:
        raise ValueError("friction_factor must be positive.")

    if abs(pressure_drop) < 1e-12:
        return 0.0

    Kf = 8.0 * friction_factor * length / (
        density * math.pi**2 * hydraulic_diameter**5
    )

    return math.copysign(math.sqrt(abs(pressure_drop) / Kf), pressure_drop)


def initial_mass_flow_guess(
    pressure_drop: float,
    density: float,
    area: float,
) -> float:
    if abs(pressure_drop) < 1e-12:
        return 0.0

    if density <= 0.0:
        raise ValueError("density must be positive.")

    if area <= 0.0:
        raise ValueError("area must be positive.")

    return math.copysign(
        area * math.sqrt(2.0 * density * abs(pressure_drop)),
        pressure_drop,
    )


def effective_area_from_mass_flow(
    mass_flow: float,
    pressure_drop: float,
    density: float,
) -> float:
    """
    Computes an equivalent CdA/effective area from:

        mdot = CdA * sqrt(2 * rho * abs(dP))
    """
    if abs(pressure_drop) < 1e-12:
        return 0.0

    if density <= 0.0:
        raise ValueError("density must be positive.")

    return abs(mass_flow) / math.sqrt(2.0 * density * abs(pressure_drop))