from System import *
from Solvers import *

from constants import *

import numpy as np


class TestDarcy(Component):

    def __init__(
        self,
        name: str,
        network: Network,
        upstream_pressure: State,
        downstream_pressure: State,
        length: float,
        inner_diameter: float,
        density: State,
        mass_flow: float | None = None,
        friction_factor: State | None = None,
    ):
        self.setup()

        if not self.mass_flow.is_assigned:
            rho = self.density.value
            dP = self.upstream_pressure.value - self.downstream_pressure.value
            area = np.pi * self.inner_diameter.value**2 / 4.0
            self.mass_flow.value = np.sign(dP) * area * np.sqrt(
                2.0 * rho * abs(dP)
            )

        self._predicted_mass_flow = None

    def evaluate_states(self):
        rho = self.density.value
        Dh = self.inner_diameter.value
        L = self.length.value
        dP = self.upstream_pressure.value - self.downstream_pressure.value
        f = self.friction_factor.value

        Kf = 8.0 * f * L / (rho * np.pi**2 * Dh**5)
        self._predicted_mass_flow = np.sign(dP) * np.sqrt(np.abs(dP) / Kf)

    @property
    def iteration_variables(self):
        return [self.mass_flow]

    @property
    def residuals(self):
        return [self.mass_flow.value - self._predicted_mass_flow]


class Churchill(Component):

    def __init__(
        self,
        name: str,
        network: Network,
        mass_flow: State,
        hydraulic_diameter: float,
        dynamic_viscosity: State,
        cross_sectional_area: float,
        roughness: float = 0.0,
        friction_factor: float | None = 0.04,
        reynolds_number: float | None = None,
    ):
        self.setup()

        if not self.reynolds_number.is_assigned:
            self.reynolds_number.value = (
                abs(self.mass_flow.value)
                * self.hydraulic_diameter.value
                / (self.dynamic_viscosity.value * self.cross_sectional_area.value)
            )

        if not self.friction_factor.is_assigned:
            Re = max(self.reynolds_number.value, 1e-12)
            eps_D = self.roughness.value / self.hydraulic_diameter.value

            A = (
                2.457
                * np.log(
                    1.0 / ((7.0 / Re) ** 0.9 + 0.27 * eps_D)
                )
            ) ** 16

            B = (37530.0 / Re) ** 16

            self.friction_factor.value = 8.0 * (
                (8.0 / Re) ** 12
                + 1.0 / (A + B) ** 1.5
            ) ** (1.0 / 12.0)

        

    def evaluate_states(self):
        self.reynolds_number.value = (
            abs(self.mass_flow.value)
            * self.hydraulic_diameter.value
            / (self.dynamic_viscosity.value * self.cross_sectional_area.value)
        )

    @property
    def iteration_variables(self):
        return [self.friction_factor]

    @property
    def residuals(self):
        f = self.friction_factor.value
        Re = max(self.reynolds_number.value, 1e-12)
        eps = self.roughness.value
        Dh = self.hydraulic_diameter.value

        relative_roughness = eps / Dh

        A = (
            2.457
            * np.log(
                1.0
                / ((7.0 / Re) ** 0.9 + 0.27 * relative_roughness)
            )
        ) ** 16

        B = (37530.0 / Re) ** 16

        f_churchill = 8.0 * (
            (8.0 / Re) ** 12
            + 1.0 / (A + B) ** 1.5
        ) ** (1.0 / 12.0)

        return [f - f_churchill]


if __name__ == "__main__":

    network = Network("Simple Darcy Test")

    atmospheric_pressure = State(14.67 * PSIA_TO_PA)

    fluid = "RP-1"

    source_fluid = FluidLookup(
        "Source Fluid",
        network,
        fluid,
        pressure=20.0 * PSIA_TO_PA,
        temperature=300.0,
    )

    drain_fluid = FluidLookup(
        "Drain Fluid",
        network,
        fluid,
        pressure=atmospheric_pressure,
        temperature=300.0,
    )

    Source = IsothermalPressureBoundary(
        "Source",
        network,
        pressure=source_fluid.pressure,
        temperature=source_fluid.temperature,
    )

    diameter = 0.5 * IN_TO_M
    area = np.pi * diameter**2 / 4.0


    Line1 = TestDarcy(
        "Line 1 Darcy",
        network,
        upstream_pressure=source_fluid.pressure,
        downstream_pressure=drain_fluid.pressure,
        length=1.0,
        inner_diameter=diameter,
        density=source_fluid.density,
    )

    Line1F = Churchill(
        "Line 1 Churchill Friction",
        network,
        mass_flow=Line1.mass_flow,
        hydraulic_diameter=diameter,
        dynamic_viscosity=source_fluid.dynamic_viscosity,
        cross_sectional_area=area,
        roughness=0.1e-4,
        friction_factor=Line1.friction_factor,
    )

    print(
        SteadyState(network).solve(
            return_type="dataframe",
            verbose=True,
            static=False,
        )
    )


    '''
    # ---- Pump Data -----
    design_rotor_speed = 20000.0       # rpm
    design_volumetric_flow = 0.0024    # m^3/s
    design_head_rise = 46.5            # m
    design_torque = 0.75               # N-m

    g = 9.80665
    density = 750.0                    # kg/m^3
    impeller_diameter = 0.040          # m


    omega_design = np.pi / 30.0 * design_rotor_speed

    design_flow_coefficient = (
        design_volumetric_flow
        / (omega_design * impeller_diameter**3)
    )

    design_head_coefficient = (
        design_head_rise
        / (omega_design**2 * impeller_diameter**2 / g)
    )

    design_torque_coefficient = (
        design_torque
        / (density * omega_design**2 * impeller_diameter**5)
    )


    normalized_flow_coefficient_map = [
        0.00,
        0.40,
        0.70,
        1.00,
        1.20,
        1.40,
    ]

    # reverse-sigmoid-ish head curve
    normalized_head_coefficient_map = [
        1.35,
        1.28,
        1.15,
        1.00,
        0.78,
        0.45,
    ]

    # torque forced to zero at zero flow and increases with flow
    normalized_torque_coefficient_map = [
        0.00,
        0.35,
        0.65,
        1.00,
        1.22,
        1.45,
    ]


    # ---- Components ----

    EPumpMap = TurboMap("E-Pump Map", MapNetwork, 
                        rotor_speed=25000,
                        volumetric_flow=2/density,
                        density=density,
                        flow_geometric_parameter=impeller_diameter,
                        head_geometric_parameter=impeller_diameter,
                        torque_geometric_parameter=impeller_diameter,
                        design_flow_coefficient=design_flow_coefficient,
                        design_head_coefficient=design_head_coefficient,
                        design_torque_coefficient=design_torque_coefficient,
                        normalized_flow_coefficient_map=normalized_flow_coefficient_map,
                        normalized_head_coefficient_map=normalized_head_coefficient_map,
                        normalized_torque_coefficient_map=normalized_torque_coefficient_map,
                        )

    solution = SteadyState(MapNetwork).solve(return_type='dataframe', verbose=True, static=False)
    print(solution.to_string(index=False))
    '''

    '''
    import numpy as np
    import matplotlib.pyplot as plt

    gamma = 1.4

    M = np.linspace(0.01, 5, 1000)

    F = M * (1 + (gamma - 1)/2 * M**2)**(
        -(gamma + 1)/(2 * (gamma - 1))
    )

    plt.figure(figsize=(8,5))
    plt.plot(M, F)

    plt.axvline(1, linestyle="--")
    plt.xlabel("Mach Number")
    plt.ylabel("F(M)")
    plt.title("Compressible Mass Flow Function")
    plt.grid(True)

    plt.show()
    '''


