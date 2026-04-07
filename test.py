from __future__ import annotations

import numpy as np
from typing import TYPE_CHECKING
from scipy.optimize import root

from System import *
from Components import DischargeCoefficient, PressureNode, Component, Variable, State
from Fluid import Fluid


class TestVolume(Component):

    def __init__(self,
                 name: str,
                 pressure: State | float,
                 density: State | float,
                 volume: State | float,
                 mass_flow_in: State | float,
                 mass_flow_out: State | float):
        super().__init__(name)

        self.p = Variable(pressure)
        self.rho = Variable(density)
        self.V = Variable(volume)
        self.mdot_in = Variable(mass_flow_in)
        self.mdot_out = Variable(mass_flow_out)

    def pre_state_evaluation(self) -> None:
        self.rho.value = 998 * (1 + (self.p.value - 101325) / (2.2e9))

    def evaluate_states(self) -> None:
        pass

    @property
    def iteration_variables(self) -> list[Variable]:
        return [self.p]

    @property
    def residuals(self) -> list[float]:
        return [self.mdot_in.value - self.mdot_out.value]


if __name__ == "__main__":

    source_pressure = State(2e5)
    drain_pressure = State(101325)

    manifold_pressure = State(2e5)
    manifold_density = State(1000)
    mdot_in = State(0)
    mdot_out = State(0)


    vol = TestVolume("manifold",
                     pressure=manifold_pressure,
                     density=manifold_density,
                     volume=0.1287,
                     mass_flow_in=mdot_in,
                     mass_flow_out=mdot_out)
    
    source = PressureNode("in",
                          pressure=source_pressure)
    
    drain = PressureNode("atm",
                         pressure=drain_pressure)
    
    line_in = DischargeCoefficient("wendell",
                                   upstream_pressure=source_pressure,
                                   downstream_pressure=manifold_pressure,
                                   density=manifold_density,
                                   discharge_coefficient=1,
                                   cross_sectional_area=0.5e-4,
                                   mass_flow=mdot_in)
    
    line_out = DischargeCoefficient("wendell 2",
                                   upstream_pressure=manifold_pressure,
                                   downstream_pressure=drain_pressure,
                                   density=manifold_density,
                                   discharge_coefficient=1,
                                   cross_sectional_area=0.5e-4,
                                   mass_flow=mdot_out)
    

    
    
    # what network will do
    components = [source, drain, vol, line_in, line_out]

    def collect_iteration_variables(component_list: list[Component]) -> list[Variable]:
        iter_vars = []
        for comp in component_list:
            iter_vars.extend(comp.iteration_variables)
        return iter_vars


    def collect_guess_values(iter_var_list: list[Variable]) -> list[float]:
        return [var.value for var in iter_var_list]


    def assign_iteration_values(iteration_values: list[float], iter_var_list: list[Variable]) -> None:
        if len(iteration_values) != len(iter_var_list):
            raise ValueError(
                f"Length mismatch: got {len(iteration_values)} iteration values "
                f"but expected {len(iter_var_list)}"
            )
        for val, var in zip(iteration_values, iter_var_list):
            var.value = val


    def collect_residuals(component_list: list[Component]) -> list[float]:
        return [
            float(r)
            for comp in component_list
            for r in comp.residuals
        ]


    # solver procedure:
    iter_vars = collect_iteration_variables(components)
    x0 = collect_guess_values(iter_vars)

    for c in components:
        c.pre_state_evaluation()

    def system_residuals(x):
        assign_iteration_values(x, iter_vars)

        for c in components:
            c.pre_state_evaluation()

        for c in components:
            c.evaluate_states()

        return collect_residuals(components)


    sol = root(system_residuals, x0)

    assign_iteration_values(sol.x, iter_vars)
    for c in components:
        c.evaluate_states()

    print("success =", sol.success)
    print("message =", sol.message)
    print("x =", sol.x)
    print("manifold_pressure =", manifold_pressure.value)
    print("manifold_density =", manifold_density.value)
    print("mdot_in =", mdot_in.value)
    print("mdot_out =", mdot_out.value)
    print("residuals =", system_residuals(sol.x))