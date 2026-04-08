from __future__ import annotations

from scipy.optimize import root

from System import *
from Solvers import *


tripleT = Network("tripleT")

source_pressure = State(2e5)
drain_pressure = State(101325)

manifold_pressure = State(2e5)
manifold_density = State(1000)
mdot_in = State(0)
mdot_out = State(0)


vol = SimpleIncompressibleVolume("manifold",
                    network=tripleT,
                    pressure=manifold_pressure,
                    density=manifold_density,
                    volume=0.1287,
                    mass_flow_in=mdot_in,
                    mass_flow_out=mdot_out)

source = PressureNode("in",
                        network=tripleT,
                        pressure=source_pressure)

drain = PressureNode("atm",
                        network=tripleT,
                        pressure=drain_pressure)

line_in = DischargeCoefficient("wendell",
                                network=tripleT,
                                upstream_pressure=source_pressure,
                                downstream_pressure=manifold_pressure,
                                density=manifold_density,
                                discharge_coefficient=1,
                                cross_sectional_area=0.5e-4,
                                mass_flow=mdot_in)

line_out = DischargeCoefficient("wendell 2",
                                network=tripleT,
                                upstream_pressure=manifold_pressure,
                                downstream_pressure=drain_pressure,
                                density=manifold_density,
                                discharge_coefficient=1,
                                cross_sectional_area=0.5e-4,
                                mass_flow=mdot_out)



if __name__ == "__main__":

    print(tripleT)
    solver = SteadyState(tripleT)
    df = solver.solve(return_type='dataframe')
    print(df)

