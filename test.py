from System import *
from Solvers import *

from constants import *

SimpleNetwork = Network("Simple Network")

fluid = "RP-1"
source_pressure = State(3e5)
source_temperature = State(300)
source_density = State()
source_mdot_in = State(0)
line1_mdot = State()

manifold_pressure = State(2e5)
manifold_density = State()

line2_mdot = State()

ambient_pressure = State(101325)

Source = SimpleIsothermalVolume("Source", SimpleNetwork,
                                fluid=fluid,
                                pressure=source_pressure,
                                temperature=source_temperature,
                                density=source_density,
                                volume=0.05,
                                mass_flow_in=source_mdot_in,
                                mass_flow_out=line1_mdot)

avg_density1 = 0.5*(source_density + manifold_density)
Line1 = DischargeCoefficient("Line 1",
                                   network=SimpleNetwork,
                                   upstream_pressure=source_pressure,
                                   downstream_pressure=manifold_pressure,
                                   density=avg_density1,
                                   discharge_coefficient=1,
                                   cross_sectional_area=0.5e-4,
                                   mass_flow=line1_mdot)

Maniold = SimpleIsothermalVolume("Manifold", SimpleNetwork,
                                fluid=fluid,
                                pressure=manifold_pressure,
                                temperature=source_temperature,
                                density=manifold_density,
                                volume=0.02,
                                mass_flow_in=line1_mdot,
                                mass_flow_out=line2_mdot)

Line2 = DischargeCoefficient("Line 2",
                                   network=SimpleNetwork,
                                   upstream_pressure=manifold_pressure,
                                   downstream_pressure=ambient_pressure,
                                   density=manifold_density,
                                   discharge_coefficient=1,
                                   cross_sectional_area=0.5e-4,
                                   mass_flow=line2_mdot)

Ambient = PressureNode("Atmosphere", 
                       network=SimpleNetwork,
                       pressure=ambient_pressure)

print(SteadyState(SimpleNetwork).solve(return_type='dataframe', filename='solution.xlsx', verbose=False))