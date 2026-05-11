from System import *
from Solvers import *


from constants import *

# --- Network Definition ---
SimpleNetwork = Network("Simple Network")

# --- State Definition ---

source_pressure = State(3e5)
source_temperature = State(300)
source_density = State()

line1_mdot = State()

manifold_pressure = State(2e5)
manifold_density = State()

line2_mdot = State()

atmospheric_pressure = State(101325)

# --- Fluid Definition ---

fluid = "RP-1"

source_fluid = DensityfromPT("Source Fluid", SimpleNetwork, "Water",
                             pressure=source_pressure,
                             temperature=source_temperature,
                             density=source_density)

manifold_fluid = DensityfromPT("Manifold Fluid", SimpleNetwork, "Water",
                             pressure=manifold_pressure,
                             temperature=source_temperature,
                             density=manifold_density)



# --- Component Definition ---

Source = IsothermalPressureBoundary("Source", SimpleNetwork,
                                    pressure=source_fluid.pressure,
                                    temperature=source_fluid.temperature,
                                    density=source_fluid.density)

avg_density = 0.5*(manifold_fluid.density + source_fluid.density)

Line1 = DischargeCoefficient("Line 1", SimpleNetwork,
                             upstream_pressure=source_fluid.pressure,
                             downstream_pressure=manifold_fluid.pressure,
                             density=avg_density,
                             discharge_coefficient=1,
                             cross_sectional_area=0.5e-4,
                             mass_flow=line1_mdot)

Manifold = IsothermalIncompressibleVolume("Manifold", SimpleNetwork,
                                          pressure=manifold_fluid.pressure,
                                          temperature=manifold_fluid.temperature,
                                          density=manifold_fluid.density,
                                          volume=0.01,
                                          mass_flow_in=line1_mdot,
                                          mass_flow_out=line2_mdot)

Line1 = DischargeCoefficient("Line 2", SimpleNetwork,
                             upstream_pressure=manifold_fluid.pressure,
                             downstream_pressure=atmospheric_pressure,
                             density=manifold_fluid.density,
                             discharge_coefficient=1,
                             cross_sectional_area=0.5e-4,
                             mass_flow=line2_mdot)

Ambient = PressureBoundary("Atmoshere", SimpleNetwork,
                           pressure=atmospheric_pressure)



print(SteadyState(SimpleNetwork).solve(return_type='dataframe', filename='solution.xlsx', verbose=True))

print(SimpleNetwork)
