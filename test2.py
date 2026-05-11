from System import *
from Solvers import *
from Utilities import Fluid


from constants import *


# --- Network Definition ---
SimpleNetwork = Network("Simple Network")

# --- State Definition ---

source_pressure = State(20 * PSIA_TO_PA)
source_temperature = State(300)

line1_mdot = State()

manifold_pressure = State(16 * PSIA_TO_PA)
manifold_density = State()

line2_mdot = State()

atmospheric_pressure = State(14.67 * PSIA_TO_PA)

# --- Fluid Definition ---

# make it so that if the user does provide a state for a fluid property, they can 
# call that value instead of just source_fluid.density
Fluid.add_alias("tone", "Acetone")
fluid = 'tone'

source_fluid = GeneralFluidLookupfromPT("Source Fluid", SimpleNetwork, fluid,
                             pressure=source_pressure,
                             temperature=source_temperature)




# --- Component Definition ---

Source = IsothermalPressureBoundary("Source", SimpleNetwork,
                                    pressure=source_fluid.pressure,
                                    temperature=source_fluid.temperature,
                                    density=source_fluid.density)

avg_density = 0.5*(manifold_density + source_fluid.density)



Line1 = DischargeCoefficient("Line 1", SimpleNetwork,
                             upstream_pressure=source_fluid.pressure,
                             downstream_pressure=manifold_pressure,
                             density=avg_density,
                             discharge_coefficient=1,
                             cross_sectional_area=0.5e-4,
                             mass_flow=line1_mdot)

manifold_fluid = GeneralFluidLookupfromPT("Manifold Fluid", SimpleNetwork, fluid,
                             pressure=manifold_pressure,
                             temperature=source_temperature,
                             density=manifold_density)




Manifold = IsothermalIncompressibleVolume("Manifold", SimpleNetwork,
                                          pressure=manifold_fluid.pressure,
                                          temperature=manifold_fluid.temperature,
                                          density=manifold_fluid.density,
                                          volume=0.01,
                                          mass_flow_in=line1_mdot,
                                          mass_flow_out=line2_mdot)

Line2 = DischargeCoefficient("Line 2", SimpleNetwork,
                             upstream_pressure=manifold_fluid.pressure,
                             downstream_pressure=atmospheric_pressure,
                             density=manifold_fluid.density,
                             discharge_coefficient=1,
                             cross_sectional_area=0.5e-4,
                             mass_flow=line2_mdot)

Ambient = PressureBoundary("Atmoshere", SimpleNetwork,
                           pressure=atmospheric_pressure)


print(SteadyState(SimpleNetwork).solve(return_type='dataframe', filename='solution.xlsx', verbose=False))