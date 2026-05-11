from System import *
from Solvers import *
from Utilities import Fluid


from constants import *

# --- Network Definition ---
SimpleNetwork = Network("Simple Network")

# --- Fluid Definition ---

Fluid.add_alias("air", "Air")
fluid = 'RP-1'

source_fluid = GeneralFluidLookupfromPT("Source Fluid", SimpleNetwork, fluid,
                             pressure=20 * PSIA_TO_PA,
                             temperature=300)

manifold_fluid = GeneralFluidLookupfromPT("Manifold Fluid", SimpleNetwork, fluid,
                             pressure= 15 * PSIA_TO_PA,
                             temperature=300)
                             #density=manifold_density)



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
                             cross_sectional_area=0.5e-4)

Manifold = IsothermalIncompressibleVolume("Manifold", SimpleNetwork,
                                          pressure=manifold_fluid.pressure,
                                          temperature=manifold_fluid.temperature,
                                          density=manifold_fluid.density,
                                          volume=0.01,
                                          mass_flow_in=Line1.mass_flow)

Line2 = DischargeCoefficient("Line 2", SimpleNetwork,
                             upstream_pressure=manifold_fluid.pressure,
                             downstream_pressure=14.67 * PSIA_TO_PA,
                             density=manifold_fluid.density,
                             discharge_coefficient=1,
                             cross_sectional_area=0.5e-4,
                             mass_flow=Manifold.mass_flow_out)

Ambient = PressureBoundary("Atmoshere", SimpleNetwork,
                           pressure=14.67 * PSIA_TO_PA)



print(SteadyState(SimpleNetwork).solve(return_type='dataframe', filename='solution.xlsx', verbose=False, static=True))
