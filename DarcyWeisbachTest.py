from System import *
from Solvers import *
from Utilities import Fluid

from constants import *

# --- Network Definition ---
SimpleNetwork = Network("Simple Network")

# these are just here to show that external states work
atmospheric_pressure = State(14.67 * PSIA_TO_PA)
manifold_density = State()

# --- Fluid Definition ---
#Fluid.add_alias("Octane", "n-Octane")
fluid = 'RP-1'

source_fluid = FluidLookup("Source Fluid", SimpleNetwork, fluid,
                             pressure=20 * PSIA_TO_PA,
                             enthalpy=Fluid(fluid,pressure=20 * PSIA_TO_PA, temperature=300).enthalpy)

manifold_fluid = FluidLookup("Manifold Fluid", SimpleNetwork, fluid,
                             pressure= 10 * PSIA_TO_PA,
                             temperature=300,
                             density=manifold_density)
# --- Component Definition ---

Source = IsothermalPressureBoundary("Source", SimpleNetwork,
                                    pressure=source_fluid.pressure,
                                    temperature=source_fluid.temperature)

avg_density = 0.5*(manifold_density + source_fluid.density)

Line1 = CircularPipeDarcyWeisbach("Line 1", SimpleNetwork,
                      upstream_pressure=Source.pressure,
                      downstream_pressure=manifold_fluid.pressure,
                      length=1,
                      inner_diameter=0.5 * IN_TO_M,
                      density=avg_density,
                      dynamic_viscosity=source_fluid.dynamic_viscosity,
                      roughness=0.1e-3)

'''
Manifold = IsothermalIncompressibleVolume("Manifold", SimpleNetwork,
                                          pressure=manifold_fluid.pressure,
                                          temperature=manifold_fluid.temperature,
                                          density=manifold_fluid.density,
                                          volume=0.01,
                                          mass_flow_in=Line1.mass_flow)
'''

Manifold = SimpleVolume("Manifold", SimpleNetwork, 
                        pressure=manifold_fluid.pressure,
                        temperature=manifold_fluid.temperature,
                        density=manifold_fluid.density,
                        volume=0.01,
                        mass_flow_in=Line1.mass_flow,
                        enthalpy_in=source_fluid.enthalpy,
                        enthalpy_out=manifold_fluid.enthalpy)
'''
Line2 = DischargeCoefficient("Line 2", SimpleNetwork,
                             upstream_pressure=manifold_fluid.pressure,
                             downstream_pressure=14.67 * PSIA_TO_PA,
                             density=manifold_fluid.density,
                             discharge_coefficient=1,
                             cross_sectional_area=0.5e-4,
                             mass_flow=Manifold.mass_flow_out)
'''
Line2 = EllipticalDuctDarcyWeisbach("Line 2", SimpleNetwork,
                      upstream_pressure=Manifold.pressure,
                      downstream_pressure=atmospheric_pressure,
                      length=2,
                      semi_major_axis=2 * IN_TO_M,
                      semi_minor_axis=1 * IN_TO_M,
                      density=manifold_density,
                      dynamic_viscosity=manifold_fluid.dynamic_viscosity,
                      roughness=0.1e-3,
                      mass_flow=Manifold.mass_flow_out)

Ambient = PressureBoundary("Atmoshere", SimpleNetwork,
                           pressure=atmospheric_pressure)

source_pressure_balance = Balance("Balance source pressure until mdot = 0.6",
                          network=SimpleNetwork,
                          variable=Line1.upstream_pressure,
                          function=Line2.mass_flow - 0.6,)
                          #bounds=(0, 3e5),
                          #keep_feasible=True)

print(SteadyState(SimpleNetwork).solve(return_type='dataframe', filename='solution.xlsx', verbose=True, static=False))

#Fluid.show_available_fluids()
