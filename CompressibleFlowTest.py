from System import *
from Solvers import *
from Utilities import Fluid

from constants import *

# --- Network Definition ---
CompressibleFlow = Network("Compressible Flow Network")

# --- Fluid Definition ---
fluid = 'Nitrogen'
supply_gas = GeneralFluidLookupfromPT("Supply Gas", CompressibleFlow, fluid,
                                      pressure=3000 * PSIA_TO_PA,
                                      temperature=300)

supply_gas_specific_heat = IdealGasSpecificHeats("Supply Gas Specific Heat", CompressibleFlow,
                                                 supply_gas.molar_mass, degrees_of_freedom=5)

tank_gas = GeneralFluidLookupfromPT("Tank Gas", CompressibleFlow, fluid,
                                      pressure=2000 * PSIA_TO_PA,
                                      temperature=300)

# --- Component Definition ---
PressReg = IsentropicGasRegulator("PR-Reg", CompressibleFlow,
                        upstream_pressure=supply_gas.pressure,
                        upstream_temperature=supply_gas.temperature,
                        set_pressure=tank_gas.pressure,
                        discharge_coefficient=1,
                        cross_sectional_area=(3.14 * 0.125**2) * IN2_TO_M2,
                        specific_gas_constant=8.314/supply_gas.molar_mass,
                        specific_heat_ratio=supply_gas_specific_heat.specific_heat_ratio)
'''
RegCd_Balance = Balance("Balance reg Cd to give a certain set pressure", CompressibleFlow,
                        variable=PressReg.discharge_coefficient,
                        function=PressReg.set_pressure - 1000*PSIA_TO_PA)
'''

#Tank = 

print(SteadyState(CompressibleFlow).solve(return_type='dataframe', filename='solution.xlsx', verbose=True, static=False))