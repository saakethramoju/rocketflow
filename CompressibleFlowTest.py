from System import *
from Solvers import *
from Utilities import FluidRegistry

from constants import *

# --- Network Definition ---
CompressibleFlowNetwork = Network("Compressible Flow Network")

# --- Fluid Definition ---
fluid = "Pressurant"
FluidRegistry.add_alias(fluid, "N2")

supply_gas = IdealGasLookup("Supply Gas", CompressibleFlowNetwork, fluid, 
                            pressure=6e5, temperature=300)

ullage_gas = IdealGasLookup("Tank gas", CompressibleFlowNetwork, fluid,
                          pressure=2e5, temperature=300)

# --- Component Definition ---
PressReg = IsentropicGasRegulator("PR-Reg", CompressibleFlowNetwork,
                        upstream_pressure=supply_gas.pressure,
                        upstream_temperature=supply_gas.temperature,
                        set_pressure=ullage_gas.pressure,
                        discharge_coefficient=1,
                        cross_sectional_area=(3.14 * 0.125**2) * IN2_TO_M2,
                        specific_gas_constant=8.314/supply_gas.molar_mass,
                        specific_heat_ratio=supply_gas.specific_heat_ratio)

Ullage = SimpleVolume("Tank", CompressibleFlowNetwork, 
                        pressure=ullage_gas.pressure,
                        temperature=ullage_gas.temperature,
                        volume=0.01,
                        mass_flow_in=PressReg.mass_flow,
                        mass_flow_out=0.01,
                        enthalpy_in=supply_gas.enthalpy,
                        enthalpy_out=ullage_gas.enthalpy)



RegCd_Balance = Balance("Balance reg Cd to give a certain set pressure", CompressibleFlowNetwork,
                        variable=PressReg.discharge_coefficient,
                        function=PressReg.set_pressure - 5e5)


print(SteadyState(CompressibleFlowNetwork).solve(return_type='dataframe', filename='solution.xlsx', verbose=True, static=False))
