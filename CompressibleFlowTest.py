from System import *
from Solvers import *
from Utilities import FluidRegistry

from constants import *

# --- Network Definition ---
CompressibleFlowNetwork = Network("Compressible Flow Network")

atmospheric_pressure = 101325

# --- Fluid Definition ---
pressurant = "Pressurant"
FluidRegistry.add_alias(pressurant, "N2")

liquid = "Water"

supply_gas = IdealGasLookup("Supply Gas", CompressibleFlowNetwork, pressurant, 
                            pressure=6e5, temperature=300)

ullage_gas = IdealGasLookup("Tank gas", CompressibleFlowNetwork, pressurant,
                          pressure=2e5, temperature=300)

tank_liquid = FluidLookup("Tank Fluid", CompressibleFlowNetwork, liquid,
                          pressure=ullage_gas.pressure, temperature=ullage_gas.temperature)

# --- Component Definition ---
PressReg = IsentropicGasRegulator("PR-Reg", CompressibleFlowNetwork,
                        upstream_pressure=supply_gas.pressure,
                        upstream_temperature=supply_gas.temperature,
                        set_pressure=ullage_gas.pressure,
                        discharge_coefficient=0.3,
                        cross_sectional_area=(3.14 * 0.125**2) * IN2_TO_M2,
                        specific_gas_constant=8.314/supply_gas.molar_mass,
                        specific_heat_ratio=supply_gas.specific_heat_ratio)


TankUllage = Boundary("Ullage", CompressibleFlowNetwork, 
                  pressure=ullage_gas.pressure, 
                  temperature=ullage_gas.temperature)


TankLiquid = Boundary("Liquid", CompressibleFlowNetwork,
                      pressure=TankUllage.pressure,
                      temperature=TankUllage.temperature)
                  

Runline = CircularPipeDarcyWeisbach("Line 1", CompressibleFlowNetwork,
                      upstream_pressure=TankLiquid.pressure,
                      downstream_pressure=atmospheric_pressure,
                      length=5,
                      inner_diameter=0.5 * IN_TO_M,
                      density=tank_liquid.density,
                      dynamic_viscosity=tank_liquid.dynamic_viscosity,
                      roughness=0.1e-3)

'''
RegCd_Balance = Balance("Balance reg Cd to give a certain set pressure", CompressibleFlowNetwork,
                        variable=PressReg.discharge_coefficient,
                        function=PressReg.set_pressure - 5e5)
'''

TankPressureBalance = Balance(
                        "Tank pressure volume balance",
                        CompressibleFlowNetwork,
                        variable=TankUllage.pressure,
                        function=PressReg.mass_flow / ullage_gas.density - Runline.mass_flow / tank_liquid.density)

print(SteadyState(CompressibleFlowNetwork).solve(return_type='dataframe', filename='solution.xlsx', verbose=True, static=False))
