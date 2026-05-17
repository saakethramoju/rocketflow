from System import *
from Solvers import *
from Utilities import FluidRegistry

from constants import *

# --- Network Definition ---
TankNetwork = Network("Pressurized Tank System")

atmospheric_pressure = 101325
chamber_pressure = State(285 * PSIA_TO_PA)


# ---- Fluids -----
pressurant = 'gn2'
oxidizer = 'lox'
fuel = 'rp-1'

PressurantSupply = FluidLookup("COPV Gas", TankNetwork, pressurant, pressure=6000*PSIA_TO_PA, temperature=300)

FuelUllageGas = IdealGasLookup("Fuel Tank Ullage Gas", TankNetwork, pressurant, pressure=400*PSIA_TO_PA, temperature=300)
OxUllageGas = IdealGasLookup("Ox Tank Ullage Gas", TankNetwork, pressurant, pressure=400*PSIA_TO_PA, temperature=300)

FuelTankFluid = FluidLookup("Tank Fuel Liquid", TankNetwork, fuel, 400*PSIA_TO_PA, temperature=300)
OxTankFluid = FluidLookup("Tank Oxidizer Liquid", TankNetwork, oxidizer, 400*PSIA_TO_PA, temperature=90)

FuelManifoldFluid = FluidLookup("Fuel Manifold Liquid", TankNetwork, fuel, pressure=300*PSIA_TO_PA, temperature=300, 
                                flash_values=("pressure", "enthalpy"))

OxManifoldFluid = FluidLookup("Ox Manifold Liquid", TankNetwork, oxidizer, pressure=300*PSIA_TO_PA, temperature=90, 
                                flash_values=("pressure", "enthalpy"))
# ----- Components ------

COPV = Boundary("COPV", TankNetwork, PressurantSupply.pressure, PressurantSupply.temperature)

FuelBangBang = IsentropicGasRegulator("FBB", TankNetwork, COPV.pressure, COPV.temperature,
                                             set_pressure=FuelUllageGas.pressure, discharge_coefficient=1,
                                             cross_sectional_area=0.0490625*IN2_TO_M2, specific_gas_constant=FuelUllageGas.gas_constant,
                                             specific_heat_ratio=OxUllageGas.specific_heat_ratio)

OxBangBang = IsentropicGasRegulator("OBB", TankNetwork, COPV.pressure, COPV.temperature,
                                             set_pressure=OxUllageGas.pressure, discharge_coefficient=1,
                                             cross_sectional_area=0.0490625*IN2_TO_M2, specific_gas_constant=OxUllageGas.gas_constant,
                                             specific_heat_ratio=OxUllageGas.specific_heat_ratio)

FuelTank = PressurizedTank("Fuel Tank", TankNetwork, pressure=FuelTankFluid.pressure, pressurant_density=FuelUllageGas.density, 
                         liquid_density=FuelTankFluid.density, collapse_factor=1, ullage_temperature=FuelUllageGas.temperature, 
                         liquid_temperature=FuelTankFluid.temperature, mass_flow_in=FuelBangBang.mass_flow)

OxTank = PressurizedTank("Ox Tank", TankNetwork, pressure=OxTankFluid.pressure, pressurant_density=OxUllageGas.density, 
                         liquid_density=OxTankFluid.density, collapse_factor=1.4, ullage_temperature=OxUllageGas.temperature, 
                         liquid_temperature=OxTankFluid.temperature, mass_flow_in=OxBangBang.mass_flow)

FuelGravitydP = GravityPressureChange("Fuel Runline Gravity dP", TankNetwork, upstream_pressure=FuelTank.pressure, density=FuelTank.liquid_density, 
                                      elevation_change=-2.5)

OxGravitydP = GravityPressureChange("Ox Runline Gravity dP", TankNetwork, upstream_pressure=OxTank.pressure, density=OxTank.liquid_density, 
                                      elevation_change=-2)

FuelRunline = GenericDarcyWeisbach("Fuel Main Line", TankNetwork, upstream_pressure=FuelGravitydP.downstream_pressure, 
                                   downstream_pressure=FuelManifoldFluid.pressure, length=2.5, cross_sectional_area=0.0490625*IN2_TO_M2,
                                   hydraulic_diameter=0.5*IN_TO_M, roughness=0.1e-3, density=FuelTank.liquid_density, 
                                   dynamic_viscosity=FuelTankFluid.dynamic_viscosity, mass_flow=FuelTank.mass_flow_out)

OxRunline = GenericDarcyWeisbach("Ox Main Line", TankNetwork, upstream_pressure=OxGravitydP.downstream_pressure, 
                                   downstream_pressure=OxManifoldFluid.pressure, length=2, cross_sectional_area=0.0490625*IN2_TO_M2,
                                   hydraulic_diameter=0.5*IN_TO_M, roughness=0.1e-3, density=OxTank.liquid_density, 
                                   dynamic_viscosity=OxTankFluid.dynamic_viscosity, mass_flow=OxTank.mass_flow_out)

FuelManifold = SimpleVolume("Fuel Injector Manifold", TankNetwork, pressure=FuelManifoldFluid.pressure, enthalpy=FuelManifoldFluid.enthalpy,
                            volume = 685*IN2_TO_M2, temperature=FuelManifoldFluid.temperature, density=FuelManifoldFluid.density, 
                            mass_flow_in=FuelRunline.mass_flow, enthalpy_in=FuelTankFluid.enthalpy)

OxManifold = SimpleVolume("Ox Injector Manifold", TankNetwork, pressure=OxManifoldFluid.pressure, enthalpy=OxManifoldFluid.enthalpy,
                            volume=685*IN2_TO_M2, temperature=OxManifoldFluid.temperature, density=OxManifoldFluid.density, 
                            mass_flow_in=OxRunline.mass_flow, enthalpy_in=OxTankFluid.enthalpy)


FuelInjector = DischargeCoefficient("Fuel Injector Orifices", TankNetwork, upstream_pressure=FuelManifold.pressure,
                                    downstream_pressure=chamber_pressure, density=FuelManifold.density, discharge_coefficient=1, 
                                    cross_sectional_area=0.5e-4, mass_flow=FuelManifold.mass_flow_out)

OxInjector = DischargeCoefficient("Ox Injector Orifices", TankNetwork, upstream_pressure=OxManifold.pressure,
                                    downstream_pressure=chamber_pressure, density=OxManifold.density, discharge_coefficient=1, 
                                    cross_sectional_area=1.0e-4, mass_flow=OxManifold.mass_flow_out)

MainChamber = MainCombustionChamber("MCC", TankNetwork, chamber_pressure=chamber_pressure, oxidizer_mass_flow=OxInjector.mass_flow,
                                    fuel_mass_flow=FuelInjector.mass_flow)

Nozzle = RocketCEAChokedNozzle("Nozzle", TankNetwork, fuel=fuel, oxidizer=oxidizer, chamber_pressure=MainChamber.chamber_pressure,
                               throat_area=5.75*IN2_TO_M2, expansion_ratio=4, 
                               mixture_ratio=MainChamber.oxidizer_mass_flow/MainChamber.fuel_mass_flow,
                               ambient_pressure=atmospheric_pressure, characterstic_velocity_efficiency=0.9,
                               thrust_coefficient_efficiency=0.95, mass_flow=MainChamber.nozzle_mass_flow)

'''
fuel_tank_balance = Balance("Fuel tank balance", TankNetwork, 
                            variable=FuelBangBang.discharge_coefficient,
                            function=FuelTank.pressure - 450*PSIA_TO_PA)

ox_tank_balance = Balance("Ox tank balance", TankNetwork, 
                            variable=OxBangBang.discharge_coefficient,
                            function=OxTank.pressure - 400*PSIA_TO_PA)
'''
Pc_balance = Balance("Chamber Pressure balance", TankNetwork,
                     variable=OxBangBang.discharge_coefficient,
                     function=MainChamber.chamber_pressure - 300*PSIA_TO_PA)

MR_balance = Balance("Mixture Ratio balance", TankNetwork,
                     variable=FuelBangBang.discharge_coefficient,
                     function=Nozzle.mixture_ratio - 1.5)

thrust_balance = Balance("Thrust balance", TankNetwork,
                     variable=Nozzle.throat_area,
                     function=Nozzle.thrust - 2300*LBF_TO_N)

solution = SteadyState(TankNetwork).solve(return_type='dataframe', verbose=True, static=False)
print(solution.to_string(index=False))


print("\n" + "="*50)
print("        ENGINE SYSTEM SUMMARY")
print("="*50)

print(
    f"Fuel Tank Pressure      : "
    f"{FuelTank.pressure.value * PA_TO_PSIA:10.3f} psia"
)

print(
    f"Ox Tank Pressure        : "
    f"{OxTank.pressure.value * PA_TO_PSIA:10.3f} psia"
)

print(
    f"Chamber Pressure        : "
    f"{MainChamber.chamber_pressure.value * PA_TO_PSIA:10.3f} psia"
)

print(
    f"Mixture Ratio (O/F)     : "
    f"{Nozzle.mixture_ratio.value:10.4f}"
)

print(
    f"Fuel Regulator Cd       : "
    f"{FuelBangBang.discharge_coefficient.value:10.5f}"
)

print(
    f"Ox Regulator Cd         : "
    f"{OxBangBang.discharge_coefficient.value:10.5f}"
)

print(
    f"Thrust                  : "
    f"{Nozzle.thrust.value*N_TO_LBF:10.3f} lbf"
)
print("="*50)