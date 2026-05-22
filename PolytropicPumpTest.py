from System import *
from Solvers import *

from constants import *


# --- Network Definition ---
PumpNetwork = Network("Pumped Sytem")

# these are useful external states,
# which is why they're at the top
fuel_shaft_speed = State(25000)
ox_shaft_speed = State(25000)

chamber_pressure = State(285 * PSIA_TO_PA)
atmospheric_pressure = 14.67 * PSIA_TO_PA


# ---- Fluids ----
pressurant = "gn2"
fuel = "rp-1"
oxidizer = "lox"


PressurantSupply = FluidLookup(
    "COPV Gas",
    PumpNetwork,
    pressurant,
    pressure=6000 * PSIA_TO_PA,
    temperature=300,
)

FuelUllageGas = IdealGasLookup(
    "Fuel Tank Ullage Gas",
    PumpNetwork,
    pressurant,
    pressure=80 * PSIA_TO_PA,
    temperature=300,
)

OxUllageGas = IdealGasLookup(
    "Ox Tank Ullage Gas",
    PumpNetwork,
    pressurant,
    pressure=80 * PSIA_TO_PA,
    temperature=300,
)

FuelTankFluid = FluidLookup(
    "Tank Fuel Liquid",
    PumpNetwork,
    fuel,
    pressure=State(80 * PSIA_TO_PA),
    temperature=300,
)

OxTankFluid = FluidLookup(
    "Tank Oxidizer Liquid",
    PumpNetwork,
    oxidizer,
    pressure=State(80 * PSIA_TO_PA),
    temperature=90,
)

FuelPumpInletFluid = FluidLookup(
    "Fuel Pump Inlet Fluid",
    PumpNetwork,
    fuel,
    pressure=State(75 * PSIA_TO_PA),
    temperature=293.15,
    flash_values=("pressure", "enthalpy"),
)

OxPumpInletFluid = FluidLookup(
    "Ox Pump Inlet Fluid",
    PumpNetwork,
    oxidizer,
    pressure=State(75 * PSIA_TO_PA),
    temperature=90,
    flash_values=("pressure", "enthalpy"),
)

FuelPumpDischargeFluid = FluidLookup(
    "Fuel Pump Outlet Fluid",
    PumpNetwork,
    fuel,
    pressure=State(430 * PSIA_TO_PA),
    temperature=293.15,
    flash_values=("pressure", "enthalpy"),
)

OxPumpDischargeFluid = FluidLookup(
    "Ox Pump Outlet Fluid",
    PumpNetwork,
    oxidizer,
    pressure=State(335 * PSIA_TO_PA),
    temperature=90,
    flash_values=("pressure", "enthalpy"),
)


FuelManifoldFluid = FluidLookup(
    "Fuel Manifold Liquid",
    PumpNetwork,
    fuel,
    pressure=State(410 * PSIA_TO_PA),
    temperature=293.15,
    flash_values=("pressure", "enthalpy"),
)

OxManifoldFluid = FluidLookup(
    "Ox Manifold Liquid",
    PumpNetwork,
    oxidizer,
    pressure=State(320 * PSIA_TO_PA),
    temperature=90,
    flash_values=("pressure", "enthalpy"),
)


# ---- Pump Design Conditions ----
FuelDesignCoefficients = TurboDesignCoefficients(
    "Fuel Pump Design Coefficients",
    PumpNetwork,
    rotor_speed=23000,
    volumetric_flow=2.7 / 1000,
    head_rise=323.4,
    torque=6.31,
    density=50.8 * LBM_FT3_TO_KG_M3,
    flow_geometric_parameter=1.56 * IN_TO_M,
    head_geometric_parameter=1.56 * IN_TO_M,
    torque_geometric_parameter=1.56 * IN_TO_M,
)

OxidizerDesignCoefficients = TurboDesignCoefficients(
    "Oxidizer Pump Design Coefficients",
    PumpNetwork,
    rotor_speed=23000,
    volumetric_flow=4.4 / 1000,
    head_rise=159.8,
    torque=4.77,
    density=1104,
    flow_geometric_parameter=1.84 * IN_TO_M,
    head_geometric_parameter=1.84 * IN_TO_M,
    torque_geometric_parameter=1.84 * IN_TO_M,
)


# ---- Pump Map Data ----
normalized_flow_coefficient_map = [
    0.00,
    0.40,
    0.70,
    1.00,
    1.15,
    1.30,
    1.45,
]

normalized_head_coefficient_map = [
    1.25,
    1.20,
    1.12,
    1.00,
    0.84,
    0.62,
    0.35,
]

normalized_torque_coefficient_map = [
    0.15,
    0.38,
    0.66,
    1.00,
    1.17,
    1.33,
    1.50,
]


# ---- Components ----
COPV = Boundary(
    "COPV",
    PumpNetwork,
    PressurantSupply.pressure,
    PressurantSupply.temperature,
)

FuelBangBang = IsentropicGasRegulator(
    "FBB",
    PumpNetwork,
    COPV.pressure,
    COPV.temperature,
    set_pressure=FuelUllageGas.pressure,
    discharge_coefficient=1,
    cross_sectional_area=np.pi / 4 * (0.25 * IN_TO_M)**2,
    specific_gas_constant=FuelUllageGas.gas_constant,
    specific_heat_ratio=FuelUllageGas.specific_heat_ratio,
)

OxBangBang = IsentropicGasRegulator(
    "OBB",
    PumpNetwork,
    COPV.pressure,
    COPV.temperature,
    set_pressure=OxUllageGas.pressure,
    discharge_coefficient=1,
    cross_sectional_area=np.pi / 4 * (0.25 * IN_TO_M)**2,
    specific_gas_constant=OxUllageGas.gas_constant,
    specific_heat_ratio=OxUllageGas.specific_heat_ratio,
)

FuelTank = PressurizedTank(
    "Fuel Tank",
    PumpNetwork,
    pressure=FuelTankFluid.pressure,
    pressurant_density=FuelUllageGas.density,
    liquid_density=FuelTankFluid.density,
    collapse_factor=1,
    ullage_temperature=FuelUllageGas.temperature,
    liquid_temperature=FuelTankFluid.temperature,
    mass_flow_in=FuelBangBang.mass_flow,
    mass_flow_out=State((2.7 / 1000) * (50.8 * LBM_FT3_TO_KG_M3)),
)

OxTank = PressurizedTank(
    "Ox Tank",
    PumpNetwork,
    pressure=OxTankFluid.pressure,
    pressurant_density=OxUllageGas.density,
    liquid_density=OxTankFluid.density,
    collapse_factor=1.4,
    ullage_temperature=OxUllageGas.temperature,
    liquid_temperature=OxTankFluid.temperature,
    mass_flow_in=OxBangBang.mass_flow,
    mass_flow_out=State(2.2 * (2.7 / 1000) * (50.8 * LBM_FT3_TO_KG_M3)),
)

OxGravitydP = GravityPressureChange(
    "Ox Runline Gravity dP",
    PumpNetwork,
    upstream_pressure=OxTank.pressure,
    density=OxTank.liquid_density,
    elevation_change=-0.5,
    mass_flow=OxTank.mass_flow_out
)

FuelRunline = DarcyWeisbach(
    "Fuel Main Line",
    PumpNetwork,
    upstream_pressure=State(),
    downstream_pressure=FuelPumpInletFluid.pressure,
    length=0.5,
    cross_sectional_area=np.pi / 4 * (0.5 * IN_TO_M)**2,
    hydraulic_diameter=0.5 * IN_TO_M,
    density=FuelTank.liquid_density,
    mass_flow=FuelTank.mass_flow_out,
    friction_factor=State(0.02),
)


FuelGravitydP = GravityPressureChange(
    "Fuel Runline Gravity dP",
    PumpNetwork,
    upstream_pressure=FuelTank.pressure,
    density=FuelTank.liquid_density,
    elevation_change=-0.5,
    mass_flow=FuelTank.mass_flow_out,
    downstream_pressure=FuelRunline.upstream_pressure
)

FuelRunlineFriction = Colebrook(
    "Fuel Main Line Friction",
    PumpNetwork,
    mass_flow=FuelRunline.mass_flow,
    hydraulic_diameter=FuelRunline.hydraulic_diameter,
    dynamic_viscosity=FuelTankFluid.dynamic_viscosity,
    cross_sectional_area=FuelRunline.cross_sectional_area,
    roughness=1.5e-6,
    friction_factor=FuelRunline.friction_factor,
)


OxRunline = DarcyWeisbach(
    "Ox Main Line",
    PumpNetwork,
    upstream_pressure=OxGravitydP.downstream_pressure,
    downstream_pressure=OxPumpInletFluid.pressure,
    length=0.5,
    cross_sectional_area=np.pi / 4 * (0.5 * IN_TO_M)**2,
    hydraulic_diameter=0.5 * IN_TO_M,
    density=OxTank.liquid_density,
    mass_flow=OxTank.mass_flow_out,
    friction_factor=State(0.02),
)

OxRunlineFriction = Colebrook(
    "Ox Main Line Friction",
    PumpNetwork,
    mass_flow=OxRunline.mass_flow,
    hydraulic_diameter=OxRunline.hydraulic_diameter,
    dynamic_viscosity=OxTankFluid.dynamic_viscosity,
    cross_sectional_area=OxRunline.cross_sectional_area,
    roughness=1.5e-6,
    friction_factor=OxRunline.friction_factor,
)


FuelPumpInlet = IsothermalIncompressibleVolume(
    "Fuel Suction Inlet",
    PumpNetwork,
    pressure=FuelPumpInletFluid.pressure,
    temperature=FuelPumpInletFluid.temperature,
    density=FuelPumpInletFluid.density,
    enthalpy=FuelPumpInletFluid.enthalpy,
    mass_flow_in=FuelRunline.mass_flow,
    mass_flow_out=State((2.7 / 1000) * (50.8 * LBM_FT3_TO_KG_M3)),
    volume=0.2 * IN3_TO_M3,
)

OxPumpInlet = IsothermalIncompressibleVolume(
    "Ox Suction Inlet",
    PumpNetwork,
    pressure=OxPumpInletFluid.pressure,
    temperature=OxPumpInletFluid.temperature,
    density=OxPumpInletFluid.density,
    enthalpy=OxPumpInletFluid.enthalpy,
    mass_flow_in=OxRunline.mass_flow,
    mass_flow_out=State(2.2 * (2.7 / 1000) * (50.8 * LBM_FT3_TO_KG_M3)),
    volume=0.2 * IN3_TO_M3,
)

FuelEPumpMap = TurboMap(
    "Fuel E-Pump Map",
    PumpNetwork,
    rotor_speed=fuel_shaft_speed,
    volumetric_flow=FuelPumpInlet.mass_flow_out / FuelPumpDischargeFluid.density,
    density=FuelPumpInletFluid.density,
    flow_geometric_parameter=1.56 * IN_TO_M,
    head_geometric_parameter=1.56 * IN_TO_M,
    torque_geometric_parameter=1.56 * IN_TO_M,
    design_flow_coefficient=FuelDesignCoefficients.flow_coefficient,
    design_head_coefficient=FuelDesignCoefficients.head_coefficient,
    design_torque_coefficient=FuelDesignCoefficients.torque_coefficient,
    normalized_flow_coefficient_map=normalized_flow_coefficient_map,
    normalized_head_coefficient_map=normalized_head_coefficient_map,
    normalized_torque_coefficient_map=normalized_torque_coefficient_map,
)

OxEPumpMap = TurboMap(
    "Ox E-Pump Map",
    PumpNetwork,
    rotor_speed=ox_shaft_speed,
    volumetric_flow=OxPumpInlet.mass_flow_out / OxPumpDischargeFluid.density,
    density=OxPumpInletFluid.density,
    flow_geometric_parameter=1.84 * IN_TO_M,
    head_geometric_parameter=1.84 * IN_TO_M,
    torque_geometric_parameter=1.84 * IN_TO_M,
    design_flow_coefficient=OxidizerDesignCoefficients.flow_coefficient,
    design_head_coefficient=OxidizerDesignCoefficients.head_coefficient,
    design_torque_coefficient=OxidizerDesignCoefficients.torque_coefficient,
    normalized_flow_coefficient_map=normalized_flow_coefficient_map,
    normalized_head_coefficient_map=normalized_head_coefficient_map,
    normalized_torque_coefficient_map=normalized_torque_coefficient_map,
)

FuelEPump = PolytropicPump(
    "Fuel Electric Pump",
    PumpNetwork,
    rotor_speed=fuel_shaft_speed,
    head_rise=FuelEPumpMap.head_rise,
    mass_flow=FuelPumpInlet.mass_flow_out,
    torque=FuelEPumpMap.torque,
    upstream_density=FuelPumpInlet.density,
    downstream_density=FuelPumpDischargeFluid.density,
    upstream_total_pressure=FuelPumpInlet.pressure,
    discharge_total_pressure=FuelPumpDischargeFluid.pressure,
    upstream_total_enthalpy=FuelPumpInlet.enthalpy,
)

OxEPump = PolytropicPump(
    "Ox Electric Pump",
    PumpNetwork,
    rotor_speed=ox_shaft_speed,
    head_rise=OxEPumpMap.head_rise,
    mass_flow=OxPumpInlet.mass_flow_out,
    torque=OxEPumpMap.torque,
    upstream_density=OxPumpInlet.density,
    downstream_density=OxPumpDischargeFluid.density,
    upstream_total_pressure=OxPumpInlet.pressure,
    discharge_total_pressure=OxPumpDischargeFluid.pressure,
    upstream_total_enthalpy=OxPumpInlet.enthalpy,
)

FuelPumpOutlet = SimpleVolume(
    "Fuel Discharge Outlet",
    PumpNetwork,
    pressure=FuelPumpDischargeFluid.pressure,
    enthalpy=FuelPumpDischargeFluid.enthalpy,
    volume=0.2 * IN3_TO_M3,
    temperature=FuelPumpDischargeFluid.temperature,
    density=FuelPumpDischargeFluid.density,
    mass_flow_in=FuelEPump.mass_flow,
    mass_flow_out=State((2.7 / 1000) * (50.8 * LBM_FT3_TO_KG_M3)),
    enthalpy_in=FuelEPump.discharge_total_enthalpy,
)


OxPumpOutlet = SimpleVolume(
    "Ox Discharge Outlet",
    PumpNetwork,
    pressure=OxPumpDischargeFluid.pressure,
    enthalpy=OxPumpDischargeFluid.enthalpy,
    volume=0.2 * IN3_TO_M3,
    temperature=OxPumpDischargeFluid.temperature,
    density=OxPumpDischargeFluid.density,
    mass_flow_in=OxEPump.mass_flow,
    mass_flow_out=State(2.2 * (2.7 / 1000) * (50.8 * LBM_FT3_TO_KG_M3)),
    enthalpy_in=OxEPump.discharge_total_enthalpy,
)

'''
FuelInjectorInletLine = DischargeCoefficient(
    "Fuel Injector Inlet",
    PumpNetwork,
    upstream_pressure=FuelPumpOutlet.pressure,
    downstream_pressure=FuelManifoldFluid.pressure,
    density=FuelPumpOutlet.density,
    discharge_coefficient=1,
    cross_sectional_area=5.0 * 0.56e-4,
    mass_flow=FuelPumpOutlet.mass_flow_out,
)

OxInjectorInletLine = DischargeCoefficient(
    "Ox Injector Inlet",
    PumpNetwork,
    upstream_pressure=OxPumpOutlet.pressure,
    downstream_pressure=OxManifoldFluid.pressure,
    density=OxPumpOutlet.density,
    discharge_coefficient=1,
    cross_sectional_area=5.0 * 1.25e-4,
    mass_flow=OxPumpOutlet.mass_flow_out,
)
'''

FuelInjectorInletLine = DarcyWeisbach(
    "Fuel Injector Inlet",
    PumpNetwork,
    upstream_pressure=FuelPumpOutlet.pressure,
    downstream_pressure=FuelManifoldFluid.pressure,
    length=0.25,
    cross_sectional_area=np.pi / 4 * (0.5 * IN_TO_M)**2,
    hydraulic_diameter=0.5 * IN_TO_M,
    density=FuelPumpOutlet.density,
    mass_flow=FuelPumpOutlet.mass_flow_out,
    friction_factor=State(0.02),
)

FuelInjectorInletFriction = Colebrook(
    "Fuel Injector Inlet Friction",
    PumpNetwork,
    mass_flow=FuelInjectorInletLine.mass_flow,
    hydraulic_diameter=FuelInjectorInletLine.hydraulic_diameter,
    dynamic_viscosity=FuelPumpDischargeFluid.dynamic_viscosity,
    cross_sectional_area=FuelInjectorInletLine.cross_sectional_area,
    roughness=1.5e-6,
    friction_factor=FuelInjectorInletLine.friction_factor,
)


OxInjectorInletLine = DarcyWeisbach(
    "Ox Injector Inlet",
    PumpNetwork,
    upstream_pressure=OxPumpOutlet.pressure,
    downstream_pressure=OxManifoldFluid.pressure,
    length=0.25,
    cross_sectional_area=np.pi / 4 * (0.5 * IN_TO_M)**2,
    hydraulic_diameter=0.5 * IN_TO_M,
    density=OxPumpOutlet.density,
    mass_flow=OxPumpOutlet.mass_flow_out,
    friction_factor=State(0.02),
)

OxInjectorInletFriction = Colebrook(
    "Ox Injector Inlet Friction",
    PumpNetwork,
    mass_flow=OxInjectorInletLine.mass_flow,
    hydraulic_diameter=OxInjectorInletLine.hydraulic_diameter,
    dynamic_viscosity=OxPumpDischargeFluid.dynamic_viscosity,
    cross_sectional_area=OxInjectorInletLine.cross_sectional_area,
    roughness=1.5e-6,
    friction_factor=OxInjectorInletLine.friction_factor,
)


FuelManifold = SimpleVolume(
    "Fuel Injector Manifold",
    PumpNetwork,
    pressure=FuelManifoldFluid.pressure,
    enthalpy=FuelManifoldFluid.enthalpy,
    volume=685 * IN3_TO_M3,
    temperature=FuelManifoldFluid.temperature,
    density=FuelManifoldFluid.density,
    mass_flow_in=FuelInjectorInletLine.mass_flow,
    enthalpy_in=FuelPumpOutlet.enthalpy,
)

OxManifold = SimpleVolume(
    "Ox Injector Manifold",
    PumpNetwork,
    pressure=OxManifoldFluid.pressure,
    enthalpy=OxManifoldFluid.enthalpy,
    volume=685 * IN3_TO_M3,
    temperature=OxManifoldFluid.temperature,
    density=OxManifoldFluid.density,
    mass_flow_in=OxInjectorInletLine.mass_flow,
    enthalpy_in=OxPumpOutlet.enthalpy,
)

FuelInjector = DischargeCoefficient(
    "Fuel Injector Orifices",
    PumpNetwork,
    upstream_pressure=FuelManifold.pressure,
    downstream_pressure=chamber_pressure,
    density=FuelManifold.density,
    discharge_coefficient=1,
    cross_sectional_area=0.56e-4,
    mass_flow=FuelManifold.mass_flow_out,
)

OxInjector = DischargeCoefficient(
    "Ox Injector Orifices",
    PumpNetwork,
    upstream_pressure=OxManifold.pressure,
    downstream_pressure=chamber_pressure,
    density=OxManifold.density,
    discharge_coefficient=1,
    cross_sectional_area=1.25e-4,
    mass_flow=OxManifold.mass_flow_out,
)

MainChamber = MainCombustionChamber(
    "MCC",
    PumpNetwork,
    chamber_pressure=chamber_pressure,
    oxidizer_mass_flow=OxInjector.mass_flow,
    fuel_mass_flow=FuelInjector.mass_flow,
)

Nozzle = RocketCEAChokedNozzle(
    "Nozzle",
    PumpNetwork,
    fuel=fuel,
    oxidizer=oxidizer,
    chamber_pressure=MainChamber.chamber_pressure,
    throat_area=10 * IN2_TO_M2,
    expansion_ratio=4,
    mixture_ratio=MainChamber.oxidizer_mass_flow / MainChamber.fuel_mass_flow,
    ambient_pressure=atmospheric_pressure,
    characterstic_velocity_efficiency=1.0,
    thrust_coefficient_efficiency=1.0,
    mass_flow=MainChamber.nozzle_mass_flow,
)


ChamberPressureBalance = Balance(
    "Chamber Pressure Balance",
    PumpNetwork,
    variable=FuelBangBang.discharge_coefficient,
    function=(chamber_pressure - 285 * PSIA_TO_PA),
)

MixtureRatioBalance = Balance(
    "Mixture Ratio Balance",
    PumpNetwork,
    variable=OxBangBang.discharge_coefficient,
    function=(OxInjector.mass_flow / FuelInjector.mass_flow - 2.2),
)

FuelPumpShaftSpeedBalance = Balance(
    "FPump Shaft Speed Balance",
    PumpNetwork,
    variable=FuelEPump.rotor_speed,
    function=FuelTank.pressure - 80*PSIA_TO_PA
)

OxPumpShaftSpeedBalance = Balance(
    "OPump Shaft Speed Balance",
    PumpNetwork,
    variable=OxEPump.rotor_speed,
    function=OxTank.pressure - 80*PSIA_TO_PA
)


solution = SteadyState(PumpNetwork).solve(
    return_type="dataframe",
    verbose=True,
    static=False,
)

print(solution.to_string(index=False))


print(f"Ox Tank:{OxTank.pressure.value * PA_TO_PSIA:.3f} psia")
print(f"Fuel Tank:{FuelTank.pressure.value * PA_TO_PSIA:.3f} psia")