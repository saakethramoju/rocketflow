from System import *
from Solvers import *
from constants import *


HeatExchanger = Network("Heat Exchanger")

L = 10 * IN_TO_M
Dw = 0.75 * IN_TO_M
Dox = 1.5 * IN_TO_M

Aw = (np.pi/4) * Dw**2
Aox = (np.pi/4) * Dox**2 - Aw

hot_fluid = "rp-1"
coolant = "water"

metal = 'c101'




# ---- Fluids ----

LiquidSource = FluidLookup(
    "Liquid Source",
    HeatExchanger,
    hot_fluid,
    pressure=6e5,
    temperature=500
)


LiquidNode1Fluid = FluidLookup(
    "Liquid Node 1 Fluid",
    HeatExchanger,
    hot_fluid,
    pressure=5e5,
    temperature=400,
    flash_values=("pressure", "enthalpy")
)

LiquidNode2Fluid = FluidLookup(
    "Liquid Node 2 Fluid",
    HeatExchanger,
    hot_fluid,
    pressure=4e5,
    temperature=300,
    flash_values=("pressure", "enthalpy")
)

LiquidDrain = FluidLookup(
    "Liquid Drain",
    HeatExchanger,
    hot_fluid,
    pressure=3e5,
    temperature=250  
)



CoolantSource = FluidLookup(
    "Coolant Source",
    HeatExchanger,
    coolant,
    pressure=6e5,
    quality=0
)


CoolantNode1Fluid = FluidLookup(
    "Coolant Node 1 Fluid",
    HeatExchanger,
    coolant,
    pressure=4e5,
    quality=0,
    flash_values=("pressure", "enthalpy")
)


CoolantNode2Fluid = FluidLookup(
    "Coolant Node 2 Fluid",
    HeatExchanger,
    coolant,
    pressure=3e5,
    quality=0,
    flash_values=("pressure", "enthalpy")
)

CoolantDrain = FluidLookup(
    "Coolant Drain",
    HeatExchanger,
    coolant,
    pressure=2e5,
    quality=0
)





# ---- Metals ----
SolidNode1Material = MaterialLookup(
    "Solid Node 1 Material",
    HeatExchanger,
    metal,
)

SolidNode2Material = MaterialLookup(
    "Solid Node 2 Material",
    HeatExchanger,
    metal,
)







# ---- Hot fluid side ----

LiquidTube1 = DarcyWeisbach(
    "Liquid Tube 1",
    HeatExchanger,
    mass_flow=20,
    upstream_pressure=LiquidSource.pressure,
    downstream_pressure=LiquidNode1Fluid.pressure,
    length=L,
    cross_sectional_area=Aw,
    hydraulic_diameter=Dw,
    density=LiquidSource.density,
    friction_factor=2e-5,
)


LiquidTube1Friction = Colebrook(
    "Liquid Tube 1 Friction",
    HeatExchanger,
    mass_flow=LiquidTube1.mass_flow,
    friction_factor=LiquidTube1.friction_factor,
    hydraulic_diameter=LiquidTube1.hydraulic_diameter,
    dynamic_viscosity=LiquidSource.dynamic_viscosity,
    cross_sectional_area=LiquidTube1.cross_sectional_area,
    roughness=1e-9,
)



LiquidNode1 = Volume(
    "Liquid Node 1",
    HeatExchanger,
    pressure=LiquidNode1Fluid.pressure,
    enthalpy=LiquidNode1Fluid.enthalpy,
    volume=1,
    total_enthalpy_in=LiquidSource.enthalpy,
    mass_flow_in=LiquidTube1.mass_flow,
    mass_flow_out=20
)




LiquidTube2 = DarcyWeisbach(
    "Liquid Tube 2",
    HeatExchanger,
    mass_flow=LiquidNode1.mass_flow_out,
    upstream_pressure=LiquidNode1.pressure,
    downstream_pressure=LiquidNode2Fluid.pressure,
    length=L,
    cross_sectional_area=Aw,
    hydraulic_diameter=Dw,
    density=LiquidNode1Fluid.density,
    friction_factor=2e-5,
)


LiquidTube2Friction = Colebrook(
    "Liquid Tube 2 Friction",
    HeatExchanger,
    mass_flow=LiquidTube2.mass_flow,
    friction_factor=LiquidTube2.friction_factor,
    hydraulic_diameter=LiquidTube2.hydraulic_diameter,
    dynamic_viscosity=LiquidNode1Fluid.dynamic_viscosity,
    cross_sectional_area=LiquidTube1.cross_sectional_area,
    roughness=1e-9,
)


LiquidNode2 = Volume(
    "Liquid Node 2",
    HeatExchanger,
    pressure=LiquidNode2Fluid.pressure,
    enthalpy=LiquidNode2Fluid.enthalpy,
    volume=1,
    total_enthalpy_in=LiquidNode1Fluid.enthalpy,
    mass_flow_in=LiquidTube2.mass_flow,
    mass_flow_out=20
)



LiquidTube3 = DarcyWeisbach(
    "Liquid Tube 3",
    HeatExchanger,
    mass_flow=LiquidNode2.mass_flow_out,
    upstream_pressure=LiquidNode2.pressure,
    downstream_pressure=LiquidDrain.pressure,
    length=L,
    cross_sectional_area=Aw,
    hydraulic_diameter=Dw,
    density=LiquidNode2Fluid.density,
    friction_factor=2e-5,
)


LiquidTube3Friction = Colebrook(
    "Liquid Tube 3 Friction",
    HeatExchanger,
    mass_flow=LiquidTube3.mass_flow,
    friction_factor=LiquidTube3.friction_factor,
    hydraulic_diameter=LiquidTube3.hydraulic_diameter,
    dynamic_viscosity=LiquidNode2Fluid.dynamic_viscosity,
    cross_sectional_area=LiquidTube3.cross_sectional_area,
    roughness=1e-9,
)










# ---- Coolant side ----


AnnulusPoiseuille = CircularAnnulus(
    "Circular Annulus Poiseuille",
    HeatExchanger,
    inner_diameter=Dw,
    outer_diameter=Dox,
)


AnnulusHydraulicDiameter = HydraulicDiameter(
    "Annulus Hydraulic Diameter",
    HeatExchanger,
    cross_sectional_area=Aox,
    wetted_perimeter=np.pi * (Dox + Dw),
)


CoolantTube1 = DarcyWeisbach(
    "Coolant Tube 1",
    HeatExchanger,
    mass_flow=10,
    upstream_pressure=CoolantSource.pressure,
    downstream_pressure=CoolantNode1Fluid.pressure,
    length=L,
    cross_sectional_area=Aox,
    hydraulic_diameter=AnnulusHydraulicDiameter.hydraulic_diameter,
    density=CoolantSource.density,
    friction_factor=2e-5,
)


CoolantTube1Friction = Colebrook(
    "Coolant Tube 1 Friction",
    HeatExchanger,
    mass_flow=CoolantTube1.mass_flow,
    friction_factor=CoolantTube1.friction_factor,
    hydraulic_diameter=CoolantTube1.hydraulic_diameter,
    dynamic_viscosity=CoolantSource.dynamic_viscosity,
    cross_sectional_area=CoolantTube1.cross_sectional_area,
    roughness=1e-9,
    poiseuille_number=AnnulusPoiseuille.poiseuille_number
)


CoolantNode1 = Volume(
    "Coolant Node 1",
    HeatExchanger,
    pressure=CoolantNode1Fluid.pressure,
    enthalpy=CoolantNode1Fluid.enthalpy,
    volume=1,
    total_enthalpy_in=CoolantSource.enthalpy,
    mass_flow_in=CoolantTube1.mass_flow,
    mass_flow_out=10
)


CoolantTube2 = DarcyWeisbach(
    "Coolant Tube 2",
    HeatExchanger,
    mass_flow=CoolantNode1.mass_flow_out,
    upstream_pressure=CoolantNode1.pressure,
    downstream_pressure=CoolantNode2Fluid.pressure,
    length=L,
    cross_sectional_area=Aox,
    hydraulic_diameter=AnnulusHydraulicDiameter.hydraulic_diameter,
    density=CoolantNode1Fluid.density,
    friction_factor=2e-5,
)


CoolantTube2Friction = Colebrook(
    "Coolant Tube 2 Friction",
    HeatExchanger,
    mass_flow=CoolantTube2.mass_flow,
    friction_factor=CoolantTube2.friction_factor,
    hydraulic_diameter=CoolantTube2.hydraulic_diameter,
    dynamic_viscosity=CoolantNode1Fluid.dynamic_viscosity,
    cross_sectional_area=CoolantTube2.cross_sectional_area,
    roughness=1e-9,
    poiseuille_number=AnnulusPoiseuille.poiseuille_number
)


CoolantNode2 = Volume(
    "Coolant Node 2",
    HeatExchanger,
    pressure=CoolantNode2Fluid.pressure,
    enthalpy=CoolantNode2Fluid.enthalpy,
    volume=1,
    total_enthalpy_in=CoolantNode1Fluid.enthalpy,
    mass_flow_in=CoolantTube2.mass_flow,
    mass_flow_out=10
)



CoolantTube3 = DarcyWeisbach(
    "Coolant Tube 3",
    HeatExchanger,
    mass_flow=CoolantNode2.mass_flow_out,
    upstream_pressure=CoolantNode2.pressure,
    downstream_pressure=CoolantDrain.pressure,
    length=L,
    cross_sectional_area=Aox,
    hydraulic_diameter=AnnulusHydraulicDiameter.hydraulic_diameter,
    density=CoolantNode2Fluid.density,
    friction_factor=2e-5,
)


CoolantTube3Friction = Colebrook(
    "Coolant Tube 3 Friction",
    HeatExchanger,
    mass_flow=CoolantTube3.mass_flow,
    friction_factor=CoolantTube3.friction_factor,
    hydraulic_diameter=CoolantTube3.hydraulic_diameter,
    dynamic_viscosity=CoolantNode2Fluid.dynamic_viscosity,
    cross_sectional_area=CoolantTube3.cross_sectional_area,
    roughness=1e-9,
    poiseuille_number=AnnulusPoiseuille.poiseuille_number
)








# ---- Heat tranfser ----






# ---- Solver ----

solution = SteadyState(HeatExchanger).solve(
    verbose=True,
    print_solution=True,
)





'''
ThermalSystem = Network("Thermal")

D = 0.5 * IN_TO_M
L = 10.0 * IN_TO_M

flow_area = (np.pi / 4.0) * D**2
wetted_area = np.pi * D * L

TubeMaterial = MaterialLookup(
    "Tube Material",
    ThermalSystem,
    "718",
)

SourceFluid = FluidLookup(
    "Fluid Source",
    ThermalSystem,
    "lox",
    pressure=5e5,
    quality=0,
)

VolFluid = FluidLookup(
    "Fluid Volume",
    ThermalSystem,
    "lox",
    pressure=1e5,
    quality=0,
)

TubeFlow = DarcyWeisbach(
    "Tube Flow",
    ThermalSystem,
    mass_flow=5,
    upstream_pressure=SourceFluid.pressure,
    downstream_pressure=VolFluid.pressure,
    length=L,
    cross_sectional_area=flow_area,
    hydraulic_diameter=D,
    density=SourceFluid.density,
    friction_factor=2e-5,
)

TubeFriction = Colebrook(
    "Tube Friction",
    ThermalSystem,
    mass_flow=TubeFlow.mass_flow,
    friction_factor=TubeFlow.friction_factor,
    hydraulic_diameter=TubeFlow.hydraulic_diameter,
    dynamic_viscosity=SourceFluid.dynamic_viscosity,
    cross_sectional_area=TubeFlow.cross_sectional_area,
    roughness=1e-9,
)

Tube = Solid(
    "Tube",
    ThermalSystem,
    temperature=TubeMaterial.temperature,
)

LOXConvection = Convection(
    "LOX Convection",
    ThermalSystem,
    surface_temperature=Tube.temperature,
    fluid_temperature=SourceFluid.temperature,
    convective_area=wetted_area,
    convection_coefficient=25,
)


LOXGnielinski = Gnielinski(
    "LOX Gnielinski",
    ThermalSystem,
    hydraulic_diameter=TubeFlow.hydraulic_diameter,
    friction_factor=TubeFriction.friction_factor,
    fluid_conductivity=SourceFluid.conductivity,
    fluid_specific_heat=SourceFluid.specific_heat_cp,
    fluid_dynamic_viscosity=SourceFluid.dynamic_viscosity,
    cross_sectional_area=TubeFlow.cross_sectional_area,
    mass_flow=TubeFlow.mass_flow,
    convection_coefficient=LOXConvection.convection_coefficient,  
)


AirConvection = Convection(
    "Air Convection",
    ThermalSystem,
    surface_temperature=Tube.temperature,
    fluid_temperature=293.15,
    convective_area=wetted_area,
    convection_coefficient=25,
)

Tube.heat_rate = LOXConvection.heat_rate + AirConvection.heat_rate

solution = SteadyState(ThermalSystem).solve(
    verbose=True,
    print_solution=True,
)
'''