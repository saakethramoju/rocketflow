from System import *
from Solvers import *
from constants import *
from thermoprop import Propellant

PropFeedSystem = Network("Prop Feed System")

#print(FluidRegistry.propellant_supported_names)

SourceProp = PropellantLookup(
    "Source Prop",
    PropFeedSystem,
    "rp-1",
    temperature=90,
    pressure=5e5
)


VolumeProp = PropellantLookup(
    "Volume Prop",
    PropFeedSystem,
    SourceProp.propellant_name,
    temperature=90,
    pressure=1e5
)


TubeIn = DischargeCoefficient(
    "Tube In",
    PropFeedSystem,
    upstream_pressure=SourceProp.pressure,
    downstream_pressure=VolumeProp.pressure,
    density=SourceProp.density,
    discharge_coefficient=1,
    cross_sectional_area=0.55e-4
)


Vol = SimpleVolume(
    "Vol",
    PropFeedSystem,
    pressure=VolumeProp.pressure,
    volume=1,
    mass_flow_in=TubeIn.mass_flow
)


TubeOut = DischargeCoefficient(
    "Tube In",
    PropFeedSystem,
    upstream_pressure=VolumeProp.pressure,
    downstream_pressure=101325,
    density=VolumeProp.density,
    discharge_coefficient=1,
    cross_sectional_area=0.55e-4,
    mass_flow=Vol.mass_flow_out
)

SteadyState(PropFeedSystem).solve(
    #model="Main Pump",
    verbose=True,
    print_solution=True,
    filename="test.xlsx",
    #evaluate_all_model_options=True,
)




# ------------------------------------------------------------------
# RocketProps state model
# ------------------------------------------------------------------
# Propellant is not a full thermodynamic EOS package like Fluid or
# IdealGas. Most RocketProps properties are functions of temperature
# alone and are evaluated from saturated-liquid correlations.
#
# Temperature-only properties typically include:
#     density
#     specific_gravity
#     viscosity
#     kinematic_viscosity
#     specific_heat
#     thermal_conductivity
#     prandtl
#     surface_tension
#     vapor_pressure
#     saturation_pressure
#     critical_pressure
#     critical_temperature
#     molecular_weight
#
# If pressure is not provided, the propellant is assumed to be a
# saturated liquid at the specified temperature and:
#
#     pressure = None
#
# The saturation pressure can still be obtained from:
#
#     vapor_pressure
#     saturation_pressure
#
# If pressure is provided, RocketProps can apply compressed-liquid
# corrections for properties that depend on pressure. The pressure
# must be greater than the saturation pressure at that temperature
#
# RocketProps does not provide a complete thermodynamic state and
# therefore does not support quantities such as:
#
#     enthalpy
#     internal_energy
#     entropy
#     gibbs_energy
#     free_energy
#     quality
#
# Those properties should instead be obtained from Fluid or IdealGas.
# ------------------------------------------------------------------