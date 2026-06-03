from System import *
from Solvers import *

from constants import *


PumpNetwork = Network("Multi Model Pump Test")

D = 3 * IN_TO_M
A = (np.pi / 4) * D**2

# ------------------------------------------------------------------
# Fluids
# ------------------------------------------------------------------

InletFluid = FluidLookup(
    "Inlet Fluid",
    PumpNetwork,
    "rp-1",
    pressure=3e5,
    temperature=300,
)

OutletFluid = FluidLookup(
    "Outlet Fluid",
    PumpNetwork,
    "rp-1",
    pressure=9e5,
    temperature=300,
    flash_values=("pressure", "enthalpy"),
)

# ------------------------------------------------------------------
# Shared pump map
# ------------------------------------------------------------------

RotorSpeed = State(22000.0)
ImpellerDiameter = State(1.6 * IN_TO_M)

ConstantPumpFlow = State(2.0e-3)
ConstantPumpHeadRise = State()
ConstantPumpTorque = State()

DesignCoefficients = TurboDesignCoefficients(
    "Pump Design Coefficients",
    PumpNetwork,
    rotor_speed=22000.0,
    volumetric_flow=2.5e-3,
    head_rise=80.0,
    torque=3.5,
    density=InletFluid.density,
    flow_geometric_parameter=ImpellerDiameter,
    head_geometric_parameter=ImpellerDiameter,
    torque_geometric_parameter=ImpellerDiameter,
)

normalized_flow_coefficient_map = [0.50, 0.75, 1.00, 1.25, 1.50]
normalized_head_coefficient_map = [1.20, 1.10, 1.00, 0.80, 0.55]
normalized_torque_coefficient_map = [0.60, 0.80, 1.00, 1.20, 1.45]

PumpMap = TurboMap(
    "Pump Map",
    PumpNetwork,
    rotor_speed=RotorSpeed,
    volumetric_flow=ConstantPumpFlow,
    density=InletFluid.density,
    flow_geometric_parameter=ImpellerDiameter,
    head_geometric_parameter=ImpellerDiameter,
    torque_geometric_parameter=ImpellerDiameter,
    design_flow_coefficient=DesignCoefficients.flow_coefficient,
    design_head_coefficient=DesignCoefficients.head_coefficient,
    design_torque_coefficient=DesignCoefficients.torque_coefficient,
    normalized_flow_coefficient_map=normalized_flow_coefficient_map,
    normalized_head_coefficient_map=normalized_head_coefficient_map,
    normalized_torque_coefficient_map=normalized_torque_coefficient_map,
    head_rise=ConstantPumpHeadRise,
    torque=ConstantPumpTorque,
)

# ------------------------------------------------------------------
# Model 1: Pump model
# Both options use the same external PumpMap outputs.
# ------------------------------------------------------------------

ConstantDensityPumpOption = ConstantDensityPump.model(
    "constant_density",
    rotor_speed=RotorSpeed,
    head_rise=PumpMap.head_rise,
    volumetric_flow=ConstantPumpFlow,
    density=InletFluid.density,
    torque=PumpMap.torque,
    upstream_total_pressure=InletFluid.pressure,
    discharge_total_pressure=OutletFluid.pressure,
    upstream_total_enthalpy=InletFluid.enthalpy,
)

PolyPumpMassFlow = State(1.5)

PolytropicPumpOption = PolytropicPump.model(
    "polytropic",
    rotor_speed=RotorSpeed,
    head_rise=PumpMap.head_rise,
    mass_flow=PolyPumpMassFlow,
    upstream_density=InletFluid.density,
    downstream_density=OutletFluid.density,
    torque=PumpMap.torque,
    upstream_total_pressure=InletFluid.pressure,
    discharge_total_pressure=OutletFluid.pressure,
    upstream_total_enthalpy=InletFluid.enthalpy,
)

PumpModel = Model(
    "Main Pump",
    PumpNetwork,
    ConstantDensityPumpOption,
    PolytropicPumpOption,
    order=[
        "constant_density",
        "polytropic",
    ],
)

# ------------------------------------------------------------------
# Model 2: Outlet model
# ------------------------------------------------------------------

OutletCdA = DischargeCoefficient.model(
    "cda",
    upstream_pressure=OutletFluid.pressure,
    downstream_pressure=101325.0,
    density=OutletFluid.density,
    discharge_coefficient=1.0,
    cross_sectional_area=A / 4,
)

OutletMassFlow = State(1.0)
OutletFrictionFactor = State(0.02)

OutletDarcy = ModelOption(
    "darcy",
    components=[
        Churchill.model(
            "Outlet Friction",
            mass_flow=OutletMassFlow,
            friction_factor=OutletFrictionFactor,
            hydraulic_diameter=D,
            dynamic_viscosity=OutletFluid.dynamic_viscosity,
            cross_sectional_area=A / 4,
        ),
        DarcyWeisbach.model(
            "Outlet Darcy",
            mass_flow=OutletMassFlow,
            upstream_pressure=OutletFluid.pressure,
            downstream_pressure=101325.0,
            length=1.0,
            cross_sectional_area=A / 4,
            hydraulic_diameter=D,
            density=OutletFluid.density,
            friction_factor=OutletFrictionFactor,
        ),
    ],
)

OutletModel = Model(
    "Outlet",
    PumpNetwork,
    OutletCdA,
    OutletDarcy,
    order=[
        "cda",
        "darcy",
    ],
)

# ------------------------------------------------------------------
# Test A: normal solve
# Uses first option from each model:
#   Main Pump -> constant_density
#   Outlet    -> cda
# ------------------------------------------------------------------
"""
SteadyState(PumpNetwork).solve(
    verbose=True,
    print_solution=True,
    filename="test.xlsx",
)
"""
# ------------------------------------------------------------------
# Test B: sweep only the pump model
# Uses the same external PumpMap for both pump options.
# ------------------------------------------------------------------

SteadyState(PumpNetwork).solve(
    model="Main Pump",
    verbose=True,
    print_solution=True,
    filename="test.xlsx",
    evaluate_all_model_options=True,
)







'''
from System import *
from Solvers import *

from constants import *


PumpNetwork = Network("Multi Model Pump Test")

D = 3 * IN_TO_M
A = (np.pi / 4) * D**2

# ------------------------------------------------------------------
# Fluids
# ------------------------------------------------------------------

InletFluid = FluidLookup(
    "Inlet Fluid",
    PumpNetwork,
    "rp-1",
    pressure=3e5,
    temperature=300,
)

OutletFluid = FluidLookup(
    "Outlet Fluid",
    PumpNetwork,
    "rp-1",
    pressure=9e5,
    temperature=300,
    flash_values=("pressure", "enthalpy"),
)

# ------------------------------------------------------------------
# Pump map setup
# ------------------------------------------------------------------

RotorSpeed = State(22000.0)
ImpellerDiameter = State(1.6 * IN_TO_M)

DesignCoefficients = TurboDesignCoefficients(
    "Pump Design Coefficients",
    PumpNetwork,
    rotor_speed=22000.0,
    volumetric_flow=2.5e-3,
    head_rise=80.0,
    torque=3.5,
    density=InletFluid.density,
    flow_geometric_parameter=ImpellerDiameter,
    head_geometric_parameter=ImpellerDiameter,
    torque_geometric_parameter=ImpellerDiameter,
)

normalized_flow_coefficient_map = [0.50, 0.75, 1.00, 1.25, 1.50]
normalized_head_coefficient_map = [1.20, 1.10, 1.00, 0.80, 0.55]
normalized_torque_coefficient_map = [0.60, 0.80, 1.00, 1.20, 1.45]

# ------------------------------------------------------------------
# Model 1: Pump model
# ------------------------------------------------------------------

ConstantPumpFlow = State(2.0e-3)
ConstantPumpHeadRise = State()
ConstantPumpTorque = State()

ConstantDensityPumpOption = ModelOption(
    "constant_density",
    components=[
        TurboMap.model(
            "Constant Density Pump Map",
            rotor_speed=RotorSpeed,
            volumetric_flow=ConstantPumpFlow,
            density=InletFluid.density,
            flow_geometric_parameter=ImpellerDiameter,
            head_geometric_parameter=ImpellerDiameter,
            torque_geometric_parameter=ImpellerDiameter,
            design_flow_coefficient=DesignCoefficients.flow_coefficient,
            design_head_coefficient=DesignCoefficients.head_coefficient,
            design_torque_coefficient=DesignCoefficients.torque_coefficient,
            normalized_flow_coefficient_map=normalized_flow_coefficient_map,
            normalized_head_coefficient_map=normalized_head_coefficient_map,
            normalized_torque_coefficient_map=normalized_torque_coefficient_map,
            head_rise=ConstantPumpHeadRise,
            torque=ConstantPumpTorque,
        ),
        ConstantDensityPump.model(
            "Constant Density Pump",
            rotor_speed=RotorSpeed,
            head_rise=ConstantPumpHeadRise,
            volumetric_flow=ConstantPumpFlow,
            density=InletFluid.density,
            torque=ConstantPumpTorque,
            upstream_total_pressure=InletFluid.pressure,
            discharge_total_pressure=OutletFluid.pressure,
            upstream_total_enthalpy=InletFluid.enthalpy,
        ),
    ],
)

PolyPumpMassFlow = State(1.5)
PolyPumpHeadRise = State()
PolyPumpTorque = State()

PolytropicPumpOption = ModelOption(
    "polytropic",
    components=[
        TurboMap.model(
            "Polytropic Pump Map",
            rotor_speed=RotorSpeed,
            volumetric_flow=PolyPumpMassFlow / InletFluid.density,
            density=InletFluid.density,
            flow_geometric_parameter=ImpellerDiameter,
            head_geometric_parameter=ImpellerDiameter,
            torque_geometric_parameter=ImpellerDiameter,
            design_flow_coefficient=DesignCoefficients.flow_coefficient,
            design_head_coefficient=DesignCoefficients.head_coefficient,
            design_torque_coefficient=DesignCoefficients.torque_coefficient,
            normalized_flow_coefficient_map=normalized_flow_coefficient_map,
            normalized_head_coefficient_map=normalized_head_coefficient_map,
            normalized_torque_coefficient_map=normalized_torque_coefficient_map,
            head_rise=PolyPumpHeadRise,
            torque=PolyPumpTorque,
        ),
        PolytropicPump.model(
            "Polytropic Pump",
            rotor_speed=RotorSpeed,
            head_rise=PolyPumpHeadRise,
            mass_flow=PolyPumpMassFlow,
            upstream_density=InletFluid.density,
            downstream_density=OutletFluid.density,
            torque=PolyPumpTorque,
            upstream_total_pressure=InletFluid.pressure,
            discharge_total_pressure=OutletFluid.pressure,
            upstream_total_enthalpy=InletFluid.enthalpy,
        ),
    ],
)

PumpModel = Model(
    "Main Pump",
    PumpNetwork,
    ConstantDensityPumpOption,
    PolytropicPumpOption,
    order=[
        "constant_density",
        "polytropic",
    ],
)

# ------------------------------------------------------------------
# Model 2: Outlet model
# ------------------------------------------------------------------

OutletCdA = DischargeCoefficient.model(
    "cda",
    upstream_pressure=OutletFluid.pressure,
    downstream_pressure=101325.0,
    density=OutletFluid.density,
    discharge_coefficient=1.0,
    cross_sectional_area=A / 4,
)

OutletMassFlow = State(1.0)
OutletFrictionFactor = State(0.02)

OutletDarcy = ModelOption(
    "darcy",
    components=[
        Churchill.model(
            "Outlet Friction",
            mass_flow=OutletMassFlow,
            friction_factor=OutletFrictionFactor,
            hydraulic_diameter=D,
            dynamic_viscosity=OutletFluid.dynamic_viscosity,
            cross_sectional_area=A / 4,
        ),
        DarcyWeisbach.model(
            "Outlet Darcy",
            mass_flow=OutletMassFlow,
            upstream_pressure=OutletFluid.pressure,
            downstream_pressure=101325.0,
            length=1.0,
            cross_sectional_area=A / 4,
            hydraulic_diameter=D,
            density=OutletFluid.density,
            friction_factor=OutletFrictionFactor,
        ),
    ],
)

OutletModel = Model(
    "Outlet",
    PumpNetwork,
    OutletCdA,
    OutletDarcy,
    order=[
        "cda",
        "darcy",
    ],
)

# ------------------------------------------------------------------
# Test A: normal solve
# Uses first option from each model:
#   Main Pump -> constant_density
#   Outlet    -> cda
# ------------------------------------------------------------------
"""
solution = SteadyState(PumpNetwork).solve(
    verbose=True,
    print_solution=True,
    filename="multi_model_default.xlsx",
)
"""
# ------------------------------------------------------------------
# Test B: sweep only the pump model
# Uncomment to test all Main Pump options.
# ------------------------------------------------------------------

SteadyState(PumpNetwork).solve(
    model="Main Pump",
    verbose=True,
    print_solution=True,
    filename="test.xlsx",
    evaluate_all_model_options=True
)
'''