from System import *
from Solvers import *

from constants import *

# --- Network Definition ---
PumpNetwork = Network("Pumped System")

# ---- Fluids -----
fuel = "rp-1"

SuctionFluid = FluidLookup(
    "Inlet Fluid",
    PumpNetwork,
    fuel,
    pressure=50 * PSIA_TO_PA,
    temperature=293.15,
)

DischargeFluid = FluidLookup(
    "Outlet Fluid",
    PumpNetwork,
    fuel,
    pressure=50 * PSIA_TO_PA,
    temperature=293.15,
)


# ----- Components ------
Source = PressureBoundary("Source", PumpNetwork, pressure=80*PSIA_TO_PA)

InletLine = DischargeCoefficient("Inlet Line", PumpNetwork, upstream_pressure=Source.pressure,
                                 downstream_pressure=SuctionFluid.pressure, density=SuctionFluid.density,
                                 discharge_coefficient=1, cross_sectional_area=1.0e-4)


Inlet = IsothermalIncompressibleVolume("Suction Inlet", PumpNetwork, pressure=SuctionFluid.pressure,
                                       temperature=SuctionFluid.temperature, mass_flow_in=InletLine.mass_flow, volume=0.2*IN3_TO_M3)
'''
FuelPumpMap = SimpleEulerCentrifugalPump(
    name="Fuel Pump Map Generator",
    network=PumpNetwork,

    rotor_speed=20000,
    volumetric_flow=0.0075,
    density=SuctionFluid.density,

    impeller_inlet_tip_radius=0.025,
    impeller_outlet_tip_radius=0.040,

    inlet_annular_flow_area=0.0025,
    outlet_annular_flow_area=0.0018,

    inlet_blade_angle=18.0,
    outlet_blade_angle=35.0,
    angle_units="degrees",

    slip_factor=0.85,

    hydraulic_efficiency=0.65,
    mechanical_efficiency=0.90,
    volumetric_efficiency=0.97,
)

FuelEPump = ConstantDensityPump(
    "Fuel E-Pump",
    PumpNetwork,
    rotor_speed=,
    head_rise=FuelPumpMap.head_rise,
    volumetric_flow=FuelPumpMap.volumetric_flow,
    density=SuctionFluid.density,
    torque=FuelPumpMap.torque,
    upstream_total_pressure=Inlet.pressure,
    discharge_total_pressure=DischargeFluid.pressure,
    mass_flow=Inlet.mass_flow_out,
)
'''

Q_design = 0.0018   # m^3/s
H_design = 47.0     # m, not 85
N_design = 20000.0  # rpm
eta_design = 0.65


def pump_head(Q, N):
    speed_ratio = N / N_design
    phi = Q / (Q_design * speed_ratio)

    H_shutoff = 1.35 * H_design
    H = H_shutoff * speed_ratio**2 * (1.0 - 0.259 * phi**2)

    return H.clip(lower=0.0)

def pump_torque(Q, N, rho):
    H = pump_head(Q, N)
    omega = np.pi / 30.0 * N

    hydraulic_power = rho * 9.80665 * Q * H
    shaft_power = hydraulic_power / eta_design

    return shaft_power / omega


shaft_speed = State(N_design)
#flowrate = State(Q_design, bounds=(0.0, 0.003), keep_feasible=True)
flowrate = State(Q_design)


FuelEPump = ConstantDensityPump(
    "Fuel E-Pump",
    PumpNetwork,
    rotor_speed=shaft_speed,
    head_rise=pump_head(flowrate, shaft_speed),
    volumetric_flow=flowrate,
    density=SuctionFluid.density,
    torque=pump_torque(flowrate, shaft_speed, SuctionFluid.density),
    upstream_total_pressure=Inlet.pressure,
    discharge_total_pressure=DischargeFluid.pressure,
    mass_flow=Inlet.mass_flow_out,
)

Outlet = IsothermalIncompressibleVolume("Discharge", PumpNetwork, pressure=DischargeFluid.pressure, temperature=DischargeFluid.temperature,
                                        volume=0.2*IN3_TO_M3, density=DischargeFluid.density, mass_flow_in=FuelEPump.mass_flow)

OutletLine = GenericDarcyWeisbach("Outlet Line", PumpNetwork, upstream_pressure=Outlet.pressure, 
                                   downstream_pressure=101325, length=1.5, cross_sectional_area=0.2*IN2_TO_M2,
                                   hydraulic_diameter=0.5*IN_TO_M, roughness=0.1e-3, density=Outlet.density, 
                                   dynamic_viscosity=DischargeFluid.dynamic_viscosity, mass_flow=Outlet.mass_flow_out)

solution = SteadyState(PumpNetwork).solve(return_type='dataframe', verbose=True, static=False)
print(solution.to_string(index=False))






# ----- Pump Summary -----
Q = FuelEPump.volumetric_flow.value
rho = FuelEPump.density.value
mdot = FuelEPump.mass_flow.value

N = FuelEPump.rotor_speed.value
omega = np.pi / 30.0 * N

H = FuelEPump.head_rise.value
T = FuelEPump.torque.value
eta = FuelEPump.efficiency.value
shaft_power = FuelEPump.shaft_power.value

po_in = FuelEPump.upstream_total_pressure.value
po_out = FuelEPump.discharge_total_pressure.value
dp = po_out - po_in

print("\n" + "=" * 60)
print("                 CONSTANT DENSITY PUMP SUMMARY")
print("=" * 60)

print(f"{'Fluid':35s}: {fuel}")
print(f"{'Density':35s}: {rho:12.3f} kg/m^3")
print(f"{'Rotor Speed':35s}: {N:12.3f} RPM")
print(f"{'Angular Speed':35s}: {omega:12.3f} rad/s")

print("\n--- Flow ---")
print(f"{'Volumetric Flow':35s}: {Q:12.6f} m^3/s")
print(f"{'Volumetric Flow':35s}: {Q / IN3_TO_M3 * 60.0:12.3f} in^3/min")
print(f"{'Mass Flow':35s}: {mdot:12.6f} kg/s")

print("\n--- Pressure / Head ---")
print(f"{'Inlet Total Pressure':35s}: {po_in:12.3f} Pa")
print(f"{'Inlet Total Pressure':35s}: {po_in / PSIA_TO_PA:12.3f} psia")
print(f"{'Discharge Total Pressure':35s}: {po_out:12.3f} Pa")
print(f"{'Discharge Total Pressure':35s}: {po_out / PSIA_TO_PA:12.3f} psia")
print(f"{'Pressure Rise':35s}: {dp:12.3f} Pa")
print(f"{'Pressure Rise':35s}: {dp / PSIA_TO_PA:12.3f} psid")
print(f"{'Head Rise':35s}: {H:12.3f} m")

print("\n--- Power / Torque ---")
print(f"{'Torque':35s}: {T:12.6f} N-m")
print(f"{'Shaft Power':35s}: {shaft_power:12.3f} W")
print(f"{'Shaft Power':35s}: {shaft_power / 1000.0:12.3f} kW")
print(f"{'Efficiency':35s}: {eta:12.6f}")
print(f"{'Efficiency':35s}: {100.0 * eta:12.3f} %")

print("=" * 60)