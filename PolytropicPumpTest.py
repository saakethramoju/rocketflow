from System import *
from Solvers import *

from constants import *

# --- Network Definition ---
PumpNetwork = Network("Pumped System")


Q_design = 0.01    # m^3/s, close to your previous converged outlet flow
H_design = 46.5      # m
N_design = 20000.0   # rpm
eta_design = 0.56745 # eta_h * eta_m * eta_v = 0.65 * 0.90 * 0.97

shaft_speed = State(30000)

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
    pressure=100 * PSIA_TO_PA,
    temperature=293.15,
    flash_values=("pressure", "enthalpy")
)


# ----- Components ------
Source = PressureBoundary("Source", PumpNetwork, pressure=1000*PSIA_TO_PA)

InletLine = DischargeCoefficient("Inlet Line", PumpNetwork, upstream_pressure=Source.pressure,
                                 downstream_pressure=SuctionFluid.pressure, density=SuctionFluid.density,
                                 discharge_coefficient=1, cross_sectional_area=1.0e-4)


Inlet = IsothermalIncompressibleVolume("Suction Inlet", PumpNetwork, pressure=SuctionFluid.pressure,
                                       temperature=SuctionFluid.temperature, mass_flow_in=InletLine.mass_flow, 
                                       mass_flow_out=1.8,
                                       volume=0.2*IN3_TO_M3)

pump_mdot = Inlet.mass_flow_out

Q_out = pump_mdot / DischargeFluid.density


def pump_head(Q_out, N):
    speed_ratio = N / N_design
    phi = Q_out / (Q_design * speed_ratio)

    H_shutoff = 1.35 * H_design
    k = 0.259

    H = H_shutoff * speed_ratio**2 * (1.0 - k * phi**2)

    return H.clip(lower=0.0)


def pump_torque(Q_out, N, rho_out):
    H = pump_head(Q_out, N)
    omega = np.pi / 30.0 * N

    mdot = rho_out * Q_out
    hydraulic_power = mdot * 9.80665 * H
    shaft_power = hydraulic_power / eta_design

    return shaft_power / omega


FuelEpump = PolytropicPump("Fuel E-Pump",
                           PumpNetwork,
                           rotor_speed=shaft_speed,
                           head_rise=pump_head(Q_out, shaft_speed),
                           mass_flow=pump_mdot,
                           torque=pump_torque(Q_out, shaft_speed, DischargeFluid.density),
                           upstream_density=SuctionFluid.density, 
                           downstream_density=DischargeFluid.density,
                           upstream_total_pressure=Inlet.pressure,
                           discharge_total_pressure=DischargeFluid.pressure,
                           upstream_total_enthalpy=SuctionFluid.enthalpy)

Outlet = SimpleVolume("Discharge",
                      PumpNetwork,
                      pressure=DischargeFluid.pressure,
                      enthalpy=DischargeFluid.enthalpy,
                      density=DischargeFluid.density,
                      mass_flow_in=FuelEpump.mass_flow,
                      enthalpy_in=FuelEpump.discharge_total_enthalpy,
                      volume=0.2*IN3_TO_M3)

OutletLine = GenericDarcyWeisbach("Outlet Line", PumpNetwork, upstream_pressure=Outlet.pressure, 
                                   downstream_pressure=101325, length=1.5, cross_sectional_area=0.2*IN2_TO_M2,
                                   hydraulic_diameter=0.5*IN_TO_M, roughness=0.1e-3, density=Outlet.density, 
                                   dynamic_viscosity=DischargeFluid.dynamic_viscosity, mass_flow=Outlet.mass_flow_out)

solution = SteadyState(PumpNetwork).solve(return_type='dataframe', verbose=True, static=False)
print(solution.to_string(index=False))






# ----- Polytropic Pump Summary -----
mdot = FuelEpump.mass_flow.value

rho1 = FuelEpump.upstream_density.value
rho2 = FuelEpump.downstream_density.value

Q1 = FuelEpump.inlet_volumetric_flow.value
Q2 = FuelEpump.outlet_volumetric_flow.value

N = FuelEpump.rotor_speed.value
omega = np.pi / 30.0 * N

H_m = FuelEpump.head_rise.value
H_specific = FuelEpump.gravitational_acceleration.value * H_m

T = FuelEpump.torque.value
eta = FuelEpump.efficiency.value
shaft_power = FuelEpump.shaft_power.value

po_in = FuelEpump.upstream_total_pressure.value
po_out_node = FuelEpump.discharge_total_pressure.value
po_out_pred = FuelEpump._predicted_discharge_total_pressure

dp = po_out_node - po_in

ho_in = FuelEpump.upstream_total_enthalpy.value
ho_out = FuelEpump.discharge_total_enthalpy.value
dho = ho_out - ho_in

pressure_ratio = po_out_node / po_in
density_ratio = rho2 / rho1

beta = 1.0 / (
    1.0 - np.log(density_ratio) / np.log(pressure_ratio)
)

print("\n" + "=" * 70)
print("                    POLYTROPIC PUMP SUMMARY")
print("=" * 70)

print(f"{'Fluid':40s}: {fuel}")
print(f"{'Rotor Speed':40s}: {N:14.3f} RPM")
print(f"{'Angular Speed':40s}: {omega:14.3f} rad/s")

print("\n--- Flow ---")
print(f"{'Mass Flow':40s}: {mdot:14.6f} kg/s")
print(f"{'Inlet Density':40s}: {rho1:14.3f} kg/m^3")
print(f"{'Outlet Density':40s}: {rho2:14.3f} kg/m^3")
print(f"{'Density Ratio rho2/rho1':40s}: {density_ratio:14.6f}")
print(f"{'Inlet Volumetric Flow':40s}: {Q1:14.6f} m^3/s")
print(f"{'Outlet Volumetric Flow':40s}: {Q2:14.6f} m^3/s")
print(f"{'Inlet Volumetric Flow':40s}: {Q1 / IN3_TO_M3 * 60.0:14.3f} in^3/min")
print(f"{'Outlet Volumetric Flow':40s}: {Q2 / IN3_TO_M3 * 60.0:14.3f} in^3/min")

print("\n--- Pressure / Head ---")
print(f"{'Inlet Total Pressure':40s}: {po_in:14.3f} Pa")
print(f"{'Inlet Total Pressure':40s}: {po_in / PSIA_TO_PA:14.3f} psia")
print(f"{'Outlet Node Total Pressure':40s}: {po_out_node:14.3f} Pa")
print(f"{'Outlet Node Total Pressure':40s}: {po_out_node / PSIA_TO_PA:14.3f} psia")
print(f"{'Pump Predicted Outlet Pressure':40s}: {po_out_pred:14.3f} Pa")
print(f"{'Pump Predicted Outlet Pressure':40s}: {po_out_pred / PSIA_TO_PA:14.3f} psia")
print(f"{'Pressure Rise':40s}: {dp:14.3f} Pa")
print(f"{'Pressure Rise':40s}: {dp / PSIA_TO_PA:14.3f} psid")
print(f"{'Pressure Ratio po2/po1':40s}: {pressure_ratio:14.6f}")
print(f"{'Head Rise':40s}: {H_m:14.3f} m")
print(f"{'Specific Head gH':40s}: {H_specific:14.3f} J/kg")
print(f"{'Polytropic Beta':40s}: {beta:14.6f}")

print("\n--- Enthalpy ---")
print(f"{'Inlet Total Enthalpy':40s}: {ho_in:14.3f} J/kg")
print(f"{'Pump Discharge Total Enthalpy':40s}: {ho_out:14.3f} J/kg")
print(f"{'Enthalpy Rise':40s}: {dho:14.3f} J/kg")
print(f"{'Ideal Head Work gH':40s}: {H_specific:14.3f} J/kg")
print(f"{'Loss / Heating Contribution':40s}: {dho - H_specific:14.3f} J/kg")

print("\n--- Power / Torque ---")
print(f"{'Torque':40s}: {T:14.6f} N-m")
print(f"{'Shaft Power':40s}: {shaft_power:14.3f} W")
print(f"{'Shaft Power':40s}: {shaft_power / 1000.0:14.3f} kW")
print(f"{'Efficiency':40s}: {eta:14.6f}")
print(f"{'Efficiency':40s}: {100.0 * eta:14.3f} %")
