from System import *
from Solvers import *
from constants import PA_PER_PSI

HETS = Network("HETS")

fuel_tank_pressure = State(400 * PA_PER_PSI)
ox_tank_pressure = State(400 * PA_PER_PSI)
fuel_inj_pressure = State(300 * PA_PER_PSI)
ox_inj_pressure = State(300 * PA_PER_PSI)
fuel_density = State(800)
ox_density = State(1104)
fuel_runline_mdot = State(0)
ox_runline_mdot = State(0)
fuel_injector_mdot = State(0)
ox_injector_mdot = State(0)
chamber_pressure = State(200 * PA_PER_PSI)
chamber_density = State(1.6)
nozzle_mdot = State(0)
atmospheric_pressure = State(0)

FuelTank = PressureNode("Fuel Tank", 
                        network=HETS,
                        pressure=fuel_tank_pressure)

OxTank = PressureNode("Fuel Tank", 
                        network=HETS,
                        pressure=ox_tank_pressure)


FuelRunline = DischargeCoefficient("Fuel Runline",
                                   network=HETS,
                                   upstream_pressure=fuel_tank_pressure,
                                   downstream_pressure=fuel_inj_pressure,
                                   density=fuel_density,
                                   discharge_coefficient=1,
                                   cross_sectional_area=0.5e-4,
                                   mass_flow=fuel_runline_mdot)

OxRunline = DischargeCoefficient("Ox Runline",
                                   network=HETS,
                                   upstream_pressure=ox_tank_pressure,
                                   downstream_pressure=ox_inj_pressure,
                                   density=ox_density,
                                   discharge_coefficient=1,
                                   cross_sectional_area=0.5e-4,
                                   mass_flow=ox_runline_mdot)

FuelInjectorManifold = SimpleIncompressibleVolume("Fuel Injector Manifold",
                                                  network=HETS,
                                                  pressure=fuel_inj_pressure,
                                                  density=fuel_density,
                                                  volume=0.1287,
                                                  mass_flow_in=fuel_runline_mdot,
                                                  mass_flow_out=fuel_injector_mdot)

OxInjectorManifold = SimpleIncompressibleVolume("Ox Injector Manifold",
                                                  network=HETS,
                                                  pressure=ox_inj_pressure,
                                                  density=ox_density,
                                                  volume=0.1287,
                                                  mass_flow_in=ox_runline_mdot,
                                                  mass_flow_out=ox_injector_mdot)

FuelInjector = DischargeCoefficient("Fuel Injector",
                                   network=HETS,
                                   upstream_pressure=fuel_inj_pressure,
                                   downstream_pressure=chamber_pressure,
                                   density=fuel_density,
                                   discharge_coefficient=1,
                                   cross_sectional_area=0.5e-4,
                                   mass_flow=fuel_injector_mdot)

OxInjector = DischargeCoefficient("Ox Injector",
                                   network=HETS,
                                   upstream_pressure=ox_inj_pressure,
                                   downstream_pressure=chamber_pressure,
                                   density=ox_density,
                                   discharge_coefficient=1,
                                   cross_sectional_area=0.5e-4,
                                   mass_flow=ox_injector_mdot)

chamber = SimpleIncompressibleVolume("Combustion Chamber",
                                        network=HETS,
                                        pressure=chamber_pressure,
                                        density=chamber_density,
                                        volume=0.1287,
                                        mass_flow_in=fuel_injector_mdot + ox_injector_mdot,
                                        mass_flow_out=nozzle_mdot,
                                    )


Nozzle = DischargeCoefficient("Nozzle",
                                   network=HETS,
                                   upstream_pressure=chamber_pressure,
                                   downstream_pressure=atmospheric_pressure,
                                   density=chamber_density,
                                   discharge_coefficient=1,
                                   cross_sectional_area=1e-4,
                                   mass_flow=nozzle_mdot)

Ambient = PressureNode("Atmosphere", 
                       network=HETS,
                       pressure=atmospheric_pressure)

print(HETS)
solver = SteadyState(HETS)
df = solver.solve(return_type='dataframe')
print(df)
