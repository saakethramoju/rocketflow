from System import *
from Solvers import *
'''
# Pure fluid: pressure-temperature
water = Fluid("Water", pressure=101325, temperature=300)
print(water)

print('----------------------------')
# Pure fluid: density-internal energy transient-style lookup
water = Fluid("Water", density=997.0, internal_energy=1.1e5)
print(water.pressure)
print(water.temperature)
print(water.enthalpy)

print('----------------------------')
# Update existing object using pair setters
water.pressure_temperature = (2e5, 350)
print(water.density)

water.density_internal_energy = (950, 2.0e5)
print(water.pressure, water.temperature)

print('----------------------------')
# Single-property setters
water.density = 900.0          # holds internal_energy fixed
print(water.pressure, water.temperature)

water.internal_energy = 2.5e5  # holds density fixed
print(water.pressure, water.temperature)

print('----------------------------')
# Saturated state: pressure-quality
sat_water = Fluid("Water", pressure=101325, quality=1.0)
print(sat_water.temperature)
print(sat_water.enthalpy)
print(sat_water.density)

print('----------------------------')
# Mixture: air-like nitrogen/oxygen/argon by mole fraction
air_mix = Fluid(
    {
        "Nitrogen": 0.78,
        "Oxygen": 0.21,
        "Argon": 0.01,
    },
    basis="mole",
    pressure=101325,
    temperature=300,
)

print(air_mix)
print(air_mix.mole_fractions)
print(air_mix.mass_fractions)


print('----------------------------')

print(Fluid.available_flash_pairs())

print('----------------------------')
# Mixture: pressure-enthalpy flash
air_ph = Fluid(
    {
        "Nitrogen": 0.78,
        "Oxygen": 0.21,
        "Argon": 0.01,
    },
    basis="mole",
    pressure=101325,
    enthalpy=air_mix.enthalpy,
)

print(air_ph.temperature)
print(air_ph.density)

print('----------------------------')
# Mixture: update using pair setter
air_mix.pressure_temperature = (2e5, 500)
print(air_mix.enthalpy)
print(air_mix.specific_heat_cp)
print(air_mix.specific_heat_ratio)
'''
'''

from Utilities import IdealGas


def section(title: str):
    print()
    print("=" * 80)
    print(title)
    print("=" * 80)


def show(gas: IdealGas):
    print(gas)
    print(f"pressure_temperature          = {gas.pressure_temperature}")
    print(f"pressure_enthalpy             = {gas.pressure_enthalpy}")
    print(f"pressure_internal_energy      = {gas.pressure_internal_energy}")


section("1. Nitrogen from temperature only")
n2_temperature = IdealGas("Nitrogen", temperature=300)
show(n2_temperature)


section("2. Nitrogen from enthalpy only")
n2_enthalpy = IdealGas("Nitrogen", enthalpy=n2_temperature.enthalpy)
show(n2_enthalpy)


section("3. Nitrogen from internal energy only")
n2_internal_energy = IdealGas(
    "Nitrogen",
    internal_energy=n2_temperature.internal_energy,
)
show(n2_internal_energy)


section("4. Nitrogen from pressure + temperature")
n2_pressure_temperature = IdealGas(
    "Nitrogen",
    pressure=101325,
    temperature=300,
)
show(n2_pressure_temperature)


section("5. Nitrogen from pressure + enthalpy")
n2_pressure_enthalpy = IdealGas(
    "Nitrogen",
    pressure=101325,
    enthalpy=n2_pressure_temperature.enthalpy,
)
show(n2_pressure_enthalpy)


section("6. Nitrogen from pressure + internal energy")
n2_pressure_internal_energy = IdealGas(
    "Nitrogen",
    pressure=101325,
    internal_energy=n2_pressure_temperature.internal_energy,
)
show(n2_pressure_internal_energy)


section("7. Nitrogen from density + temperature")
rho = n2_pressure_temperature.density

n2_density_temperature = IdealGas(
    "Nitrogen",
    density=rho,
    temperature=300,
)
show(n2_density_temperature)


section("8. Nitrogen from density + enthalpy")
n2_density_enthalpy = IdealGas(
    "Nitrogen",
    density=rho,
    enthalpy=n2_pressure_temperature.enthalpy,
)
show(n2_density_enthalpy)


section("9. Nitrogen from density + internal energy")
n2_density_internal_energy = IdealGas(
    "Nitrogen",
    density=rho,
    internal_energy=n2_pressure_temperature.internal_energy,
)
show(n2_density_internal_energy)


section("10. Nitrogen from pressure + density")
n2_pressure_density = IdealGas(
    "Nitrogen",
    pressure=101325,
    density=rho,
)
show(n2_pressure_density)


section("11. Test setters")
gas = IdealGas("Nitrogen", temperature=300)

gas.pressure = 101325
print("After pressure set:")
show(gas)

gas.pressure_temperature = (202650, 400)
print("After pressure_temperature set:")
show(gas)

gas.pressure_enthalpy = (101325, n2_pressure_temperature.enthalpy)
print("After pressure_enthalpy set:")
show(gas)

gas.pressure_internal_energy = (
    101325,
    n2_pressure_temperature.internal_energy,
)
print("After pressure_internal_energy set:")
show(gas)

gas.density_temperature = (rho, 300)
print("After density_temperature set:")
show(gas)

gas.density_enthalpy = (rho, n2_pressure_temperature.enthalpy)
print("After density_enthalpy set:")
show(gas)

gas.density_internal_energy = (
    rho,
    n2_pressure_temperature.internal_energy,
)
print("After density_internal_energy set:")
show(gas)

gas.pressure_density = (101325, rho)
print("After pressure_density set:")
show(gas)


section("12. Mixture from temperature only")
air_temperature = IdealGas(
    {
        "Nitrogen": 0.78,
        "Oxygen": 0.21,
        "Argon": 0.01,
    },
    basis="mole",
    temperature=300,
)
show(air_temperature)


section("13. Mixture from pressure + temperature")
air_pressure_temperature = IdealGas(
    {
        "Nitrogen": 0.78,
        "Oxygen": 0.21,
        "Argon": 0.01,
    },
    basis="mole",
    pressure=101325,
    temperature=300,
)
show(air_pressure_temperature)


section("14. Mixture from enthalpy only")
air_enthalpy = IdealGas(
    {
        "Nitrogen": 0.78,
        "Oxygen": 0.21,
        "Argon": 0.01,
    },
    basis="mole",
    enthalpy=air_pressure_temperature.enthalpy,
)
show(air_enthalpy)


section("15. Mixture from internal energy only")
air_internal_energy = IdealGas(
    {
        "Nitrogen": 0.78,
        "Oxygen": 0.21,
        "Argon": 0.01,
    },
    basis="mole",
    internal_energy=air_pressure_temperature.internal_energy,
)
show(air_internal_energy)


section("16. Mixture from density + internal energy")
air_density_internal_energy = IdealGas(
    {
        "Nitrogen": 0.78,
        "Oxygen": 0.21,
        "Argon": 0.01,
    },
    basis="mole",
    density=air_pressure_temperature.density,
    internal_energy=air_pressure_temperature.internal_energy,
)
show(air_density_internal_energy)


section("17. Confirm matching states")
print("Nitrogen T from temperature only      :", n2_temperature.temperature)
print("Nitrogen T from enthalpy only         :", n2_enthalpy.temperature)
print("Nitrogen T from internal_energy only  :", n2_internal_energy.temperature)
print("Nitrogen T from density+u             :", n2_density_internal_energy.temperature)
print("Nitrogen P from density+u             :", n2_density_internal_energy.pressure)
print("Nitrogen rho from P+T                 :", n2_pressure_temperature.density)
print("Nitrogen rho from density+u           :", n2_density_internal_energy.density)

print()
print("Air T from P+T                        :", air_pressure_temperature.temperature)
print("Air T from enthalpy only              :", air_enthalpy.temperature)
print("Air T from internal_energy only       :", air_internal_energy.temperature)
print("Air T from density+u                  :", air_density_internal_energy.temperature)
print("Air P from P+T                        :", air_pressure_temperature.pressure)
print("Air P from density+u                  :", air_density_internal_energy.pressure)

'''
from Utilities import FluidRegistry

FluidRegistry.add_alias("pressurant", "Helium")

SimpleNetwork = Network("Simple Network")

supply_gas = IdealGasLookup("Supply Gas", SimpleNetwork, "Pressurant", enthalpy=1553814.3007257418)

f = FluidLookup("Fluid Lookup", SimpleNetwork, "pressurant", pressure=101325, enthalpy=1553814.3007257418)

print(supply_gas.temperature)
print(f.temperature)

print(FluidRegistry.aliases)