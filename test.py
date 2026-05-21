
'''
# ---- Pump Data -----
design_rotor_speed = 20000.0       # rpm
design_volumetric_flow = 0.0024    # m^3/s
design_head_rise = 46.5            # m
design_torque = 0.75               # N-m

g = 9.80665
density = 750.0                    # kg/m^3
impeller_diameter = 0.040          # m


omega_design = np.pi / 30.0 * design_rotor_speed

design_flow_coefficient = (
    design_volumetric_flow
    / (omega_design * impeller_diameter**3)
)

design_head_coefficient = (
    design_head_rise
    / (omega_design**2 * impeller_diameter**2 / g)
)

design_torque_coefficient = (
    design_torque
    / (density * omega_design**2 * impeller_diameter**5)
)


normalized_flow_coefficient_map = [
    0.00,
    0.40,
    0.70,
    1.00,
    1.20,
    1.40,
]

# reverse-sigmoid-ish head curve
normalized_head_coefficient_map = [
    1.35,
    1.28,
    1.15,
    1.00,
    0.78,
    0.45,
]

# torque forced to zero at zero flow and increases with flow
normalized_torque_coefficient_map = [
    0.00,
    0.35,
    0.65,
    1.00,
    1.22,
    1.45,
]


# ---- Components ----

EPumpMap = TurboMap("E-Pump Map", MapNetwork, 
                    rotor_speed=25000,
                    volumetric_flow=2/density,
                    density=density,
                    flow_geometric_parameter=impeller_diameter,
                    head_geometric_parameter=impeller_diameter,
                    torque_geometric_parameter=impeller_diameter,
                    design_flow_coefficient=design_flow_coefficient,
                    design_head_coefficient=design_head_coefficient,
                    design_torque_coefficient=design_torque_coefficient,
                    normalized_flow_coefficient_map=normalized_flow_coefficient_map,
                    normalized_head_coefficient_map=normalized_head_coefficient_map,
                    normalized_torque_coefficient_map=normalized_torque_coefficient_map,
                    )

solution = SteadyState(MapNetwork).solve(return_type='dataframe', verbose=True, static=False)
print(solution.to_string(index=False))
'''


import numpy as np
import matplotlib.pyplot as plt

gamma = 1.4

M = np.linspace(0.01, 5, 1000)

F = M * (1 + (gamma - 1)/2 * M**2)**(
    -(gamma + 1)/(2 * (gamma - 1))
)

plt.figure(figsize=(8,5))
plt.plot(M, F)

plt.axvline(1, linestyle="--")
plt.xlabel("Mach Number")
plt.ylabel("F(M)")
plt.title("Compressible Mass Flow Function")
plt.grid(True)

plt.show()