from __future__ import annotations

import numpy as np
from typing import TYPE_CHECKING

from System import Component, State

if TYPE_CHECKING:
    from System import Network



class ConstantDensityPump(Component):

    def __init__(self,
                 name: str, 
                 network: Network,
                 rotor_speed: State,
                 head_rise: State,
                 volumetric_flow: State,
                 density: State,
                 torque: State,
                 upstream_total_pressure: State,
                 discharge_total_pressure: State,
                 upstream_total_enthalpy: State | None = None,
                 discharge_total_enthalpy: State | None = None,
                 gravitational_acceleration: float = 9.80665,
                 efficiency: State | None = None,
                 shaft_power: State | None = None,
                 mass_flow: State | None = None):
        self.setup()

        self._predicted_discharge_total_pressure = None
        self._predicted_discharge_total_enthalpy = None

    def evaluate_states(self):
        H = self.head_rise.value 
        Q = self.volumetric_flow.value
        g = self.gravitational_acceleration.value
        T = self.torque.value
        rho = self.density.value
        po_in = self.upstream_total_pressure.value
        N = self.rotor_speed.value

        omega = (np.pi / 30.0) * N

        if abs(Q) < 1e-12:
            raise ValueError(f"{self.name}: volumetric_flow is too close to zero.")

        if abs(rho * Q) < 1e-12:
            raise ValueError(f"{self.name}: mass flow is too close to zero.")

        if abs(T * omega) < 1e-12:
            raise ValueError(f"{self.name}: shaft power is too close to zero.")

        shaft_power = T * omega
        hydraulic_power = rho * g * H * Q

        po_out = po_in + rho * g * H
        eta = hydraulic_power / shaft_power
        mdot = rho * Q

        self._predicted_discharge_total_pressure = po_out

        self.efficiency.value = eta
        self.shaft_power.value = shaft_power
        self.mass_flow.value = mdot

        if self.upstream_total_enthalpy.is_assigned:
            ho_in = self.upstream_total_enthalpy.value
            dho = shaft_power / mdot
            ho_out = ho_in + dho

            self._predicted_discharge_total_enthalpy = ho_out

            if self.discharge_total_enthalpy.is_assigned:
                self.discharge_total_enthalpy.value = ho_out
        else:
            self._predicted_discharge_total_enthalpy = None

    @property
    def iteration_variables(self) -> list[State]:
        return [self.volumetric_flow]
    
    @property
    def residuals(self) -> list[float]:
        return [
            self._predicted_discharge_total_pressure
            - self.discharge_total_pressure.value
        ]


class SimpleEulerCentrifugalPump(Component):

    def __init__(self, 
                 name: str, 
                 network: Network,
                 rotor_speed: State, # rpm
                 volumetric_flow: State,
                 density: State,
                 impeller_inlet_tip_radius: float,
                 impeller_outlet_tip_radius: float,
                 inlet_annular_flow_area: float,
                 outlet_annular_flow_area: float,
                 inlet_blade_angle: float,
                 outlet_blade_angle: float,
                 angle_units: str = "degrees",
                 gravitational_acceleration: float = 9.80665,
                 slip_factor: float = 1.0,
                 hydraulic_efficiency: float | State = 1.0,
                 mechanical_efficiency: float | State = 1.0,
                 volumetric_efficiency: float | State = 1.0,
                 head_rise: State | None = None,
                 torque: State | None = None,
                 shaft_power: State | None = None,
                 stagnation_pressure_rise: State | None = None,
                 mass_flow: State | None = None,):
        self.setup()

    def pre_evaluation(self):
        self.evaluate_states()

    def evaluate_states(self):
        N = self.rotor_speed.value              # rpm
        r1 = self.impeller_inlet_tip_radius.value
        r2 = self.impeller_outlet_tip_radius.value
        A1 = self.inlet_annular_flow_area.value
        A2 = self.outlet_annular_flow_area.value
        Q = self.volumetric_flow.value          # delivered external flow rate
        rho = self.density.value

        eta_h = self.hydraulic_efficiency.value
        eta_v = self.volumetric_efficiency.value
        eta_m = self.mechanical_efficiency.value

        beta1 = self.inlet_blade_angle.value
        beta2 = self.outlet_blade_angle.value
        sigma = self.slip_factor.value
        g = self.gravitational_acceleration.value

        if self.angle_units.lower() in {"degree", "degrees", "deg"}:
            beta1 = np.deg2rad(beta1)
            beta2 = np.deg2rad(beta2)
        elif self.angle_units.lower() in {"radian", "radians", "rad"}:
            pass
        else:
            raise ValueError(
                "angle_units must be 'degrees' or 'radians'."
            )

        if abs(N) < 1e-12:
            raise ValueError("shaft_speed must be nonzero for pump torque calculation.")

        if eta_h <= 0.0 or eta_v <= 0.0 or eta_m <= 0.0:
            raise ValueError("Pump efficiencies must be positive.")

        omega = np.pi * N / 30.0

        U1 = omega * r1
        U2 = omega * r2

        # The impeller internally processes more flow than the delivered flow
        # if volumetric efficiency is below 1.
        Q_impeller = Q / eta_v

        Cm1 = Q_impeller / A1
        Cm2 = Q_impeller / A2

        V_theta1 = U1 - Cm1 / np.tan(beta1)
        V_theta2 = sigma * U2 - Cm2 / np.tan(beta2)

        specific_work_euler = U2 * V_theta2 - U1 * V_theta1

        H_euler = specific_work_euler / g
        H_actual = eta_h * H_euler

        hydraulic_power = rho * g * H_actual * Q

        # Since hydraulic efficiency already reduced the delivered head,
        # only mechanical efficiency remains in shaft power.
        shaft_power = hydraulic_power / eta_m

        torque = shaft_power / omega

        self.head_rise.value = H_actual
        self.shaft_power.value = shaft_power
        self.torque.value = torque
    



class PolytropicPump(Component): pass