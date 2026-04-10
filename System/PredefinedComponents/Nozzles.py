from __future__ import annotations

from typing import TYPE_CHECKING

from System import Component, State
from Utilities import create_SI_CEA_object

if TYPE_CHECKING:
    from System import Network


class RocketCEAChokedNozzle(Component):

    def __init__(self, 
                 name: str,
                 network: Network,
                 fuel: str,
                 oxidizer: str,
                 chamber_pressure: State,
                 mixture_ratio: State,
                 throat_area: float,
                 expansion_ratio: float,
                 ambient_pressure: State,
                 characterstic_velocity_efficiency: float,
                 thrust_coefficient_efficiency: float,
                 thrust: State,
                 mass_flow: State):
        self.initialize_component(name, network)

        self.fuel = fuel
        self.ox = oxidizer
        self.Pc = chamber_pressure
        self.MR = mixture_ratio
        self.At = State(throat_area)
        self.eps = State(expansion_ratio)
        self.Pamb = ambient_pressure
        self.F = thrust
        self.mdot = mass_flow
        self.eta_cstar = State(characterstic_velocity_efficiency)
        self.eta_Cf = State(thrust_coefficient_efficiency)
        self._cea_obj = create_SI_CEA_object(self.fuel, self.ox)


    def evaluate_states(self) -> None:
         Pc = self.Pc.value
         MR = self.MR.value
         At = self.At.value

         cstar_ideal = self._cea_obj.get_Cstar(Pc, MR)
         _, Cf_ideal, _ = self._cea_obj.get_PambCf(self.Pamb.value, Pc, MR, self.eps.value)

         self.mdot.value = Pc * At / (self.eta_cstar.value * cstar_ideal)
         self.F.value = self.eta_Cf.value * Cf_ideal * Pc * At

    @property
    def iteration_variables(self) -> list[State]:
        return []

    @property
    def residuals(self) -> list[float]:
        return []