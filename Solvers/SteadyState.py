from __future__ import annotations
from typing import TYPE_CHECKING

import numpy as np
from scipy.optimize import least_squares, Bounds

if TYPE_CHECKING:
    from System import Network


class SteadyState:

    def __init__(self,
                 network: Network,
                 method: str | None = None):

        self.network = network
        self.method = method

    def residual(self, x: np.ndarray) -> np.ndarray:
        self.network.assign_iteration_values(list(x))
        self.network.pre_evaluation()
        self.network.evaluate_states()
        return np.array(self.network.residuals, dtype=float)

    def solve(self, filename: str | None = None, return_type: str = "dict"):
        x0 = np.array(self.network.iteration_values, dtype=float)

        kwargs = {}

        # least_squares methods: 'trf', 'dogbox', 'lm'
        if self.method is not None:
            kwargs["method"] = self.method

        bounds = Bounds(self.network.lower_bounds, 
                        self.network.upper_bounds,
                        self.network.keep_feasible)

        sol = least_squares(self.residual, x0, bounds=bounds, **kwargs)

        # assign final solution back into network
        self.network.assign_iteration_values(list(sol.x))
        self.network.pre_evaluation()
        self.network.evaluate_states()

        solution = self.network.save(filename=filename, return_type=return_type)

        return solution