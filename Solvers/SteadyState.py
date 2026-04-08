from __future__ import annotations
from typing import TYPE_CHECKING

import numpy as np
from scipy.optimize import root

if TYPE_CHECKING:
    from System import Network

class SteadyState:

    def __init__(self, 
                 network: Network, 
                 tolerance: float | None = None,
                 method: str | None = None,
                 options: dict | None = None):
        
        self.network = network
        self.tolerance = tolerance
        self.method = method
        self.options = options

    def residual(self, x: np.ndarray) -> np.ndarray:
        self.network.assign_iteration_values(list(x))
        self.network.pre_evaluation()
        self.network.evaluate_states()
        return np.array(self.network.residuals, dtype=float)
        

    def solve(self, filename: str | None = None, return_type: str = "dict"):
        x0 = np.array(self.network.iteration_values, dtype=float)

        kwargs = {}

        if self.method is not None:
            kwargs["method"] = self.method
        if self.tolerance is not None:
            kwargs["tol"] = self.tolerance
        if self.options is not None:
            kwargs["options"] = self.options

        sol = root(self.residual, x0, **kwargs)

        # assign final solution back into network
        self.network.assign_iteration_values(list(sol.x))
        self.network.pre_evaluation()
        self.network.evaluate_states()

        solution = self.network.save(filename=filename, return_type=return_type)

        return solution