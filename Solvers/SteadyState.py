from __future__ import annotations
from typing import TYPE_CHECKING

import numpy as np
from scipy.optimize import least_squares, Bounds, root

if TYPE_CHECKING:
    from System import Network


class SteadyState:

    def __init__(self, network: Network):
        self.network = network
    '''
    def residual(self, x: np.ndarray) -> np.ndarray:
        self.network.assign_iteration_values(list(x))
        # self.network.pre_evaluation()
        self.network.evaluate_states()
        return np.array(self.network.scaled_residuals, dtype=float)
    '''

    def residual(self, x: np.ndarray) -> np.ndarray:
        self.network.assign_iteration_values(list(x))
        n_residuals = len(self.network.residuals)
        try:
            self.network.evaluate_states()
            return np.array(self.network.residuals, dtype=float)
        except Exception:
            penalty = np.ones(n_residuals, dtype=float) * 1e6
            # add x-dependence so finite-difference Jacobian is nonzero
            penalty += 1e-3 * np.linalg.norm(x)
            return penalty
        

    def static_evaluate(
        self,
        filename: str | None = None,
        return_type: str = "dict",
        verbose: bool = True,
    ):
        """
        Evaluate a static network without nonlinear solving.

        Runs pre_evaluation() and evaluate_states(), then saves/returns
        the network state.
        """
        self.network.pre_evaluation()
        self.network.evaluate_states()

        if verbose:
            print("\n" + "=" * 50)
            print("        STATIC NETWORK EVALUATION")
            print("=" * 50)
            print("  No nonlinear solve was performed.")
            print("=" * 50 + "\n")

        return self.network.save(filename=filename, return_type=return_type)


    def solve(
        self,
        filename: str | None = None,
        return_type: str = "dict",
        verbose: bool = True,
        static: bool = False,
    ):
        """
        Solve the network steady-state nonlinear system.

        Behavior
        --------
        - If `static=True`, the network is only evaluated once
        without nonlinear iteration.

        - If the network has no iteration variables and no residual
        equations, the network is only evaluated once without nonlinear
        iteration.

        - `root()` is used for square, unconstrained systems:
            number of residuals == number of iteration variables

        - `least_squares()` is automatically used when:
            - any variable bounds exist
            - the system is overdetermined
                (more residuals than iteration variables)

        - If `root()` fails to converge, the solver automatically
        falls back to `least_squares()`.
        """

        if static:
            return self.static_evaluate(
                filename=filename,
                return_type=return_type,
                verbose=verbose,
            )

        self.network.pre_evaluation()

        x0 = np.array(self.network.iteration_values, dtype=float)

        self.network.evaluate_states()

        r0 = np.array(self.network.residuals, dtype=float)

        no_iteration_variables = len(x0) == 0
        no_residuals = len(r0) == 0

        if no_iteration_variables and no_residuals:
            return self.static_evaluate(
                filename=filename,
                return_type=return_type,
                verbose=verbose,
            )

        if len(r0) < len(x0):
            raise ValueError(
                f"SteadyState solve requires at least as many residuals as "
                f"iteration variables. Got {len(x0)} iteration variables "
                f"and {len(r0)} residuals."
            )

        bounds = Bounds(
            self.network.lower_bounds,
            self.network.upper_bounds,
            self.network.keep_feasible,
        )

        switched_to_least_squares = False

        if self.network.has_bounds or len(r0) > len(x0):
            solver_name = "scipy.optimize.least_squares"
            sol = least_squares(self.residual, x0, bounds=bounds)

        else:
            solver_name = "scipy.optimize.root"
            sol = root(self.residual, x0)

            if not sol.success:
                switched_to_least_squares = True
                solver_name = "scipy.optimize.least_squares (fallback)"
                sol = least_squares(self.residual, x0, bounds=bounds)

        if verbose and switched_to_least_squares:
            print(
                "\n[Warning] root() failed, so the solver switched to "
                "least_squares().\n"
            )

        if verbose:
            self._verbose_print(sol, solver_name)

        final_residual = np.array(sol.fun, dtype=float)

        if not sol.success or np.max(np.abs(final_residual)) > 1e3:
            raise RuntimeError(
                "Steady-state solve failed or converged to penalty residuals.\n"
                f"success = {sol.success}\n"
                f"message = {sol.message}\n"
                f"max |residual| = {np.max(np.abs(final_residual)):.3e}"
            )

        self.network.assign_iteration_values(list(sol.x))
        self.network.evaluate_states()

        solution = self.network.save(filename=filename, return_type=return_type)

        return solution



    def _verbose_print(self, sol, solver_name: str) -> None:
        print("\n" + "=" * 50)
        print("        STEADY-STATE SOLVER SUMMARY")
        print("=" * 50)

        print("\n[Solver]")
        print(f"  Method         : {solver_name}")

        print("\n[Convergence]")
        print(f"  Success        : {sol.success}")
        print(f"  Status         : {sol.status}")
        print(f"  Message        : {sol.message}")

        print("\n[Performance]")
        print(f"  Function evals : {sol.nfev}")
        if hasattr(sol, "njev"):
            print(f"  Jacobian evals : {sol.njev}")

        print("\n[Optimality]")
        if hasattr(sol, "cost"):
            print(f"  Cost (½‖r‖²)   : {sol.cost:.6e}")
        if hasattr(sol, "optimality"):
            print(f"  Optimality     : {sol.optimality:.3e}")

        print("\n[Solution Variables]")
        iter_vars = self.network.iteration_variables
        if isinstance(iter_vars, str):
            iter_vars = [iter_vars]

        if len(iter_vars) == len(sol.x):
            for name, val in zip(iter_vars, sol.x):
                print(f"  {str(name):<40} = {val:.6e}")
        else:
            for i, val in enumerate(sol.x):
                print(f"  x[{i:<2}] {'':<34} = {val:.6e}")

        print("\n[Residuals]")
        for i, r in enumerate(sol.fun):
            print(f"  r[{i:<2}] {'':<34} = {r:.3e}")

        max_resid = np.max(np.abs(sol.fun))
        rms_resid = np.sqrt(np.mean(sol.fun**2))

        print("\n[Residual Summary]")
        print(f"  Max |residual| : {max_resid:.3e}")
        print(f"  RMS residual   : {rms_resid:.3e}")

        print("=" * 50 + "\n")