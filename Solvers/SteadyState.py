# SteadyState.py
from __future__ import annotations
from typing import TYPE_CHECKING

import numpy as np
from scipy.optimize import least_squares, Bounds, root
from rich.console import Console
from rich.table import Table
from rich import box
import time

if TYPE_CHECKING:
    from System import Network


class SteadyState:
    """
    Steady-state nonlinear solver.

    Overall solving process
    -----------------------
    1. Run network.pre_evaluation()
       Components can initialize internal helper states.

    2. Collect the current iteration variables into a vector x0.
       These are solver-owned States such as pressures, enthalpies,
       friction factors, pump speeds, and Balance variables.

    3. For each solver residual call:
          a. Assign the solver vector x into the network States.
          b. Evaluate component states repeatedly until derived states settle.
          c. Restore iteration variables during settling so components cannot
             overwrite the solver's current iterate.
          d. Collect component residuals and balance residuals.

    4. scipy.optimize.least_squares finite-differences the residual function,
       builds a Jacobian, and computes a trust-region correction.

    5. The process repeats until scipy reports convergence and the maximum
       residual is below rtol.

    """

    def __init__(self, network: Network):
        self.network = network
        self.console = Console()

    # ------------------------------------------------------------------
    # Residual function passed to scipy
    # ------------------------------------------------------------------

    def residual(self, x: np.ndarray) -> np.ndarray:
        """
        Map solver vector x to residual vector r.

        scipy calls this many times while solving.
        """
        # Assign current solver iterate to network State objects.
        self.network.assign_iteration_values(list(x))

        try:
            # Propagate explicit component states to a consistent state.
            self.evaluate_network_states()

            # Collect component and balance residuals.
            return np.array(self.network.residuals, dtype=float)

        except Exception as e:
            raise RuntimeError(
                "Solver encountered an error while evaluating the network "
                "inside evaluate_network_states().\n\n"
                f"Original error:\n{type(e).__name__}: {e}"
            ) from e

    # ------------------------------------------------------------------
    # Order-independent state evaluation for steady-state solving
    # ------------------------------------------------------------------

    def evaluate_network_states(
        self,
        max_passes: int = 20,
        tolerance: float = 1e-10,
    ) -> None:
        """
        Repeatedly evaluate components until derived State values settle.

        This reduces dependence on user component order. For example:

            composition -> density -> mass_flow

        may require several passes before the final mass flow is consistent
        with the final composition.

        Iteration variables are protected because they belong to the nonlinear
        solver, not to explicit component propagation.
        """
        # Save current solver-owned variables.
        iteration_snapshot = self._snapshot_iteration_variables()

        for _ in range(max_passes):
            # Track non-iteration states before this pass.
            old_values = self._collect_state_values()

            for c in self.network.component_list:
                # Restore solver variables before each component.
                self._restore_iteration_variables(iteration_snapshot)

                # Let the component update its derived states.
                c.evaluate_states()

                # Restore solver variables again in case the component changed them.
                self._restore_iteration_variables(iteration_snapshot)

            # Track non-iteration states after this pass.
            new_values = self._collect_state_values()

            # Stop when another pass would not materially change the network.
            if self._max_state_change(old_values, new_values) < tolerance:
                return

    def _snapshot_iteration_variables(self) -> dict[int, float]:
        """
        Save current values of all solver-owned iteration variables.

        These values represent the current solver iterate x.
        """
        snapshot = {}

        for var in self.network.collect_all_iteration_variables():
            if var.is_assigned:
                snapshot[id(var)] = float(var.value)

        return snapshot

    def _restore_iteration_variables(
        self,
        snapshot: dict[int, float],
    ) -> None:
        """
        Restore solver-owned iteration variables.

        This prevents evaluate_states() from permanently overwriting x.
        """
        for var in self.network.collect_all_iteration_variables():
            key = id(var)

            if key in snapshot:
                var.value = snapshot[key]

    def _collect_state_values(self) -> dict[int, float]:
        """
        Collect assigned non-iteration State values.

        These values are used only to decide whether explicit state propagation
        has settled.
        """
        values = {}

        # Iteration variables are ignored because they are fixed by the solver.
        iteration_ids = {
            id(var)
            for var in self.network.collect_all_iteration_variables()
        }

        def collect(attr_value):
            # Composition-like object.
            if hasattr(attr_value, "fraction"):
                for _, state in attr_value:
                    if id(state) in iteration_ids:
                        continue

                    if state.is_assigned:
                        try:
                            values[id(state)] = float(state.value)
                        except Exception:
                            pass
                return

            # State-like object.
            if hasattr(attr_value, "is_assigned"):
                if id(attr_value) in iteration_ids:
                    return

                if attr_value.is_assigned:
                    try:
                        values[id(attr_value)] = float(attr_value.value)
                    except Exception:
                        pass

        # Scan component attributes.
        for comp in self.network.component_list:
            for attr_value in comp.__dict__.values():
                collect(attr_value)

        # Scan balance attributes too.
        for bal in self.network.balance_list:
            for attr_value in bal.__dict__.values():
                collect(attr_value)

        return values

    def _max_state_change(
        self,
        old: dict[int, float],
        new: dict[int, float],
    ) -> float:
        """
        Return the largest normalized change between two state snapshots.

        A small value means the explicit state propagation has settled.
        """
        max_change = 0.0

        for key, new_value in new.items():
            old_value = old.get(key)

            # Newly assigned state.
            if old_value is None:
                max_change = max(max_change, abs(new_value))
                continue

            # Relative change, protected from division by tiny values.
            scale = max(abs(new_value), 1.0)

            max_change = max(
                max_change,
                abs(new_value - old_value) / scale,
            )

        return max_change

    # ------------------------------------------------------------------
    # Static evaluation path
    # ------------------------------------------------------------------

    def static_evaluate(
        self,
        filename: str | None = None,
        return_type: str = "dict",
        verbose: bool = False,
        print_solution: bool = False,
    ):
        """
        Evaluate a network without nonlinear solving.

        This still uses steady-state state-settling so component order is less
        important even for static evaluations.
        """
        start_time = time.perf_counter()

        self.network.pre_evaluation()
        self.evaluate_network_states()

        elapsed_time = time.perf_counter() - start_time

        if verbose:
            self._verbose_static_print(
                elapsed_time=elapsed_time,
            )

        solution = self.network.save(
            filename=filename,
            return_type=return_type,
        )

        if print_solution:
            self.print_solution()

        return solution

    # ------------------------------------------------------------------
    # Main nonlinear solve
    # ------------------------------------------------------------------

    def solve(
        self,
        filename: str | None = None,
        return_type: str = "dict",
        verbose: bool = False,
        static: bool = False,
        print_solution: bool = False,
        solver_method: str = "trf",
        jacobian_method: str = "3-point",
        ftol: float = 1e-8,
        xtol: float = 1e-8,
        gtol: float = 1e-8,
        rtol: float = 1e-2,
    ):
        """
        Solve the steady-state nonlinear system.

        The unknowns are all component iteration variables plus Balance
        variables. The residuals are all component residuals plus Balance
        residuals.
        """

        # Static mode skips nonlinear solving.
        if static:
            return self.static_evaluate(
                filename=filename,
                return_type=return_type,
                verbose=verbose,
                print_solution=print_solution,
            )

        # Validate selected least_squares method.
        method = solver_method.lower()

        valid_solver_methods = ("trf", "dogbox", "lm")

        if method not in valid_solver_methods:
            raise ValueError(
                "solver_method must be one of "
                f"{valid_solver_methods}. "
                f"Got '{solver_method}'."
            )

        # Validate selected finite-difference Jacobian method.
        jac = jacobian_method.lower()

        valid_jac_methods = ("2-point", "3-point")

        if jac not in valid_jac_methods:
            raise ValueError(
                "jacobian_method must be one of "
                f"{valid_jac_methods}. "
                f"Got '{jacobian_method}'."
            )

        # Validate residual acceptance tolerance.
        if rtol <= 0.0:
            raise ValueError(
                f"Residual tolerance (rtol) must be positive. "
                f"Got {rtol}"
            )

        # One-time component setup.
        self.network.pre_evaluation()

        # Initial solver vector.
        x0 = np.array(
            self.network.iteration_values,
            dtype=float,
        )

        # Initial state propagation and residual vector.
        self.evaluate_network_states()

        r0 = np.array(
            self.network.residuals,
            dtype=float,
        )

        # If there is nothing to solve, just export current network state.
        no_iteration_variables = len(x0) == 0
        no_residuals = len(r0) == 0

        if no_iteration_variables and no_residuals:
            return self.static_evaluate(
                filename=filename,
                return_type=return_type,
                verbose=verbose,
                print_solution=print_solution,
            )

        # least_squares requires at least as many residuals as unknowns.
        if len(r0) < len(x0):
            raise ValueError(
                "SteadyState solve requires at least as many "
                "residuals as iteration variables. "
                f"Got {len(x0)} iteration variables "
                f"and {len(r0)} residuals."
            )

        overconstrained = len(r0) > len(x0)

        # Build scipy Bounds object from State metadata.
        bounds = Bounds(
            self.network.lower_bounds,
            self.network.upper_bounds,
            self.network.keep_feasible,
        )

        # start timing
        start_time = time.perf_counter()

        # Main nonlinear solve.
        sol = least_squares(
            fun=self.residual,
            x0=x0,
            method=method,
            bounds=bounds,
            x_scale="jac",
            jac=jac,
            ftol=ftol,
            xtol=xtol,
            gtol=gtol,
        )

        # end time
        elapsed_time = time.perf_counter() - start_time

        # Optional rich solver summary.
        if verbose:
            self._verbose_print(
                sol=sol,
                x0=x0,
                method=method,
                jac=jac,
                ftol=ftol,
                xtol=xtol,
                gtol=gtol,
                rtol=rtol,
                overconstrained=overconstrained,
                elapsed_time=elapsed_time,
            )

        # Check final residual quality.
        final_residual = np.array(
            sol.fun,
            dtype=float,
        )

        max_residual = np.max(np.abs(final_residual))

        if (
            not sol.success
            or max_residual > rtol
        ):
            raise RuntimeError(
                "Steady-state solve failed or converged "
                "to unacceptable residuals.\n"
                f"success = {sol.success}\n"
                f"message = {sol.message}\n"
                f"max |residual| = {max_residual:.3e}\n"
                f"residual tolerance = {rtol:.3e}"
            )

        # Load final solver result into network.
        self.network.assign_iteration_values(
            list(sol.x)
        )

        # Re-evaluate final derived states using final solver variables.
        self.evaluate_network_states()

        # Save/export final solution.
        solution = self.network.save(
            filename=filename,
            return_type=return_type,
        )

        if print_solution:
            self.print_solution()

        return solution

    # ------------------------------------------------------------------
    # Solution printing
    # ------------------------------------------------------------------

    def print_solution(self) -> None:
        """Print exported network state as a rich table."""
        records = self.network.save(return_type="dict")

        table = Table(
            title=f"{self.network.name} Solution",
            box=box.SIMPLE_HEAVY,
            show_header=True,
            header_style="bold",
        )

        table.add_column("Component", style="#D84135", no_wrap=True)
        table.add_column("Type", style="#3B629E", no_wrap=True)
        table.add_column("Attribute", style="#fdf0d5", no_wrap=True)
        table.add_column("Value", justify="right")

        for record in records:
            value = record["value"]

            if isinstance(value, float):
                value_text = f"{value:.6g}"
            else:
                value_text = str(value)

            if value_text == "<uninitialized>":
                value_text = "[dim]<uninitialized>[/dim]"

            elif value_text == "<unavailable>":
                value_text = "[red]<unavailable>[/red]"

            table.add_row(
                str(record["component_name"]),
                str(record["component_type"]),
                str(record["attribute"]),
                value_text,
            )

        self.console.print()
        self.console.print(table)
        self.console.print()

    # ------------------------------------------------------------------
    # Verbose static evaluation printing
    # ------------------------------------------------------------------

    def _verbose_static_print(
        self,
        elapsed_time: float,
    ) -> None:
        """Print a short summary for static evaluation mode."""
        table = Table(
            title="Static Network Evaluation",
            box=box.SIMPLE_HEAVY,
            show_header=True,
            header_style="bold",
        )

        table.add_column("Quantity", style="bold")
        table.add_column("Value", justify="right")

        table.add_row("Mode", "Static evaluation")
        table.add_row("Nonlinear solve", "Not performed")
        table.add_row("Evaluation time",f"{elapsed_time:.3f} s")
        table.add_row("Components", str(len(self.network.components)))
        table.add_row(
            "Iteration variables",
            str(len(self.network.iteration_variables)),
        )
        table.add_row("Residuals", str(len(self.network.residuals)))

        self.console.print()
        self.console.print(table)
        self.console.print()

    # ------------------------------------------------------------------
    # Verbose nonlinear solve printing
    # ------------------------------------------------------------------

    def _verbose_print(
        self,
        sol,
        x0: np.ndarray,
        method: str,
        jac: str,
        ftol: float,
        xtol: float,
        gtol: float,
        rtol: float,
        overconstrained: bool = False,
        elapsed_time: float = 0.0,
    ) -> None:
        """Print solver summary, final variables, and final residuals."""

        max_residual = np.max(np.abs(sol.fun))
        rms_residual = np.sqrt(np.mean(sol.fun**2))

        dx = np.array(sol.x, dtype=float) - np.array(x0, dtype=float)

        normalized_correction = np.max(
            np.abs(dx) / np.maximum(np.abs(sol.x), 1.0)
        )

        # --------------------------------------------------------------
        # Solver summary table.
        # --------------------------------------------------------------
        summary = Table(
            title="Steady-State Solver Summary",
            box=box.SIMPLE_HEAVY,
            show_header=True,
            header_style="bold",
        )

        summary.add_column("Quantity", style="bold")
        summary.add_column("Value", justify="right")

        summary.add_row("Success", str(sol.success))
        summary.add_row("Status", str(sol.status))
        summary.add_row("Message", str(sol.message))

        if overconstrained:
            summary.add_row(
                "Warning",
                "System is overconstrained",
                style="yellow",
            )

        summary.add_row("Solver method", method)
        summary.add_row("Jacobian method", jac)
        summary.add_row("Function evaluations", str(sol.nfev))
        summary.add_row("Solve time",f"{elapsed_time:.3f} s")

        if hasattr(sol, "njev") and sol.njev is not None:
            summary.add_row("Jacobian evaluations", str(sol.njev))

        if hasattr(sol, "cost"):
            summary.add_row("Cost", f"{sol.cost:.6e}")

        if hasattr(sol, "optimality"):
            summary.add_row("Optimality", f"{sol.optimality:.3e}")

        summary.add_row("Max |residual|", f"{max_residual:.3e}")
        summary.add_row("RMS residual", f"{rms_residual:.3e}")
        summary.add_row("Max normalized correction", f"{normalized_correction:.3e}")
        summary.add_row("Residual tolerance", f"{rtol:.3e}")
        summary.add_row("ftol", f"{ftol:.3e}")
        summary.add_row("xtol", f"{xtol:.3e}")
        summary.add_row("gtol", f"{gtol:.3e}")

        # --------------------------------------------------------------
        # Iteration variable table.
        # --------------------------------------------------------------
        variables = Table(
            title="Solution Variables",
            box=box.SIMPLE,
            show_header=True,
            header_style="bold",
        )

        variables.add_column("Index", justify="right", style="dim")
        variables.add_column("Variable", style="#fdf0d5")
        variables.add_column("Value", justify="right", style="#D84135")

        def find_variable_labels(target):
            """Find all component/balance attributes that reference a State."""
            labels = []

            for component in self.network.component_list:
                for attr_name, attr_value in component.__dict__.items():
                    if attr_value is target:
                        labels.append(
                            f"{component.name}.{attr_name}"
                        )

                    if hasattr(attr_value, "fraction"):
                        for species, state in attr_value:
                            if state is target:
                                labels.append(
                                    f"{component.name}.{attr_name}.{species}"
                                )

            for balance in self.network.balance_list:
                for attr_name, attr_value in balance.__dict__.items():
                    if attr_value is target:
                        labels.append(
                            f"{balance.name}.{attr_name}"
                        )

            if labels:
                return labels

            return [str(target)]

        variable_labels = [
            find_variable_labels(var)
            for var in self.network.collect_all_iteration_variables()
        ]

        for i, val in enumerate(sol.x):
            label = (
                "\n".join(variable_labels[i])
                if i < len(variable_labels)
                else "<unlabeled>"
            )

            variables.add_row(
                f"x[{i}]",
                label,
                f"{val:.6e}",
            )

        # --------------------------------------------------------------
        # Residual table.
        # --------------------------------------------------------------
        residuals = Table(
            title="Residuals",
            box=box.SIMPLE,
            show_header=True,
            header_style="bold",
        )

        residuals.add_column("Index", justify="right", style="dim")
        residuals.add_column("Residual", style="#fdf0d5")
        residuals.add_column("Value", justify="right", style="#3B629E")

        residual_labels = []

        for component in self.network.component_list:
            component_residuals = component.residuals

            if isinstance(component_residuals, (list, tuple)):
                for i in range(len(component_residuals)):
                    residual_labels.append(f"{component.name}.residual[{i}]")
            else:
                residual_labels.append(f"{component.name}.residual")

        for balance in self.network.balance_list:
            balance_residuals = balance.residuals

            if isinstance(balance_residuals, (list, tuple)):
                for i in range(len(balance_residuals)):
                    residual_labels.append(f"{balance.name}.residual[{i}]")
            else:
                residual_labels.append(f"{balance.name}.residual")

        for i, r in enumerate(sol.fun):
            label = (
                residual_labels[i]
                if i < len(residual_labels)
                else "<unlabeled>"
            )

            residuals.add_row(
                f"r[{i}]",
                label,
                f"{r:.6e}",
            )

        self.console.print()
        self.console.print(summary)
        self.console.print(variables)
        self.console.print(residuals)
        self.console.print()