from __future__ import annotations
from typing import TYPE_CHECKING

import numpy as np
from scipy.optimize import least_squares, Bounds, root
from rich.console import Console
from rich.table import Table
from rich import box

from Solvers.SparsityDetector import SparsityDetector

if TYPE_CHECKING:
    from System import Network


class SteadyState:

    def __init__(self, network: Network):
        self.network = network
        self.console = Console()

    def residual(self, x: np.ndarray) -> np.ndarray:
        self.network.assign_iteration_values(list(x))

        try:
            self.network.evaluate_states()
            return np.array(self.network.residuals, dtype=float)

        except Exception as e:

            raise RuntimeError(
                "Solver encountered an error while evaluating the network "
                "inside network.evaluate_states().\n\n"
                f"Original error:\n{type(e).__name__}: {e}"
            ) from e

    def static_evaluate(
        self,
        filename: str | None = None,
        return_type: str = "dict",
        verbose: bool = False,
        print_solution: bool = False,
    ):
        """
        Evaluate a static network without nonlinear solving.

        Runs pre_evaluation() and evaluate_states(), then saves/returns
        the network state.
        """
        self.network.pre_evaluation()
        self.network.evaluate_states()

        if verbose:
            self._verbose_static_print()

        solution = self.network.save(
            filename=filename,
            return_type=return_type,
        )

        if print_solution:
            self.print_solution()

        return solution

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
        auto_sparsity: bool = True,
        sparsity_relative_step: float = 1e-5,
        sparsity_absolute_step: float = 1e-7,
        sparsity_residual_tolerance: float = 1e-12,
        plot_sparsity: bool = False,
        sparsity_plot_filename: str = "jacobian_sparsity.png",
    ):

        if static:
            return self.static_evaluate(
                filename=filename,
                return_type=return_type,
                verbose=verbose,
                print_solution=print_solution,
            )

        method = solver_method.lower()

        valid_solver_methods = ("trf", "dogbox", "lm")

        if method not in valid_solver_methods:
            raise ValueError(
                "solver_method must be one of "
                f"{valid_solver_methods}. "
                f"Got '{solver_method}'."
            )

        jac = jacobian_method.lower()

        valid_jac_methods = ("2-point", "3-point")

        if jac not in valid_jac_methods:
            raise ValueError(
                "jacobian_method must be one of "
                f"{valid_jac_methods}. "
                f"Got '{jacobian_method}'."
            )

        if rtol <= 0.0:
            raise ValueError(
                f"Residual tolerance (rtol) must be positive. "
                f"Got {rtol}"
            )

        self.network.pre_evaluation()

        x0 = np.array(
            self.network.iteration_values,
            dtype=float,
        )

        self.network.evaluate_states()

        r0 = np.array(
            self.network.residuals,
            dtype=float,
        )

        no_iteration_variables = len(x0) == 0
        no_residuals = len(r0) == 0

        if no_iteration_variables and no_residuals:
            return self.static_evaluate(
                filename=filename,
                return_type=return_type,
                verbose=verbose,
                print_solution=print_solution,
            )

        if len(r0) < len(x0):
            raise ValueError(
                "SteadyState solve requires at least as many "
                "residuals as iteration variables. "
                f"Got {len(x0)} iteration variables "
                f"and {len(r0)} residuals."
            )

        overconstrained = len(r0) > len(x0)

        bounds = Bounds(
            self.network.lower_bounds,
            self.network.upper_bounds,
            self.network.keep_feasible,
        )

        jac_sparsity = None

        if auto_sparsity:
            if method == "lm":
                raise ValueError("auto_sparsity cannot be used with solver_method='lm'.")

            detector = SparsityDetector(self.network)

            jac_sparsity = detector.detect(
                relative_step=sparsity_relative_step,
                absolute_step=sparsity_absolute_step,
                residual_tolerance=sparsity_residual_tolerance,
            )

            if plot_sparsity:
                detector.plot(
                    jac_sparsity,
                    filename=sparsity_plot_filename,
                )

        sol = least_squares(
            fun=self.residual,
            x0=x0,
            method=method,
            bounds=bounds,
            x_scale="jac",
            jac=jac,
            jac_sparsity=jac_sparsity,
            #tr_solver="lsmr" if jac_sparsity is not None else None,
            #tr_options={"regularize": True} if jac_sparsity is not None else None,
            ftol=ftol,
            xtol=xtol,
            gtol=gtol,
        )

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
            )

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

        self.network.assign_iteration_values(
            list(sol.x)
        )

        self.network.evaluate_states()

        solution = self.network.save(
            filename=filename,
            return_type=return_type,
        )

        if print_solution:
            self.print_solution()

        return solution

    def print_solution(self) -> None:

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

    def _verbose_static_print(self) -> None:

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
        table.add_row("Components", str(len(self.network.components)))
        table.add_row(
            "Iteration variables",
            str(len(self.network.iteration_variables)),
        )
        table.add_row("Residuals", str(len(self.network.residuals)))

        self.console.print()
        self.console.print(table)
        self.console.print()


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
    ) -> None:

        max_residual = np.max(np.abs(sol.fun))
        rms_residual = np.sqrt(np.mean(sol.fun**2))

        dx = np.array(sol.x, dtype=float) - np.array(x0, dtype=float)

        normalized_correction = np.max(
            np.abs(dx) / np.maximum(np.abs(sol.x), 1.0)
        )

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

            labels = []

            for component in self.network.component_list:

                for attr_name, attr_value in component.__dict__.items():

                    if attr_value is target:

                        labels.append(
                            f"{component.name}.{attr_name}"
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


'''
class SteadyStateOld:

    def __init__(self, network: Network):
        self.network = network
    """
    def residual(self, x: np.ndarray) -> np.ndarray:
        self.network.assign_iteration_values(list(x))
        # self.network.pre_evaluation()
        self.network.evaluate_states()
        return np.array(self.network.scaled_residuals, dtype=float)
    """

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
'''