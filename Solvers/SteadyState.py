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
    Steady-state nonlinear network solver.

    Overview
    --------
    This solver finds a set of iteration variables x that drives all network
    residual equations to zero.

    Unknowns may include:

        - pressures
        - enthalpies
        - temperatures
        - mass flows
        - friction factors
        - pump speeds
        - user Balance variables

    Residuals may include:

        - mass conservation
        - energy conservation
        - momentum equations
        - friction-factor equations
        - user-defined component residuals
        - Balance residuals

    Mathematical Form
    -----------------
    The solver seeks:

        F(x) = 0

    where:

        x = iteration variables
        F = network residual vector

    The nonlinear system is solved using scipy.optimize.least_squares()
    with finite-difference Jacobians and a trust-region algorithm.

    State Propagation Architecture
    ------------------------------
    Most network quantities are not iteration variables.

    Examples:

        density
        viscosity
        Reynolds number
        species composition
        fluid properties
        heat-transfer coefficients

    These quantities are derived from the current iteration variables.

    A challenge in large networks is that derived quantities may depend on
    other derived quantities:

        composition
            -> density
                -> mass flow
                    -> composition

    and users may define components in any order.

    To make component ordering largely irrelevant, each residual evaluation
    contains a state-settling phase.

    For a given solver iterate x:

        1. Assign x to the network.
        2. Evaluate all components.
        3. Measure changes in non-iteration States.
        4. Repeat until the network state stops changing.

    This produces a self-consistent set of derived quantities before
    residuals are evaluated.

    Iteration Variable Protection
    -----------------------------
    Iteration variables belong exclusively to the nonlinear solver.

    During state settling, component evaluations are allowed to update
    derived States but are not allowed to permanently modify solver-owned
    variables.

    To enforce this:

        - iteration variables are snapshotted
        - component evaluations are performed
        - iteration variables are restored

    before and after every component evaluation pass.

    This prevents explicit state calculations from corrupting the solver's
    current iterate.

    Convergence Process
    -------------------
    For each nonlinear iteration:

        x_k
        ↓
        assign iteration variables
        ↓
        settle derived states
        ↓
        evaluate residuals
        ↓
        finite-difference Jacobian
        ↓
        trust-region correction
        ↓
        x_(k+1)

    The process repeats until scipy reports convergence and the final
    residual magnitude satisfies the requested residual tolerance.

    Benefits
    --------
    This architecture provides:

        - user-defined component ordering
        - automatic propagation of derived quantities
        - support for algebraic balances
        - support for overdetermined systems
        - separation between solver variables and derived states
        - compatibility with future transient solvers
    """

    def __init__(self, network: Network):
        self.network = network
        self.console = Console()

        # Debug snapshots used when scipy fails before returning a solution.
        # These let verbose=True still print the last attempted variables,
        # variable adjustments, and any residuals that were available.
        self._last_debug_x = None
        self._last_debug_residual = None
        self._last_debug_error = None

        # Default state-settling controls used inside residual(), since scipy
        # calls residual(x) without passing solve() keyword arguments.
        self._state_max_passes = 20
        self._state_tolerance = 1e-10

    # ------------------------------------------------------------------
    # Residual function passed to scipy
    # ------------------------------------------------------------------

    def residual(self, x: np.ndarray) -> np.ndarray:
        """
        Map solver vector x to residual vector r.

        scipy calls this many times while solving.
        """
        # Save the current attempted iterate for failure diagnostics.
        self._last_debug_x = np.array(x, dtype=float)
        self._last_debug_residual = None
        self._last_debug_error = None

        # Assign current solver iterate to network State objects.
        self.network.assign_iteration_values(list(x))

        try:
            # Propagate explicit component states to a consistent state.
            self.evaluate_network_states(
                max_passes=self._state_max_passes,
                tolerance=self._state_tolerance,
            )

            # Collect component and balance residuals.
            r = np.array(self.network.residuals, dtype=float)
            self._last_debug_residual = r

            return r

        except Exception as e:
            # Preserve the original error and whatever residuals can still
            # be collected, then re-raise so failures are not masked.
            self._last_debug_error = e

            try:
                self._last_debug_residual = np.array(
                    self.network.residuals,
                    dtype=float,
                )
            except Exception:
                self._last_debug_residual = None

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
    

    def _get_model(self, model):
        """
        Return a Model object from either a Model instance or model name.
        """

        if model is None:
            return None

        if hasattr(model, "build") and hasattr(model, "available_options"):
            return model

        for candidate in self.network.model_list:
            if candidate.name == model:
                return candidate

        raise ValueError(
            f"Unknown model {model!r}. "
            f"Available models are {[m.name for m in self.network.model_list]}."
        )


    def _safe_sheet_name(self, name: str) -> str:
        """
        Return an Excel-safe worksheet name.
        """

        invalid_chars = ["\\", "/", "*", "[", "]", ":", "?"]

        for char in invalid_chars:
            name = name.replace(char, "_")

        return name[:31]


    def _save_model_option_results(
        self,
        results: dict,
        filename: str,
    ) -> None:
        """
        Save static model-option evaluations.

        For .xlsx/.xls, each model option is written to a separate sheet.
        For .json, all model options are written as one keyed dictionary.
        For .csv, each model option is written to a separate CSV file.
        """

        ext = filename.split(".")[-1].lower()
        base = ".".join(filename.split(".")[:-1])

        if ext == "json":
            import json

            with open(filename, "w") as f:
                json.dump(results, f, indent=4)

            return

        if ext in {"xlsx", "xls"}:
            import pandas as pd

            with pd.ExcelWriter(filename) as writer:
                for option_name, records in results.items():
                    df = pd.DataFrame(records)
                    df.to_excel(
                        writer,
                        sheet_name=self._safe_sheet_name(option_name),
                        index=False,
                    )

            return

        if ext == "csv":
            import pandas as pd

            for option_name, records in results.items():
                df = pd.DataFrame(records)
                df.to_csv(
                    f"{base}_{option_name}.csv",
                    index=False,
                )

            return

        raise ValueError(
            "Unsupported file extension. Use .csv, .json, or .xlsx"
        )
        
    def _build_unbuilt_models(self) -> None:
        """
        Build the first option for every unbuilt registered Model.
        """

        for model in self.network.model_list:
            if model.active_component is None:
                model.build()


    def static_evaluate(
        self,
        model: str | None = None,
        filename: str | None = None,
        return_type: str = "dict",
        verbose: bool = False,
        print_solution: bool = False,
        state_max_passes: int = 20,
        state_tolerance: float = 1e-10,
    ):
        """
        Evaluate a network without nonlinear solving.

        If model is None, the first option from every unbuilt Model is built
        and the network is evaluated once.

        If model is provided, every option in that Model is evaluated. The return
        value is keyed by model option name. Printed output shows only the final
        active option.
        """

        start_time = time.perf_counter()

        selected_model = self._get_model(model)

        # --------------------------------------------------------------
        # Evaluate every option in one requested Model.
        # --------------------------------------------------------------
        if selected_model is not None:
            results = {}

            for option_name in selected_model.order:
                # Replace the active model option before each evaluation.
                selected_model.replace(option_name)

                self.network.pre_evaluation()
                self.evaluate_network_states(
                    max_passes=state_max_passes,
                    tolerance=state_tolerance,
                )

                # Always collect raw records first so file export is consistent.
                records = self.network.save(
                    filename=None,
                    return_type="dict",
                )

                if return_type == "dict":
                    results[option_name] = records

                elif return_type == "dataframe":
                    import pandas as pd

                    results[option_name] = pd.DataFrame(records)

                else:
                    raise ValueError("return_type must be 'dict' or 'dataframe'")

            elapsed_time = time.perf_counter() - start_time

            if filename is not None:
                # File export always uses raw record dictionaries.
                raw_results = {
                    option_name: (
                        result.to_dict(orient="records")
                        if hasattr(result, "to_dict")
                        else result
                    )
                    for option_name, result in results.items()
                }

                self._save_model_option_results(
                    raw_results,
                    filename,
                )

            # Print only the final active option.
            if verbose:
                self._verbose_static_print(
                    elapsed_time=elapsed_time,
                )

            if print_solution:
                self.print_solution()

            return results

        # --------------------------------------------------------------
        # Normal static evaluation.
        # --------------------------------------------------------------

        # Static evaluation should include model components by default.
        # If no model is specified, build the first option for every unbuilt model.
        self._build_unbuilt_models()

        self.network.pre_evaluation()
        self.evaluate_network_states(
            max_passes=state_max_passes,
            tolerance=state_tolerance,
        )

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
        model: str | None = None,
        filename: str | None = None,          # Optional file path/name for saving the exported solution.
        return_type: str = "dict",            # Solution export format passed to network.save().
        verbose: bool = False,                # If True, print solver summary, variables, and residuals.
        static: bool = False,                 # If True, skip nonlinear solving and only evaluate the network.
        print_solution: bool = False,         # If True, print the final exported network state table.
        solver_method: str = "trf",           # scipy least_squares method: "trf", "dogbox", or "lm".
        jacobian_method: str = "3-point",     # Finite-difference Jacobian method: "2-point" or "3-point".
        ftol: float = 1e-8,                   # scipy cost-function termination tolerance.
        xtol: float = 1e-8,                   # scipy solution-step termination tolerance.
        gtol: float = 1e-8,                   # scipy gradient/optimality termination tolerance.
        rtol: float = 1e-2,                   # Maximum acceptable final absolute residual.
        state_max_passes: int = 20,           # Max derived-state settling passes per residual evaluation.
        state_tolerance: float = 1e-10,       # Normalized derived-state settling tolerance.
    ):

        """
        Solve the network's steady-state nonlinear system.

        This method finds values for all iteration variables such that the
        network residuals are driven to zero (or as close to zero as
        possible within the requested tolerances).

        The solve process is:

            1. Run component pre-evaluation.
            2. Collect all iteration variables.
            3. Propagate derived states until the network settles.
            4. Evaluate residuals.
            5. Use scipy.optimize.least_squares() to solve the system.
            6. Re-evaluate the network using the converged solution.
            7. Export and optionally print the final results.

        Parameters
        ----------
        filename : str | None, optional
            Optional file path passed to network.save().

        return_type : str, optional
            Format returned by network.save().
            Common options are "dict" and "dataframe".

        verbose : bool, optional
            If True, print detailed solver diagnostics including:

                - convergence information
                - solve time
                - iteration variables
                - variable adjustments
                - residual values

            Diagnostic information is also printed when the solve fails.

        static : bool, optional
            If True, skip nonlinear solving and perform only a network
            evaluation using the current State values.

        print_solution : bool, optional
            If True, print the final exported network solution table.

        solver_method : str, optional
            scipy.optimize.least_squares() algorithm.

            Supported options:

                "trf"     Trust Region Reflective (recommended)
                "dogbox"  Dogleg trust region
                "lm"      Levenberg-Marquardt

        jacobian_method : str, optional
            Finite-difference Jacobian approximation.

            Supported options:

                "2-point"  Faster
                "3-point"  More accurate

        ftol : float, optional
            Cost-function convergence tolerance passed to scipy.

        xtol : float, optional
            Iteration-variable convergence tolerance passed to scipy.

        gtol : float, optional
            Gradient/optimality convergence tolerance passed to scipy.

        rtol : float, optional
            Maximum acceptable absolute residual after convergence.

            The solve is considered unsuccessful if:

                max(abs(residual)) > rtol

            even if scipy reports convergence.

        state_max_passes : int, optional
            Maximum number of state-settling passes performed during each
            residual evaluation.

            Higher values may improve convergence for strongly coupled
            networks but increase solve time.

        state_tolerance : float, optional
            Normalized convergence tolerance used by the state-settling
            algorithm.

            Smaller values enforce tighter consistency between derived
            states but increase solve time.

        Returns
        -------
        Any
            Result returned by network.save(return_type=...).

        Raises
        ------
        ValueError
            If solver settings are invalid or the system is underdetermined.

        RuntimeError
            If the nonlinear solve fails, encounters an invalid network
            state, or converges to residuals larger than rtol.

        Notes
        -----
        This solver supports:

            - user-defined components
            - user-defined balances
            - overdetermined systems
            - variable bounds
            - derived-state propagation
            - order-independent component evaluation

        The reported "Variable Adjustment" is:

            x_final - x_initial

        and the reported "Normalized Variable Adjustment" is:

            abs(x_final - x_initial) / max(abs(x_final), 1)

        which indicates how far each iteration variable moved from its
        initial guess during the solve.
        """

        # Store state-settling controls for residual(), since scipy only calls
        # residual(x) and does not pass solve() keyword arguments through.
        self._state_max_passes = state_max_passes
        self._state_tolerance = state_tolerance

        # Static mode skips nonlinear solving.
        if static:
            return self.static_evaluate(
                filename=filename,
                return_type=return_type,
                verbose=verbose,
                print_solution=print_solution,
                state_max_passes=state_max_passes,
                state_tolerance=state_tolerance,
                model=model,
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

        # Validate state-settling controls.
        if state_max_passes <= 0:
            raise ValueError(
                f"state_max_passes must be positive. "
                f"Got {state_max_passes}"
            )

        if state_tolerance <= 0.0:
            raise ValueError(
                f"state_tolerance must be positive. "
                f"Got {state_tolerance}"
            )

        # One-time component setup.
        self.network.pre_evaluation()

        # Initial solver vector.
        x0 = np.array(
            self.network.iteration_values,
            dtype=float,
        )

        # Initial state propagation and residual vector.
        self.evaluate_network_states(
            max_passes=state_max_passes,
            tolerance=state_tolerance,
        )

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
                state_max_passes=state_max_passes,
                state_tolerance=state_tolerance,
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

        try:
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

        except Exception:
            elapsed_time = time.perf_counter() - start_time

            if verbose:
                self._verbose_failure_print(
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

            raise

        # end time
        elapsed_time = time.perf_counter() - start_time

        # Check final residual quality.
        final_residual = np.array(
            sol.fun,
            dtype=float,
        )

        max_residual = np.max(np.abs(final_residual))

        # Optional rich solver summary.
        # This is intentionally printed before the residual-quality failure
        # check so verbose=True still reports failed convergences.
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
        self.evaluate_network_states(
            max_passes=state_max_passes,
            tolerance=state_tolerance,
        )

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
    # Shared verbose-label helpers
    # ------------------------------------------------------------------

    def _find_variable_labels(self, target):
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

    def _collect_residual_labels(self) -> list[str]:
        """Collect printed labels for component and balance residuals."""
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

        return residual_labels

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

        # Variable adjustment from initial guess to final solution.
        normalized_variable_adjustment = (
            np.abs(dx) / np.maximum(np.abs(sol.x), 1.0)
        )

        max_normalized_variable_adjustment = np.max(
            normalized_variable_adjustment
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
        summary.add_row("Solve time",f"{elapsed_time:.3f} s")
        summary.add_row("Function evaluations", str(sol.nfev))

        if hasattr(sol, "njev") and sol.njev is not None:
            summary.add_row("Jacobian evaluations", str(sol.njev))

        if hasattr(sol, "cost"):
            summary.add_row("Cost", f"{sol.cost:.6e}")

        if hasattr(sol, "optimality"):
            summary.add_row("Optimality", f"{sol.optimality:.3e}")

        summary.add_row("Max |residual|", f"{max_residual:.3e}")
        summary.add_row("RMS residual", f"{rms_residual:.3e}")
        summary.add_row(
            "Max normalized variable adjustment",
            f"{max_normalized_variable_adjustment:.3e}",
        )
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
        variables.add_column("Variable Adjustment", justify="right", style="#3B629E")
        variables.add_column("Normalized Variable Adjustment", justify="right", style="#3B629E")

        variable_labels = [
            self._find_variable_labels(var)
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
                f"{dx[i]:+.6e}",
                f"{normalized_variable_adjustment[i]:.3e}",
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

        residual_labels = self._collect_residual_labels()

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

    def _verbose_failure_print(
        self,
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
        """
        Print the last attempted solver state when scipy fails before
        returning a solution object.
        """
        x = self._last_debug_x
        r = self._last_debug_residual
        error = self._last_debug_error

        summary = Table(
            title="Steady-State Solver Failure",
            box=box.SIMPLE_HEAVY,
            show_header=True,
            header_style="bold",
        )

        summary.add_column("Quantity", style="bold")
        summary.add_column("Value", justify="right")

        summary.add_row("Success", "False")
        summary.add_row("Solver method", method)
        summary.add_row("Jacobian method", jac)
        summary.add_row("Solve time", f"{elapsed_time:.3f} s")

        if overconstrained:
            summary.add_row(
                "Warning",
                "System is overconstrained",
                style="yellow",
            )

        if error is not None:
            summary.add_row("Error type", type(error).__name__)
            summary.add_row("Error", str(error))

        summary.add_row("Residual tolerance", f"{rtol:.3e}")
        summary.add_row("ftol", f"{ftol:.3e}")
        summary.add_row("xtol", f"{xtol:.3e}")
        summary.add_row("gtol", f"{gtol:.3e}")

        self.console.print()
        self.console.print(summary)

        if x is not None:
            dx = np.array(x, dtype=float) - np.array(x0, dtype=float)

            # Variable adjustment from initial guess to last attempted iterate.
            normalized_variable_adjustment = (
                np.abs(dx) / np.maximum(np.abs(x), 1.0)
            )

            variables = Table(
                title="Last Attempted Solution Variables",
                box=box.SIMPLE,
                show_header=True,
                header_style="bold",
            )

            variables.add_column("Index", justify="right", style="dim")
            variables.add_column("Variable", style="#fdf0d5")
            variables.add_column("Value", justify="right", style="#D84135")
            variables.add_column("Variable Adjustment", justify="right", style="#3B629E")
            variables.add_column("Normalized Variable Adjustment", justify="right", style="#3B629E")

            variable_labels = [
                self._find_variable_labels(var)
                for var in self.network.collect_all_iteration_variables()
            ]

            for i, val in enumerate(x):
                label = (
                    "\n".join(variable_labels[i])
                    if i < len(variable_labels)
                    else "<unlabeled>"
                )

                variables.add_row(
                    f"x[{i}]",
                    label,
                    f"{val:.6e}",
                    f"{dx[i]:+.6e}",
                    f"{normalized_variable_adjustment[i]:.3e}",
                )

            self.console.print(variables)

        if r is not None:
            residuals = Table(
                title="Last Available Residuals",
                box=box.SIMPLE,
                show_header=True,
                header_style="bold",
            )

            residuals.add_column("Index", justify="right", style="dim")
            residuals.add_column("Residual", style="#fdf0d5")
            residuals.add_column("Value", justify="right", style="#3B629E")

            residual_labels = self._collect_residual_labels()

            for i, value in enumerate(r):
                label = (
                    residual_labels[i]
                    if i < len(residual_labels)
                    else "<unlabeled>"
                )

                residuals.add_row(
                    f"r[{i}]",
                    label,
                    f"{value:.6e}",
                )

            self.console.print(residuals)

        self.console.print()