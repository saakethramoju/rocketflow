from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from System import Component, State, Balance


class Network:

    def __init__(self, name: str):
        self.name = name
        self.component_list = []
        self.balance_list = []

    def add_component(self, component: Component) -> None:
        self.component_list.append(component)

    def add_balance(self, balance: Balance) -> None:
        self.balance_list.append(balance)

    @property
    def components(self) -> list[str]:
        return [c.name for c in self.component_list]

    @property
    def balances(self) -> list[str]:
        return [b.name for b in self.balance_list]


    # -------------- STEADY-STATE -------------- #

    @property
    def lower_bounds(self) -> list[float]:
        return [var.lower_bound for var in self.collect_all_iteration_variables()]
    
    @property
    def upper_bounds(self) -> list[float]:
        return [var.upper_bound for var in self.collect_all_iteration_variables()]
    
    @property
    def has_bounds(self) -> bool:
        return any([var.has_bounds for var in self.collect_all_iteration_variables()])
    
    @property
    def keep_feasible(self) -> list[float]:
        return [var.keep_feasible for var in self.collect_all_iteration_variables()]

            

    @property
    def iteration_variables(self) -> str:
        names = []

        def find_name(obj, target):
            for k, v in obj.__dict__.items():
                if v is target:
                    return f"{obj.name}:{k}"
            return f"{obj.name}:<unknown>"

        # components
        for comp in self.component_list:
            for var in comp.iteration_variables:
                names.append(find_name(comp, var))

        # balances
        for bal in self.balance_list:
            for var in bal.iteration_variables:
                names.append(find_name(bal, var))

        return "\n".join(names)
            

    @property
    def iteration_values(self) -> list[float]:
        self._validate_no_iteration_overlap()
        return [var.value for var in self.collect_all_iteration_variables()]


    @property
    def residuals(self) -> list[float]:
        residuals = []

        for comp in self.component_list:
            try:
                comp_residuals = comp.residuals
            except Exception as e:
                original_msg = str(e).splitlines()[0]

                raise RuntimeError(
                    f"Failed while evaluating residuals for component `{comp.name}` "
                    f"of type `{type(comp).__name__}`.\n\n"
                    "A State used inside this component's residual equations is probably "
                    "unassigned.\n\n"
                    "Likely fixes:\n"
                    "  - Give the missing State an initial value\n"
                    "  - Connect it to another component that computes it\n"
                    "  - Make it an iteration variable\n"
                    "  - If this is a static/transient-only quantity, do not use it in "
                    "steady-state residuals\n\n"
                    f"Original error: {original_msg}"
                ) from None

            if isinstance(comp_residuals, (list, tuple)):
                residuals.extend(comp_residuals)
            else:
                residuals.append(comp_residuals)

        balance_residuals = [
            r for balance in self.balance_list for r in balance.residuals
        ]

        residuals.extend(balance_residuals)

        return residuals



    def pre_evaluation(self) -> None:
        for c in self.component_list:
            c.pre_evaluation()


    def evaluate_states(
        self,
        max_passes: int = 20,
        tolerance: float = 1e-10,
    ) -> None:
        """
        Evaluate all component states until explicit state updates stop changing.

        This removes dependence on component ordering. For example:

            composition -> density -> mass_flow

        may require several passes to propagate through the network.

        Iteration variables are frozen during settling so components cannot
        accidentally overwrite solver-owned states (pressure, enthalpy,
        balance variables, etc.).

        Examples
        --------
        FluidLookup -> Darcy -> Volume

        or

        Composition -> FluidLookup -> FlowSplitter -> Darcy

        may require multiple passes before all derived states are consistent.

        """

        # Save current solver variables.
        iteration_snapshot = self._snapshot_iteration_variables()

        for _ in range(max_passes):

            # Record all non-iteration state values before this pass.
            old_values = self._collect_state_values()

            for c in self.component_list:

                # Protect solver-owned variables.
                self._restore_iteration_variables(iteration_snapshot)

                c.evaluate_states()

                # Protect solver-owned variables again in case a component
                # modified them.
                self._restore_iteration_variables(iteration_snapshot)

            # Record all non-iteration state values after this pass.
            new_values = self._collect_state_values()

            # Stop once explicit state updates have settled.
            if self._max_state_change(old_values, new_values) < tolerance:
                return
                    
    def _snapshot_iteration_variables(self) -> dict[int, float]:
        """
        Store the current values of all iteration variables.

        These values belong to the nonlinear solver and should not be
        modified by evaluate_states() during settling.
        """

        snapshot = {}

        for var in self.collect_all_iteration_variables():
            if var.is_assigned:
                snapshot[id(var)] = float(var.value)

        return snapshot


    def _restore_iteration_variables(
        self,
        snapshot: dict[int, float],
    ) -> None:
        """
        Restore solver-owned iteration variables to their original values.

        This prevents component evaluation from overwriting the solver's
        current iterate.
        """

        for var in self.collect_all_iteration_variables():
            key = id(var)

            if key in snapshot:
                var.value = snapshot[key]
            

    def _collect_state_values(self) -> dict[int, float]:
        """
        Collect all non-iteration State values currently assigned in the
        network.

        These values are used only to determine whether explicit state
        propagation has converged.

        Iteration variables are intentionally ignored because they are
        controlled by the nonlinear solver, not by state-settling passes.
        """

        values = {}

        iteration_ids = {
            id(var)
            for var in self.collect_all_iteration_variables()
        }

        def collect(attr_value):

            # Composition container.
            if hasattr(attr_value, "fraction"):
                for _, state in attr_value:

                    # Ignore solver variables.
                    if id(state) in iteration_ids:
                        continue

                    if state.is_assigned:
                        try:
                            values[id(state)] = float(state.value)
                        except Exception:
                            pass
                return

            # Ordinary State.
            if hasattr(attr_value, "is_assigned"):

                # Ignore solver variables.
                if id(attr_value) in iteration_ids:
                    return

                if attr_value.is_assigned:
                    try:
                        values[id(attr_value)] = float(attr_value.value)
                    except Exception:
                        pass

        for comp in self.component_list:
            for attr_value in comp.__dict__.values():
                collect(attr_value)

        for bal in self.balance_list:
            for attr_value in bal.__dict__.values():
                collect(attr_value)

        return values
    


    def _max_state_change(
        self,
        old: dict[int, float],
        new: dict[int, float],
    ) -> float:
        """
        Compute the largest normalized change between two state snapshots.

        Used to determine whether explicit state propagation has settled.

        A value near zero means another evaluation pass would produce
        essentially the same network state.
        """

        max_change = 0.0

        for key, new_value in new.items():
            old_value = old.get(key)

            # Newly-created state.
            if old_value is None:
                max_change = max(max_change, abs(new_value))
                continue

            # Relative change with protection against division by small values.
            scale = max(abs(new_value), 1.0)

            max_change = max(
                max_change,
                abs(new_value - old_value) / scale,
            )

        return max_change


    def collect_iteration_variables(self) -> list[State]:
        iter_vars = []
        for comp in self.component_list:
            iter_vars.extend(comp.iteration_variables)

        return iter_vars
    
    def collect_balance_iteration_variables(self) -> list[State]:
        iter_vars = []
        for bal in self.balance_list:
            iter_vars.extend(bal.iteration_variables)

        return iter_vars
    
    def collect_all_iteration_variables(self) -> list[State]:
        return (self.collect_iteration_variables() + self.collect_balance_iteration_variables())


    def assign_iteration_values(self, iteration_values: list[float]) -> None:
        iter_var_list = self.collect_all_iteration_variables()
        if len(iteration_values) != len(iter_var_list):
            raise ValueError(
                f"Length mismatch: got {len(iteration_values)} iteration values "
                f"but expected {len(iter_var_list)}"
            )
        for val, var in zip(iteration_values, iter_var_list):
            var.value = val

    

    def _validate_no_iteration_overlap(self) -> None:
        comp_ids = {id(v) for v in self.collect_iteration_variables()}
        bal_ids = {id(v) for v in self.collect_balance_iteration_variables()}
        overlap_ids = comp_ids & bal_ids

        if overlap_ids:
            raise ValueError(self._format_iteration_overlap_error(overlap_ids))
    

    def _format_iteration_overlap_error(self, overlap_ids: set[int]) -> str:
        component_names = []
        balance_names = []

        for comp in self.component_list:
            comp_iter_ids = {id(v) for v in comp.iteration_variables}
            shared_ids = comp_iter_ids & overlap_ids
            if not shared_ids:
                continue

            for attr_name, attr_value in comp.__dict__.items():
                if id(attr_value) in shared_ids:
                    component_names.append(f"{comp.name}:{attr_name}")

        for bal in self.balance_list:
            bal_iter_ids = {id(v) for v in bal.iteration_variables}
            shared_ids = bal_iter_ids & overlap_ids
            if not shared_ids:
                continue

            # show both the balance name and its variable attribute
            for attr_name, attr_value in bal.__dict__.items():
                if id(attr_value) in shared_ids:
                    balance_names.append(f"{bal.name}:{attr_name}")

        component_names = sorted(set(component_names))
        balance_names = sorted(set(balance_names))

        lines = [
            "Iteration variable overlap detected.",
            "",
            "Balance iteration variables cannot be the same as component iteration variables.",
            "A State used as a Balance solve variable must not also appear in a Component iteration_variables list.",
            "",
            "Overlapping component variables:",
            *[f"  - {name}" for name in component_names],
            "",
            "Conflicting balance variables:",
            *[f"  - {name}" for name in balance_names],
        ]

        return "\n".join(lines)

    # -------------- SOLUTION EXPORTING -------------- #

    def __str__(self) -> str:
        lines = []

        lines.append(f"Network: {self.name}")
        lines.append(f"Components ({len(self.component_list)}):")

        for comp in self.component_list:
            lines.append(
                f"  ├─ [{comp.__class__.__name__}]: {comp.name}"
            )

        return "\n".join(lines)


    def save(
        self,
        filename: str | None = None,
        return_type: str = "dict",
    ):
        return_type = return_type.lower()

        # ---------- build base data ----------
        records = []
        for comp in self.component_list:
            ignored_attributes = {"name", "network"} | comp.ignored_export_attributes

            for attr_name, attr_value in comp.__dict__.items():
                if attr_name in ignored_attributes or attr_name.startswith("_"):
                    continue

                if hasattr(attr_value, "fraction"):
                    if attr_value.is_assigned:
                        for species, state in attr_value:
                            try:
                                value = state.value
                            except Exception:
                                value = "<unavailable>"

                            records.append({
                                "component_name": comp.name,
                                "component_type": comp.__class__.__name__,
                                "attribute": f"{attr_name}.{species}",
                                "value": value,
                            })
                    else:
                        records.append({
                            "component_name": comp.name,
                            "component_type": comp.__class__.__name__,
                            "attribute": attr_name,
                            "value": "<uninitialized>",
                        })

                    continue

                if hasattr(attr_value, "is_assigned"):
                    if attr_value.is_assigned:
                        try:
                            value = attr_value.value
                        except Exception:
                            value = "<unavailable>"
                    else:
                        value = "<uninitialized>"
                else:
                    value = attr_value

                records.append({
                    "component_name": comp.name,
                    "component_type": comp.__class__.__name__,
                    "attribute": attr_name,
                    "value": value,
                })

        # ---------- balances ----------
        for bal in self.balance_list:
            ignored_attributes = {"name", "network"}

            for attr_name, attr_value in bal.__dict__.items():
                if attr_name in ignored_attributes or attr_name.startswith("_"):
                    continue

                if hasattr(attr_value, "fraction"):
                    if attr_value.is_assigned:
                        for species, state in attr_value:
                            records.append({
                                "component_name": bal.name,
                                "component_type": bal.__class__.__name__,
                                "attribute": f"{attr_name}.{species}",
                                "value": state.value,
                            })
                    else:
                        records.append({
                            "component_name": bal.name,
                            "component_type": bal.__class__.__name__,
                            "attribute": attr_name,
                            "value": "<uninitialized>",
                        })

                    continue

                if hasattr(attr_value, "is_assigned"):
                    if attr_value.is_assigned:
                        try:
                            value = attr_value.value
                        except Exception:
                            value = "<unavailable>"
                    else:
                        value = "<uninitialized>"
                else:
                    value = attr_value

                records.append({
                    "component_name": bal.name,
                    "component_type": bal.__class__.__name__,
                    "attribute": attr_name,
                    "value": value,
                })

        # ---------- return object ----------
        if return_type == "dict":
            result = records

        elif return_type == "dataframe":
            import pandas as pd
            result = pd.DataFrame(records)

        else:
            raise ValueError("return_type must be 'dict' or 'dataframe'")

        # ---------- file export ----------
        if filename is not None:
            ext = filename.split(".")[-1].lower()

            # ensure we have dataframe for file writing
            import pandas as pd
            df = pd.DataFrame(records)

            if ext == "csv":
                df.to_csv(filename, index=False)

            elif ext == "json":
                import json
                with open(filename, "w") as f:
                    json.dump(records, f, indent=4)

            elif ext in {"xlsx", "xls"}:
                df.to_excel(filename, index=False)

            else:
                raise ValueError(
                    "Unsupported file extension. Use .csv, .json, or .xlsx"
                )

        return result