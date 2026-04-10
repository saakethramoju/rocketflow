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
        return (
            [r for comp in self.component_list for r in comp.residuals]
            + [r for bal in self.balance_list for r in bal.residuals]
        )


    def pre_evaluation(self) -> None:
        for c in self.component_list:
            c.pre_evaluation()

    def evaluate_states(self) -> None:
        for c in self.component_list:
            c.evaluate_states()


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
            lines.append(f"  ├─ [{comp.__class__.__name__}]: {comp.name}")

            for attr_name, attr_value in comp.__dict__.items():
                # skip noisy/internal attributes
                if attr_name in {"name", "network"} or attr_name.startswith("_"):
                    continue

                # pretty formatting for State-like objects vs other types
                if hasattr(type(attr_value), "value"):
                    try:
                        val = attr_value.value
                        val_str = f"{val:.4g}"
                    except (ValueError, AttributeError):
                        val_str = "—"
                else:
                    val_str = str(attr_value)

                lines.append(f"  │   {attr_name:<12}: {val_str}")

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
            for attr_name, attr_value in comp.__dict__.items():
                if attr_name in {"name", "network"} or attr_name.startswith("_"):
                    continue

                if hasattr(attr_value, "value"):
                    value = attr_value.value
                else:
                    value = attr_value

                records.append({
                    "component_name": comp.name,
                    "component_type": comp.__class__.__name__,
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