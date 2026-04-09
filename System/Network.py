from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from System import Component, Variable


class Network:

    def __init__(self, name: str):
        self.name = name
        self.component_list = []
        self._iter_var_cache = None # iteration variable list[str]
        self._iter_var_name_cache = None # iteration variable list[Variable]

    def _clear_cache(self) -> None:
        self._iter_var_cache = None
        self._iter_var_name_cache = None

    def add_component(self, component: Component) -> None:
        self.component_list.append(component)
        self._clear_cache()

    @property
    def components(self) -> list[str]:
        return [c.name for c in self.component_list]



    # -------------- STEADY-STATE -------------- #

    @property
    def lower_bounds(self) -> list[float]:
        return [var.lower_bound for var in self.collect_iteration_variables()]
    
    @property
    def upper_bounds(self) -> list[float]:
        return [var.upper_bound for var in self.collect_iteration_variables()]
    
    @property
    def keep_feasible(self) -> list[bool]:
        return [var.keep_feasible for var in self.collect_iteration_variables()]
    
    @property
    def iteration_variables(self) -> list[str]:
        if self._iter_var_name_cache is None:
            names = []
            for comp in self.component_list:
                comp_iter_vars = comp.iteration_variables
                comp_iter_var_ids = {id(var) for var in comp_iter_vars}

                for attr_name, attr_value in comp.__dict__.items():
                    if id(attr_value) in comp_iter_var_ids:
                        names.append(f"{comp.name}:{attr_name}")

            self._iter_var_name_cache = names

        return self._iter_var_name_cache
    

    @property
    def iteration_values(self) -> list[float]:
        return [var.value for var in self.collect_iteration_variables()]
    

    @property
    def residuals(self) -> list[float]:
        residuals = []
        for comp in self.component_list:
            residuals.extend(comp.residuals)
        return residuals


    def pre_evaluation(self) -> None:
        for c in self.component_list:
            c.pre_evaluation()

    def evaluate_states(self) -> None:
        for c in self.component_list:
            c.evaluate_states()


    def collect_iteration_variables(self) -> list[Variable]:
        if self._iter_var_cache is None:
            iter_vars = []
            for comp in self.component_list:
                iter_vars.extend(comp.iteration_variables)
            self._iter_var_cache = iter_vars

        return self._iter_var_cache


    def assign_iteration_values(self, iteration_values: list[float]) -> None:
        iter_var_list = self.collect_iteration_variables()
        if len(iteration_values) != len(iter_var_list):
            raise ValueError(
                f"Length mismatch: got {len(iteration_values)} iteration values "
                f"but expected {len(iter_var_list)}"
            )
        for val, var in zip(iteration_values, iter_var_list):
            var.value = val


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

                # pretty formatting for Variables vs other types
                if hasattr(attr_value, "value"):
                    val_str = f"{attr_value.value:.4g}"
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