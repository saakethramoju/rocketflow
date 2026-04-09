import numpy as np


class Variable:
    def __init__(self, obj):
        if hasattr(type(obj), "value"):
            self._state = obj
            self._is_state = True
        elif isinstance(obj, (int, float)):
            self._value = float(obj)
            self._is_state = False
        else:
            raise TypeError(
                f"Variable must wrap a State or float, got {type(obj)}"
            )

    @property
    def value(self) -> float:
        if self._is_state:
            return self._state.value
        return self._value

    @value.setter
    def value(self, v: float) -> None:
        if self._is_state:
            self._state.value = v
        else:
            self._value = v

    @property
    def lower_bound(self) -> float:
        if self._is_state and hasattr(self._state, "lower_bound"):
            return self._state.lower_bound
        return -np.inf

    @property
    def upper_bound(self) -> float:
        if self._is_state and hasattr(self._state, "upper_bound"):
            return self._state.upper_bound
        return np.inf

    @property
    def bounds(self) -> tuple[float, float]:
        return (self.lower_bound, self.upper_bound)

    @property
    def has_bounds(self) -> bool:
        return not (self.lower_bound == -np.inf and self.upper_bound == np.inf)

    @property
    def keep_feasible(self) -> bool:
        if self._is_state and hasattr(self._state, "keep_feasible"):
            return self._state.keep_feasible
        return False