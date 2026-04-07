class Variable:
    def __init__(self, obj):
        if hasattr(obj, "value"):
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