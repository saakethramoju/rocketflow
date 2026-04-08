import math
from typing import Callable


class State:
    def __init__(self, value: float | None = None, expr: Callable[[], float] | None = None):
        self._value = value
        self._expr = expr

    @property
    def value(self) -> float:
        if self._expr is not None:
            return self._expr()
        if self._value is None:
            raise ValueError("State has neither a stored value nor an expression.")
        return self._value
        
    @property
    def is_derived(self) -> bool:
        return self._expr is not None

    @value.setter
    def value(self, v: float) -> None:
        if self._expr is not None:
            raise AttributeError("Cannot assign to a derived State.")
        self._value = v

    def __str__(self) -> str:
        return f"State(value={self.value})"

    def __repr__(self) -> str:
        return f"State({self.value})"

    @staticmethod
    def _coerce(other) -> "State":
        if isinstance(other, State):
            return other
        return State(float(other))

    def __add__(self, other):
        other = self._coerce(other)
        return State(expr=lambda: self.value + other.value)

    def __radd__(self, other):
        other = self._coerce(other)
        return State(expr=lambda: other.value + self.value)

    def __sub__(self, other):
        other = self._coerce(other)
        return State(expr=lambda: self.value - other.value)

    def __rsub__(self, other):
        other = self._coerce(other)
        return State(expr=lambda: other.value - self.value)

    def __mul__(self, other):
        other = self._coerce(other)
        return State(expr=lambda: self.value * other.value)

    def __rmul__(self, other):
        other = self._coerce(other)
        return State(expr=lambda: other.value * self.value)

    def __truediv__(self, other):
        other = self._coerce(other)
        return State(expr=lambda: self.value / other.value)

    def __rtruediv__(self, other):
        other = self._coerce(other)
        return State(expr=lambda: other.value / self.value)

    def __pow__(self, other):
        other = self._coerce(other)
        return State(expr=lambda: self.value ** other.value)

    def __rpow__(self, other):
        other = self._coerce(other)
        return State(expr=lambda: other.value ** self.value)

    def __neg__(self):
        return State(expr=lambda: -self.value)

    def __abs__(self):
        return State(expr=lambda: abs(self.value))

    def sqrt(self):
        return State(expr=lambda: math.sqrt(self.value))