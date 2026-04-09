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

    # ---------- arithmetic ----------
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

    # ---------- basic math ----------
    def sqrt(self):
        return State(expr=lambda: math.sqrt(self.value))

    def exp(self):
        return State(expr=lambda: math.exp(self.value))

    def expm1(self):
        return State(expr=lambda: math.expm1(self.value))

    def log(self):
        return State(expr=lambda: math.log(self.value))

    def log10(self):
        return State(expr=lambda: math.log10(self.value))

    def log2(self):
        return State(expr=lambda: math.log2(self.value))

    def log1p(self):
        return State(expr=lambda: math.log1p(self.value))

    # ---------- trig ----------
    def sin(self):
        return State(expr=lambda: math.sin(self.value))

    def cos(self):
        return State(expr=lambda: math.cos(self.value))

    def tan(self):
        return State(expr=lambda: math.tan(self.value))

    def asin(self):
        return State(expr=lambda: math.asin(self.value))

    def acos(self):
        return State(expr=lambda: math.acos(self.value))

    def atan(self):
        return State(expr=lambda: math.atan(self.value))

    # ---------- hyperbolic ----------
    def sinh(self):
        return State(expr=lambda: math.sinh(self.value))

    def cosh(self):
        return State(expr=lambda: math.cosh(self.value))

    def tanh(self):
        return State(expr=lambda: math.tanh(self.value))

    def asinh(self):
        return State(expr=lambda: math.asinh(self.value))

    def acosh(self):
        return State(expr=lambda: math.acosh(self.value))

    def atanh(self):
        return State(expr=lambda: math.atanh(self.value))

    # ---------- angle conversions ----------
    def degrees(self):
        return State(expr=lambda: math.degrees(self.value))

    def radians(self):
        return State(expr=lambda: math.radians(self.value))

    # ---------- rounding / integer ----------
    def floor(self):
        return State(expr=lambda: math.floor(self.value))

    def ceil(self):
        return State(expr=lambda: math.ceil(self.value))

    def trunc(self):
        return State(expr=lambda: math.trunc(self.value))

    # ---------- misc ----------
    def modf(self):
        return (
            State(expr=lambda: math.modf(self.value)[0]),
            State(expr=lambda: math.modf(self.value)[1]),
        )

    def fmod(self, other):
        other = self._coerce(other)
        return State(expr=lambda: math.fmod(self.value, other.value))

    def hypot(self, other):
        other = self._coerce(other)
        return State(expr=lambda: math.hypot(self.value, other.value))

    def copysign(self, other):
        other = self._coerce(other)
        return State(expr=lambda: math.copysign(self.value, other.value))