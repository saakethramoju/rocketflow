import math
from typing import Callable


class State:
    def __init__(
        self,
        value: float | None = None,
        expr: Callable[[], float] | None = None,
        bounds: tuple[float | None, float | None] | None = None,
    ):
        """
        A scalar state that can either:
        - store a direct numeric value, or
        - be derived from an expression via `expr`.

        This class also supports "placeholder" states via `State()`.
        That is useful when a State must exist up front so multiple
        components can share the same reference, even though its numeric
        value will only be assigned later.

        Parameters
        ----------
        value
            Stored numeric value for a non-derived state.
        expr
            Callable used to compute a derived state's value on demand.
        bounds
            Optional tuple of the form:
                (lower, upper)
                (lower, None)
                (None, upper)
            Bounds are enforced only when assigning to a non-derived state.

        Notes
        -----
        - `State()` is allowed and creates an uninitialized, non-derived state.
        - Accessing `.value` on an uninitialized non-derived state raises ValueError.
        - Derived states cannot be assigned to directly.
        """
        self._expr = expr
        self._lower_bound, self._upper_bound = self._normalize_bounds(bounds)
        self._value = None

        if value is not None:
            value = float(value)
            self._validate_bounds(value)
            self._value = value

    @staticmethod
    def _normalize_bounds(
        bounds: tuple[float | None, float | None] | None,
    ) -> tuple[float | None, float | None]:
        if bounds is None:
            return None, None

        if not isinstance(bounds, tuple) or len(bounds) != 2:
            raise ValueError("bounds must be None or a tuple of the form (lower, upper).")

        lower, upper = bounds

        if lower is not None:
            lower = float(lower)
        if upper is not None:
            upper = float(upper)

        if lower is not None and upper is not None and lower > upper:
            raise ValueError(
                f"Invalid bounds: lower bound {lower} is greater than upper bound {upper}."
            )

        return lower, upper

    def _validate_bounds(self, v: float) -> None:
        if self._lower_bound is not None and v < self._lower_bound:
            raise ValueError(
                f"Value {v} is below the lower bound of {self._lower_bound}."
            )
        if self._upper_bound is not None and v > self._upper_bound:
            raise ValueError(
                f"Value {v} is above the upper bound of {self._upper_bound}."
            )

    @property
    def value(self) -> float:
        """
        Return the current numeric value.

        Raises
        ------
        ValueError
            If this is a non-derived State that has not been assigned yet.
        """
        if self._expr is not None:
            return self._expr()

        if self._value is None:
            raise ValueError(
                "State has no assigned value yet. "
                "Initialize it with State(value=...) or assign state.value = ..."
            )

        return self._value

    @value.setter
    def value(self, v: float) -> None:
        if self._expr is not None:
            raise AttributeError("Cannot assign to a derived State.")

        v = float(v)
        self._validate_bounds(v)
        self._value = v

    @property
    def is_derived(self) -> bool:
        return self._expr is not None

    @property
    def is_assigned(self) -> bool:
        """
        True if this State currently has a usable value.

        - For derived states, this is always True because their value is computed.
        - For non-derived states, this is True only after assignment.
        """
        return self._expr is not None or self._value is not None

    @property
    def bounds(self) -> tuple[float | None, float | None] | None:
        if self._lower_bound is None and self._upper_bound is None:
            return None
        return (self._lower_bound, self._upper_bound)

    @property
    def lower_bound(self) -> float | None:
        return self._lower_bound

    @property
    def upper_bound(self) -> float | None:
        return self._upper_bound

    @property
    def has_bounds(self) -> bool:
        return self._lower_bound is not None or self._upper_bound is not None

    def is_within_bounds(self, v: float | None = None) -> bool:
        """
        Check whether a value lies within bounds.

        If `v` is None, checks the current state value.
        """
        if v is None:
            v = self.value

        if self._lower_bound is not None and v < self._lower_bound:
            return False
        if self._upper_bound is not None and v > self._upper_bound:
            return False
        return True

    def __str__(self) -> str:
        if self.is_derived:
            value_str = str(self.value)
        elif self._value is None:
            value_str = "<uninitialized>"
        else:
            value_str = str(self._value)

        if self.has_bounds:
            return f"State(value={value_str}, bounds={self.bounds})"
        return f"State(value={value_str})"

    def __repr__(self) -> str:
        if self.is_derived:
            value_str = str(self.value)
        elif self._value is None:
            value_str = "<uninitialized>"
        else:
            value_str = str(self._value)

        if self.has_bounds:
            return f"State({value_str}, bounds={self.bounds})"
        return f"State({value_str})"

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