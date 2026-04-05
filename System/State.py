class State:
    def __init__(self, 
                 value: float):
        self._value = value

    @property
    def value(self) -> float:
        return self._value

    @value.setter
    def value(self, v: float) -> None:
        self._value = v

    def __str__(self) -> str:
        return f"State: {self._value}"