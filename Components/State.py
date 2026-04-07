class State:
    def __init__(self, value: float):
        self.value = value

    def __str__(self) -> str:
        return f"State (value={self.value})"