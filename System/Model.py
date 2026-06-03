class ModelOption:
    """
    Deferred component definition used by Model.

    A ModelOption stores the component type and constructor
    arguments required to build a component later.

    Unlike a normal Component, a ModelOption does not register
    itself with a Network when created.

    Examples
    --------
    >>> Darcy = DischargeCoefficient.model(
    ...     "darcy",
    ...     upstream_pressure=P1,
    ...     downstream_pressure=P2,
    ... )

    >>> Darcy.build("Outlet 1", network)
    """

    def __init__(
        self,
        name: str,
        component_class: type,
        kwargs: dict,
    ):
        self.name = name
        self.component_class = component_class
        self.kwargs = kwargs

    def build(
        self,
        component_name: str,
        network,
    ):
        """
        Construct and return the component represented
        by this ModelOption.
        """
        return self.component_class(
            component_name,
            network,
            **self.kwargs,
        )
        
    @property
    def component_name(self) -> str:
        return self.component_class.__name__

    def __str__(self) -> str:
        return self.name

    def __repr__(self) -> str:
        return (
            f"ModelOption("
            f"name='{self.name}', "
            f"component_class={self.component_class.__name__})"
        )