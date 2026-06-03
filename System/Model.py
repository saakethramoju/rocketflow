class ModelOption:
    """
    Deferred component definition used by Model.

    A ModelOption stores the component class and constructor
    arguments needed to build a real component later.

    Unlike a normal Component, a ModelOption does not register
    itself with a Network when created.
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
        Construct and return the real component represented by this option.
        """

        return self.component_class(
            component_name,
            network,
            **self.kwargs,
        )

    @property
    def component_name(self) -> str:
        """
        Name of the component class this option will build.
        """

        return self.component_class.__name__

    def __str__(self) -> str:
        return self.name

    def __repr__(self) -> str:
        return (
            f"ModelOption("
            f"name={self.name!r}, "
            f"component_class={self.component_class.__name__})"
        )


class Model:
    """
    Collection of alternative component implementations.

    A Model stores one or more ModelOptions and builds one selected
    option into the network when build() is called.

    Notes
    -----
    1) Model does not build automatically during initialization.
    2) Only the active option is converted into a real component.
    3) Switching options should happen between solve attempts, not during
       a Newton iteration.
    """

    def __init__(
        self,
        name: str,
        network,
        components: list[ModelOption],
        order: list[str] | None = None,
    ):
        """
        Parameters
        ----------
        name:
            Name assigned to the real component when it is built.

        network:
            Network the selected component will be added to.

        components:
            Available ModelOptions.

        order:
            Optional list of option names defining the try order.
            Defaults to the declaration order in components.
        """

        self.name = name
        self.network = network

        # Store options by user-facing option name.
        self.components = {
            component.name: component
            for component in components
        }

        # Default to the same order the options were provided in.
        self.order = (
            order
            if order is not None
            else [component.name for component in components]
        )

        # Validate early so typos fail before solving.
        self._validate()

        self.active_option_name = None
        self.active_component = None

    def _validate(self) -> None:
        """
        Validate component names and order entries.
        """

        if len(self.components) == 0:
            raise ValueError(f"{self.name}: Model requires at least one option.")

        if len(self.components) != len(set(self.components)):
            raise ValueError(f"{self.name}: duplicate model option names found.")

        invalid_options = [
            option_name
            for option_name in self.order
            if option_name not in self.components
        ]

        if invalid_options:
            raise ValueError(
                f"{self.name}: order contains invalid options "
                f"{invalid_options}. Valid options are "
                f"{list(self.components)}."
            )

    def build(
        self,
        option_name: str | None = None,
    ):
        """
        Build and return the selected component.

        If no option name is supplied, the first option in self.order is used.
        """

        if self.active_component is not None:
            raise RuntimeError(
                f"{self.name}: model has already built "
                f"option {self.active_option_name!r}."
            )

        # Use the first option by default.
        option_name = option_name or self.order[0]

        if option_name not in self.components:
            raise ValueError(
                f"{self.name}: unknown model option {option_name!r}. "
                f"Valid options are {list(self.components)}."
            )

        option = self.components[option_name]

        self.active_option_name = option_name

        self.active_component = option.build(
            self.name,
            self.network,
        )

        return self.active_component

    def next(self) -> str:
        """
        Return the name of the next model option.

        This does not build the next option. It only tells you what the
        next option would be.
        """

        if self.active_option_name is None:
            return self.order[0]

        current_index = self.order.index(self.active_option_name)
        next_index = current_index + 1

        if next_index >= len(self.order):
            raise RuntimeError(f"{self.name}: no remaining model options.")

        return self.order[next_index]

    def build_next(self):
        """
        Build the next model option in the order list.
        """

        return self.build(self.next())

    @property
    def active_option(self):
        """
        Currently selected ModelOption.
        """

        if self.active_option_name is None:
            return None

        return self.components[self.active_option_name]
    

    def clear(self) -> None:
        """
        Remove the currently active component from the network.

        After clearing, this Model can build another option.
        """

        if self.active_component is None:
            return

        self.network.remove_component(self.active_component)

        self.active_component = None
        self.active_option_name = None



    def replace(
        self,
        option_name: str,
    ):
        """
        Replace the active component with another model option.
        """

        self.clear()
        return self.build(option_name)



    def build_next(self):
        """
        Replace the active component with the next option in the order list.
        """

        next_option_name = self.next()
        return self.replace(next_option_name)



    @property
    def available_options(self) -> list[str]:
        """
        Available option names.
        """

        return list(self.components)
    

    @property
    def has_next(self) -> bool:
        """
        True if another model option remains.
        """

        if self.active_option_name is None:
            return len(self.order) > 0

        current_index = self.order.index(self.active_option_name)

        return current_index < len(self.order) - 1


    def __str__(self) -> str:
        return self.name

    def __repr__(self) -> str:
        return (
            f"Model("
            f"name={self.name!r}, "
            f"options={list(self.components)}, "
            f"order={self.order}, "
            f"active={self.active_option_name!r})"
        )