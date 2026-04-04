from Physics import NodeMethod, BranchMethod

class Component:
    
    def __init__(self, 
                 name: str):
        self.name = name


class Node(Component):
    def __init__(self, 
                 name: str,
                 method: NodeMethod):
        super().__init__(name)


class Branch(Component):
    def __init__(self, 
                 name: str,
                 method: BranchMethod):
        super().__init__(name)
