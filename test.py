from System import *
from Solvers import *

from constants import *


MixtureNetwork = Network("Mixture Flow")

SourceFluid = FluidLookup(
    "Source Fluid",
    MixtureNetwork,
    {"gn2": 0.75, "O2": 0.25},
    #"o2",
    pressure=3e5,
    temperature=300,
)

D = 6 * IN_TO_M
A = (np.pi/4) * D**2



'''
solution = SteadyState(MixtureNetwork).solve(
    return_type="dataframe",
    verbose=True,
    static=False,
    print_solution=True,
    filename='test.xlsx'
)
'''


c1 = Composition({"N2": 0.75, "O2": 0.25})
c2 = Composition({"Ar": 1.0})

c3 = Composition({"N2": 0.5, "O2": 0.25, "Ar": 0.25})
c4 = Composition({"N2": 1.0})

print((c1 | c2) <= set(c4.species))
print((c1 | c2) - set(c4.species))