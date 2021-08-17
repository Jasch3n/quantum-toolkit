import numpy as np
from numpy import polyval
from scipy.special import genlaguerre
import matplotlib.pyplot as plt

from schrodinger_solve import *

# def associated_laguerre(alpha, n):
#     """Given alpha and n, return a ordpoly representing 
#         the associated laguerre polynomial L^(alpha)_n

#     Args:
#         alpha ([int]): [> -1]
#         n ([int]): [> -1]

#     Returns:
#         [ordpoly]: [a numpy polynomial representating an associated laguerre polynomial]
#     """
#     if (alpha < 0) or (n < 0):
#         raise("invalid input for p and q")
#     return genlaguerre(n, alpha)

# test_laguerre = associated_laguerre(2, 10)
# domain = np.linspace(-100, 100, 100)
# lag_values = list(map(lambda x: polyval(test_laguerre, x), domain))

# plt.figure()
# plt.plot(domain, lag_values)
# plt.plot()
# plt.show()

def free_particle_pot(x):
    return 0

def harmonic_oscillator_pot(x, omega=1):
    return 0.5 * omega**2 * x ** 2

def step_potential(x, height=20):
    return height * np.heaviside(x, 1)

def random_potential(x):
    return (2.71 ** x)

# step_potential_eq = SchrodingerEquation(step_potential)
# step_potential_eq.graph_nth_eigenstate(4)
# step_potential_eq.display_graph()

test_eq = SchrodingerEquation(harmonic_oscillator_pot)
test_eq.graph_nth_eigenstate(3)
test_eq.graph_superimposed([1,2,3,4,5,6,7,8,9])
# test_equation.display_graph()
wavefunction_time_evolution_anim(test_eq, 5)