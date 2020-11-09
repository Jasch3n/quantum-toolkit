import numpy as np
from numpy import polyval
from scipy.special import genlaguerre
import matplotlib.pyplot as plt

def associated_laguerre(alpha, n):
    """Given alpha and n, return a ordpoly representing 
        the associated laguerre polynomial L^(alpha)_n

    Args:
        alpha ([int]): [> -1]
        n ([int]): [> -1]

    Returns:
        [ordpoly]: [a numpy polynomial representating an associated laguerre polynomial]
    """
    if (alpha < 0) or (n < 0):
        raise("invalid input for p and q")
    return genlaguerre(n, alpha)

test_laguerre = associated_laguerre(2, 10)
domain = np.linspace(-100, 100, 100)
lag_values = list(map(lambda x: polyval(test_laguerre, x), domain))

plt.figure()
plt.plot(domain, lag_values)
plt.plot()
plt.show()