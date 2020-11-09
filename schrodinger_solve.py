import numpy as np
from numpy import polyval
from scipy.special import genlaguerre
import matplotlib.pyplot as plt


class SchrodingerEquation:
    """
        takes in a potential, and (if any) values for m and hbar
        to construct an object containing information about the 
        time independent Schrodinger Equation obtained from 
        these information
    """
    def __init__(self, potential, m=1, hbar=1, min_x=-4, max_x=4, intervals=1000):
        self.x_values = np.linspace(min_x, max_x, intervals)
        self.potential_values = list(map(potential, self.x_values))
        self.discrete_potential_matrix = np.diag(self.potential_values)
        self.dimension = intervals
        self.ddx2_matrix = self.ddx2(intervals)
        self.hamiltonian = - (hbar**2 / (2 * m)) * self.ddx2_matrix + self.discrete_potential_matrix
        self.eigensolution = self.solve_eigenproblem()
        #todo

    def ddx2(self, N):
        """
            provides a discrete approximation to the second derivative in 
            one dimension for a function approximated by a N-dimendional vector
        """
        delta = self.x_values[0] - self.x_values[1]
        ddx2_matrix = (-2 * np.diag(np.ones(N, np.float32), 0) + np.diag(np.ones(N-1, np.float32), -1) + np.diag(np.ones(N-1, np.float32), 1))/(delta**2)
        return ddx2_matrix

    def solve_eigenproblem(self):
        """
            uses the discrete Hamiltonian matrix obtained from the constructor
            to compute the eigenvectors and eigenvalues of it
        """
        e_vals, e_vecs = np.linalg.eig(self.hamiltonian)
        sorted_e_vals_indices = np.argsort(e_vals)
        return e_vals[sorted_e_vals_indices], e_vecs[:, sorted_e_vals_indices]
    
    def superimpose(self, states):
        """
            takes in an array of states to superimpose, and returns a vector 
            that approximates the (normalized) superimposed state

            states is an array of length at least one
        """
        ans = np.zeros(self.dimension)
        for i in states:
            ans += self.eigensolution[1][:, i - 1]
        return ans/np.linalg.norm(ans)


#test code
test_equation = SchrodingerEquation(lambda x: 20 * np.heaviside(x, 1)) # omega=1, m=1, hbar=1

# print("\nThe discrete apprixmation to V(x):\n", test_equation.discrete_potential_matrix)
# print("\nThe discrete approximation to ddx2:\n", test_equation.ddx2_matrix)


#The first ten eigenvalues are relatively nicely approximated.. after that the linear growth in the approximated eigenvalues gets dominated
# for n in range(1, 40): # the _th eigenvalue
#     print("\nThe " + str(n) + "th energy eigenvalue for this setup is:", len(test_equation.eigensolution[0][n-1]), "\nThe analytic solution for this setup is:", (n - 1) + 0.5)
#test_superimposed = test_equation.superimpose([1, 2, 3])

plt.figure()
plt.xlabel("$x$")
plt.ylabel("$\psi(x)$")
plt.title("Wavefunction vs. Position")
# plt.plot(test_equation.x_values, test_equation.eigensolution[1][:, 4], linestyle="--")
# plt.plot(test_equation.x_values, test_equation.eigensolution[1][:, 6], linestyle="--")
plt.plot(test_equation.x_values, test_equation.eigensolution[1][:, 10], linestyle="--")
plt.plot(test_equation.x_values, test_equation.potential_values/np.linalg.norm(test_equation.potential_values), linestyle=":", label="potential")
#plt.plot(test_equation.x_values, test_superimposed, label="Superimoposed (0, 1, 2) State")
plt.legend()
plt.show()