import numpy as np
from numpy import polyval
from scipy.special import genlaguerre
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib import animation

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
        self.hamiltonian = - (hbar ** 2 / (2 * m)) * self.ddx2_matrix + self.discrete_potential_matrix
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
            to compute the eigenvectors and eigenvalues of it and sort them by
            ascending order
        """
        e_vals, e_vecs = np.linalg.eig(self.hamiltonian)
        sorted_e_vals_indices = np.argsort(e_vals)
        return e_vals[sorted_e_vals_indices], e_vecs[:, sorted_e_vals_indices]

    def graph_nth_eigenstate(self, n):
        """
            Given an integer n, graphs the n-th energy eigenstate
            along with the potential

        Args:
            n ([int]): [index of energy eigenstate, with 1 being the first]
        """
        plt.figure()
        plt.xlabel("$x$")
        plt.ylabel("$\psi(x)$")
        title = "Wavefunciton vs. Position"
        plt.title(title)
        plt.plot(self.x_values, self.eigensolution[1][:, n - 1], label="$\psi_{"+ str(n) + "}$")
        plt.plot(self.x_values, self.potential_values/np.linalg.norm(self.potential_values), label="V(x)")
        plt.legend()
    
    def graph_superimposed(self, states):
        """
            takes in an array of states to superimpose, and returns a vector 
            that approximates the (normalized) superimposed state
            states is an array of length at least one
        """
        plt.figure()
        plt.title("Superimposed State vs. Position")
        superimposed = np.zeros(self.dimension)
        for i in states:
            superimposed += self.eigensolution[1][:, i - 1]
        label = ""
        length = len(states)
        for i in range(length):
            if i == length - 1: 
                label += "$\psi_" + str(states[i]) + "$"
            else:
                label += "$\psi_" + str(states[i]) + "$ +"
        plt.plot(self.x_values, superimposed/np.linalg.norm(superimposed), label = label)
        plt.plot(self.x_values, self.potential_values/np.linalg.norm(self.potential_values), label="V(x)")

    def graph_nth_prob_density(self, n):
        """
            takes in a state index and graphs its probability density against position

        Args:
            n (int): Index of the State to be graphed
        """
        plt.figure()
        plt.xlabel("$x$")
        plt.ylabel("$|\psi(x)|^2$")
        title = "Probability Density vs. Position"
        plt.title(title)
        plt.plot(self.x_values, self.eigensolution[1][:, n - 1] * self.eigensolution[1][:, n - 1], label="$|\psi_{"+ str(n) + "}|^2$")
        plt.legend()

    def time_evolved_state(self, n, t, hbar=1):
        """Given a state n and a time t, evaluates the vector approximating the state at time t

        Args:
            n (int): index of state
            t (float): time
        """
        energy_eigenstate = self.eigensolution[1][:, n - 1]
        energy_eigenvalue = self.eigensolution[0][n - 1]
        phase_factor_from_time = np.exp(-1j * energy_eigenvalue / hbar * t)
        return energy_eigenstate * phase_factor_from_time

    def display_graph(self):
        plt.show()


def time_evolution_anim(ase, n):
    """[summary]

    Args:
        ase ([type]): [description]
        n ([type]): [description]
        t ([type]): [description]
    """
    fig, ax = plt.subplots()

    plt.title("Time Evolution")
    def anim_init():
        line.set_data(ase.x_values, ase.eigensolutions[1][n - 1])
        return line,

    def animate(t):
        time_scale = 0.01
        psi_of_t = ase.time_evolved_state(n, t * 0.01)
        line.set_data(ase.x_values, psi_of_t)
        return line
    
    anim = animation.FuncAnimation(fig, animate, init_func=anim_init,
                               frames=200, interval=20, blit=True)
    plt.show()   
#test code
test_equation = SchrodingerEquation(lambda x: x ** 2) # omega=1, m=1, hbar=1

# test_equation.graph_nth_eigenstate(10)
# test_equation.graph_nth_prob_density(10)
# test_equation.graph_superimposed([1, 2, 3])

def wavefunction_time_evolution_anim(ase, n):
    """
        Takes in a Schrodinger Equation datum and an integer,
        and creates an animation demonstrating the time evolution of the 
        n-th excited state of the solution to the Schrodinger Equation

        Args:
            ase (SchrodingerEquation): a SchrodingerEquation datum
            n (int): the index of the excited state
    """
    fig = plt.figure()
    ax = plt.axes(xlim=(-4, 4), ylim=(-0.1, 0.1))
    line, = ax.plot([], [], lw=3)

    def init():
        line.set_data([], [])
        return line,
    def animate(i):
        x = np.linspace(-4, 4, 1000)
        y = ase.time_evolved_state(n, i * 0.01)
        line.set_data(x, y)
        return line,

    anim = FuncAnimation(fig, animate, init_func=init, interval=10, blit=True)
    plt.show()

# print("\nThe discrete apprixmation to V(x):\n", test_equation.discrete_potential_matrix)
# print("\nThe discrete approximation to ddx2:\n", test_equation.ddx2_matrix)


#The first ten eigenvalues are relatively nicely approximated.. after that the linear growth in the approximated eigenvalues gets dominated
# for n in range(1, 40): # the _th eigenvalue
#     print("\nThe " + str(n) + "th energy eigenvalue for this setup is:", len(test_equation.eigensolution[0][n-1]), "\nThe analytic solution for this setup is:", (n - 1) + 0.5)
#test_superimposed = test_equation.superimpose([1, 2, 3])

# plt.figure()
# plt.xlabel("$x$")
# plt.ylabel("$\psi(x)$")
# plt.title("Wavefunction vs. Position")
# # plt.plot(test_equation.x_values, test_equation.eigensolution[1][:, 4], linestyle="--")
# # plt.plot(test_equation.x_values, test_equation.eigensolution[1][:, 6], linestyle="--")
# plt.plot(test_equation.x_values, test_equation.eigensolution[1][:, 10], linestyle="--")
# plt.plot(test_equation.x_values, test_equation.potential_values/np.linalg.norm(test_equation.potential_values), linestyle=":", label="potential")
# #plt.plot(test_equation.x_values, test_superimposed, label="Superimoposed (0, 1, 2) State")
# plt.legend()
# plt.show()