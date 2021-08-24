# quantum-toolkit
 
This is a program that (currently) numerically approximates solutions to the Schrodinger Equation for any given potential function in one dimension. The main functionalities are provided by the `numpy` and `matplotlib` libraries. The core of the approximation relies on a crude numerical estimation of the Hessian matrix, and such an approximation will be relatively well-behaved for around the first 10 energy levels or so (tested with the harmonic oscillator potential). Look at `test.py` for some examples of code usage. 

The applcation can be started by running the command `python schrodinger_solve.py`, or for the interactive application, `python -i schrodinger_solve.py`.

All usage should start by creating a `SchroedingerEquation`, whose constructor takes in a function representing the potential of the Hamiltonian. 

Currently, a `SchroedingerEquation` object supports the following method calls:
- `graph_nth_eigenstate`: takes in an integer `n` and graphs the approximated `n`-th eigenstate for the given `SchroedingerEquation`
- `graph_superimposed`: takes in a list of integers and graphs the superimposed states of the states given
- `graph_nth_probability`: takes in an integer `n` and graphs the probability density function of the `n`-th eigenstate
- `time_evolved_state`: takes in an integer `n` and a float `t` and graphs the `n`-th eigenstate at time `t`

In addition, calling `time_evolve_anim` on a `SchroedingerEquation` and an integer `n` creates an animation of the time evolution of the `n`-th eigenstate. 

Stay tuned for future improvements and additions to functionalities!
