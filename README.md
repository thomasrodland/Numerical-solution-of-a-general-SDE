# Numerical-solution-of-a-general-SDE
This repo contains Euler-Maruyama and Milstein methods for solving a general SDE of the form $$d X(t)=a(t, X(t)) d t+b(t, X(t)) d B(t)$$

The repo contains the following files:

 - Normal_generator.m, provides the random number needed.

 - Euler_Maruyama.m, implements the Euler-Maruyama method.

 - Milstein.m, implements the Milstein method.

 - EU_M_and_milstein.m, implements a solution to $d S(t)=\mu S(t) d t+\sigma S(t) d B(t)$ with previously defined numerical methoods and looks at the strong and weak convergence of both methods with $500000$ simulations with different values of h.

 - plot_solution_sde.m, plots a solution using the specified numerical method.
