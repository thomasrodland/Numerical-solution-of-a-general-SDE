# Numerical-solution-of-a-general-SDE
Euler-Maruyama and Milstein methods for solving a general SDE of the form $d X(t)=a(t, X(t)) d t+b(t, X(t)) d B(t), 0<t<T$
Normal_generator.m provides the random number needed
Euler_Maruyama.m implements the Euler - Maruyama method 
Milstein.m implements the Milstein method
EU_M_and_milstein implements a solution to $d S(t)=\mu S(t) d t+\sigma S(t) d B(t), 0<t<T$ with previously define numerical methoods and looks at the string and weak convergence of both methods with $500000$ simulations with different values of h

