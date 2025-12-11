# OptimalResetting
Stochastic resetting using optimal control for return protocols

The python code OptimalResetting.py simulates Brownian motion with Poissonian resetting with a rate r; each resetting event is implemented as an optimal transport protocol towards the origin.
The code is structured in loops, computing the mean work and first passage time to a target l for an ensemble of independent trajectories using various potential stiffness \kappa and protocol duration t_f.
The matlab code plot.m simply plots the restulting data.

Physical parameters are those of the paper, and correspond to a micro-sizerd particle immersed in water in a typicall optical tweezer experiment.

The resutls of this code correspond to the figure shown in the Appendix of the paper, but therefore also contains the main time-energy plot of the main text, for one specific case of stiffness \kappa.
