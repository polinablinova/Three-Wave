# Three-Wave Resonant Interaction Solver
Solves equations modelling laser pulse compression via stimulated Brillouin scattering (SBS) in the strongly and finitely damped regimes, as in https://arxiv.org/abs/2505.17956
SBS is implemented in a counter-propagating geometry, where a pump laser (a) and a seed laser (b) interact resonantly via the excited acoustic field (f). Energy is transferred from the pump to a shorter counter-propagating seed. 
## jihoon_no_extra_term_in_acou.m
Implements a function based on RK4 to solve the strongly damped SBS (see Eq. 3 in the paper referenced above). Simulation is done in the frame of the seed, propagating from right to left, and interacting with a pump that is being continuously fed from the left edge of the simulations domain.
## run_two_wave.m
Runs jihoon_no_extra_term_in_acou.m
## ringing_back.m
Implements a function to solve the finitely damped SBS (Eqs. 1,2) based on forward Euler method for the pump and seed lasers and backward Euler and RK4 for the acoustic field. 
## run_ringing_back.m
Runs ringing_back.m
## sbs_params.nb
Converts physical parameters into simulation inputs for jihoon_no_extra_term_in_acou.m and ringing_back.m
