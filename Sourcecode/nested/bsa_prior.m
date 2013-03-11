#############################################################################################################
# Part of the Bayesian Spectrum Analysis (BSA) package in Octave, by Emma
# Granqvist & Richard J. Morris, June 2011
# This code is only proof-of-principle, it is released freely and without any warranty
# For more information see Granqvist et al., BMC Systems Biology 2011,
# and Bretthorst, Springer-Verlag 1988
# The Nested Sampling routine, following closely Skilling 2004, Sivia & Skilling 2005
# Please cite the above papers when using this code
#############################################################################################################
# bsa_prior, sets up a uniform sample of n omegas between min and max
# Can be easily adapted for non-uniform prior
################################################
function x = bsa_prior (omega_min, omega_max, dim, n)
  x = omega_min + (omega_max-omega_min)*rand(n,dim); # uniform prior
endfunction