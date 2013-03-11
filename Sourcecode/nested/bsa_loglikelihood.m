#############################################################################################################
# Part of the Bayesian Spectrum Analysis (BSA) package in Octave, by Emma
# Granqvist & Richard J. Morris, June 2011
# This code is only proof-of-principle, it is released freely and without any warranty
# For more information see Granqvist et al., BMC Systems Biology 2011,
# and Bretthorst, Springer-Verlag 1988
# The Nested Sampling routine, following closely Skilling 2004, Sivia & Skilling 2005
# Please cite the above papers when using this code
#############################################################################################################
# log likelihood function - needs changing for different problems
#############################################################################
function ll = bsa_loglikelihood (data,omega,interval,nbackg,dim)
  [ll] = bsa_post_nested(data,omega,interval,nbackg,dim,0); 
endfunction