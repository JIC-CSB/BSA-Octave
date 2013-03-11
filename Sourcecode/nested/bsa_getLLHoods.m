#############################################################################################################
# Part of the Bayesian Spectrum Analysis (BSA) package in Octave, by Emma
# Granqvist & Richard J. Morris, June 2011
# This code is only proof-of-principle, it is released freely and without any warranty
# For more information see Granqvist et al., BMC Systems Biology 2011,
# and Bretthorst, Springer-Verlag 1988
# The Nested Sampling routine, following closely Skilling 2004, Sivia & Skilling 2005
# Please cite the above papers when using this code
#############################################################################################################
# Computes the log likelihood values for a set of objects
#################################################################################
function llvalues = bsa_getLLHoods (data,x,interval,nbackg,dim)
  for i=1:length(x)
    llvalues(i)=bsa_loglikelihood(data,x(i,:),interval,nbackg,dim);
  endfor
endfunction