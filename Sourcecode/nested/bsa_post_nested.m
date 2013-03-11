#############################################################################################################
# Part of the Bayesian Spectrum Analysis (BSA) package in Octave, by Emma
# Granqvist & Richard J. Morris, June 2011
# This code is only proof-of-principle, it is released freely and without any warranty
# For more information see Granqvist et al., BMC Systems Biology 2011,
# and Bretthorst, Springer-Verlag 1988
# The Nested Sampling routine, following closely Skilling 2004, Sivia & Skilling 2005
# Please cite the above papers when using this code
#############################################################################################################
# bsa_post_nested for nested sampling routine
# Input(5): data, nsamples (number of samples),nbackg (Legendre polynomials expansion order), dim (dimensions, in this case = 1), normp (in this case = 0)
# Output(1):
# Value of logp for point omega
#####################################################################################################
function [logp] = bsa_post_nested (data,omega,interval,nbackg,dim,normp) 
  tpoints = linspace(interval,size(data,1)*interval,size(data,1));
  fpoints= bsa_samplepoints(tpoints,omega,nbackg);
  [a,b,c,d,f] = bsa_prob (transpose(data),fpoints,normp);
  logp=a;
endfunction