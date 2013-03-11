#############################################################################################################
# Part of the Bayesian Spectrum Analysis (BSA) package in Octave, by Emma
# Granqvist & Richard J. Morris, June 2011
# This code is only proof-of-principle, it is released freely and without any warranty
# For more information see Granqvist et al., BMC Systems Biology 2011,
# and Bretthorst, Springer-Verlag 1988
# The Nested Sampling routine, following closely Skilling 2004, Sivia & Skilling 2005
# Please cite the above papers when using this code
#############################################################################################################
# bsa_temp, tests mean/variance
###################################
function [m st] = bsa_temp(a)
  m=0.0;
  st=0.0;
  for i=1:length(a)
    m+=a(i);
    st+=a(i)*a(i);
  endfor
  m/=length(a);
  st=sqrt(st/length(a)-m*m);
endfunction