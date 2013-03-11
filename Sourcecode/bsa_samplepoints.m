#############################################################################################################
# Part of the Bayesian Spectrum Analysis (BSA) package in Octave, by Emma
# Granqvist & Richard J. Morris, June 2011
# This code is only proof-of-principle, it is released freely and without any warranty
# For more information see Granqvist et al., BMC Systems Biology 2011,
# and Bretthorst, Springer-Verlag 1988
# Please cite the above papers when using this code
#############################################################################################################
# bsa_samplepoints, samples points from model functions
# Input(3): 
# tpoints= vector of timepoints
# omega= single value 
# nbackg= no. of background functions, here Legendre polynomials used
# Output(1):
# f= vector of function values sampled
#####################################################################################################
function f = bsa_samplepoints(tpoints,omega,nbackg)
  f = [];
  tscale = tpoints((columns(tpoints)))-tpoints(1);
  for (i = 1:columns(tpoints))
    f(i,1) = 1;
    for (n = 1:nbackg)
      f(i,n+1) =  bsa_loworderlegendre(n, tpoints(i)/tscale-0.5);
    endfor
    f(i,nbackg+2) = sin (tpoints(i)*omega);
    f(i,nbackg+3) = cos (tpoints(i)*omega);
  endfor
endfunction