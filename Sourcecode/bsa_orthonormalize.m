#############################################################################################################
# Part of the Bayesian Spectrum Analysis (BSA) package in Octave, by Emma
# Granqvist & Richard J. Morris, June 2011
# This code is only proof-of-principle, it is released freely and without any warranty
rth
# For more information see Granqvist et al., BMC Systems Biology 2011,
# and Bretthorst, Springer-Verlag 1988
# Please cite the above papers when using this code
#############################################################################################################
# bsa_orthonormalize, orthonormalizes the model functions. A part of bsa_prob.
# a = matrix
#####################################################################################################
function [ortha,evalues] = bsa_orthonormalize(a)
  [evectors,evalues] = eig(a'*a);
  ortha = a*evectors;
  anorm = sqrt(sumsq(ortha));
  for j = 1:columns(a)
      ortha(:,j) = ortha(:,j)/anorm(j);
  endfor
endfunction
