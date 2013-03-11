#############################################################################################################
# Part of the Bayesian Spectrum Analysis (BSA) package in Octave, by Emma
# Granqvist & Richard J. Morris, June 2011
# This code is only proof-of-principle, it is released freely and without any warranty
# For more information see Granqvist et al., BMC Systems Biology 2011,
# and Bretthorst, Springer-Verlag 1988
# Please cite the above papers when using this code
#############################################################################################################
# bsa_prob_short, part of the bsa_ratiocal function
# shorter version of bsa_prob
################################################
function [meandsq,meanhsq] = bsa_prob_short (data,fvalues,maxlogST)
  ndata=rows(fvalues);
  nfunc=columns(fvalues);
  [orthofvalues,evalues] = bsa_orthonormalize (fvalues);
  h = data*orthofvalues;
  meanhsq = sumsq(h)/nfunc;
  meandsq = sumsq(data)/ndata;
endfunction
