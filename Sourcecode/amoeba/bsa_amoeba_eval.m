#############################################################################################################
# Part of the Bayesian Spectrum Analysis (BSA) package in Octave, by Emma
# Granqvist & Richard J. Morris, June 2011
# This code is only proof-of-principle, it is released freely and without any warranty
# For more information see Granqvist et al., BMC Systems Biology 2011,
# and Bretthorst, Springer-Verlag 1988
# Please cite the above papers when using this code
#############################################################################################################
# bsa_amoeba_eval, ranks the two points (P and R) for fit, for Downhill-Simplex optimisation 
#####################################################################
function [Pnew,Rnew] = bsa_amoeba_eval(P,R,tpoints,nbackg,data)
  if ((bsa_amoeba_function(P,tpoints,nbackg,data)(1)) >= (bsa_amoeba_function(R,tpoints,nbackg,data)(1)))
    Rnew = P;
    Pnew = R;
    else
      Pnew = P;
      Rnew = R;
  endif
endfunction