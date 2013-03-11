#############################################################################################################
# Part of the Bayesian Spectrum Analysis (BSA) package in Octave, by Emma
# Granqvist & Richard J. Morris, June 2011
# This code is only proof-of-principle, it is released freely and without any warranty
# For more information see Granqvist et al., BMC Systems Biology 2011,
# and Bretthorst, Springer-Verlag 1988
# Please cite the above papers when using this code
#############################################################################################################
# bsa_amoeba_function, function to maximise, for Downhill-Simplex optimisation 
#########################################################
function [logp,omega] = bsa_amoeba_function(omega,tpoints,nbackg,data)
    fpoints = [];
    fpoints=  bsa_samplepoints(tpoints,omega,nbackg);
    [logp] = bsa_prob (transpose(data),fpoints,0);
endfunction
