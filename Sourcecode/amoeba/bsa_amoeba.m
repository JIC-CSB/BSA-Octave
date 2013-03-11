#############################################################################################################
# Part of the Bayesian Spectrum Analysis (BSA) package in Octave, by Emma
# Granqvist & Richard J. Morris, June 2011
# This code is only proof-of-principle, it is released freely and without any warranty
# For more information see Granqvist et al., BMC Systems Biology 2011,
# and Bretthorst, Springer-Verlag 1988
# Please cite the above papers when using this code
#############################################################################################################
# bsa_amoeba, 1D Downhill-Simplex optimisation. P and R are the two
# moving points. xini is the initial guess, lambda is the initial
# distance between the two points, limit is the covergence criteria. The
# output is the optimal logp value and its associated omega
#########################################################
function [logp,omega1] = bsa_amoeba (x,xini,lambda,tpoints,nbackg,data)
  P = xini;
  R = xini + lambda;
  limit = 0.000001;
  [P,R] = bsa_amoeba_eval(P,R,tpoints,nbackg,data);# Identify the two key points
  for (i = 1:x) # Main loop
    [P,R] = bsa_amoeba_eval(P,R,tpoints,nbackg,data);# Identify the two key points
    if ((bsa_amoeba_function(R,tpoints,nbackg,data)(1) - bsa_amoeba_function(P,tpoints,nbackg,data)(1)) < limit) # Convergence?
      break;
    else
      [P,R] = bsa_amoeba_reflect(P,R); # Reflect P
      if (bsa_amoeba_function((P)(1),tpoints,nbackg,data) > \
	  bsa_amoeba_function(R,tpoints,nbackg,data)(1)) # P better?
	[P] = bsa_amoeba_extrapol(P,R); # Good move, do more
      else # P still the worst?
	[P,R] = bsa_amoeba_contract(P,R); # contraction away from P
      endif
    endif
  endfor
  [P,R] = bsa_amoeba_eval(P,R,tpoints,nbackg,data);# Final ID of the two key points
  [logp,omega1] = bsa_amoeba_function(R,tpoints,nbackg,data);
  if ((bsa_amoeba_function(R,tpoints,nbackg,data)(1) - bsa_amoeba_function(P,tpoints,nbackg,data)(1)) < limit) # Convergence?
  #  disp("converged!");
  else
  #  disp("not converged, but the looping finished");
  endif
endfunction