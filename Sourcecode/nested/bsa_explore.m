#############################################################################################################
# Part of the Bayesian Spectrum Analysis (BSA) package in Octave, by Emma
# Granqvist & Richard J. Morris, June 2011
# This code is only proof-of-principle, it is released freely and without any warranty
# For more information see Granqvist et al., BMC Systems Biology 2011,
# and Bretthorst, Springer-Verlag 1988
# The Nested Sampling routine, following closely Skilling 2004, Sivia & Skilling 2005
# Please cite the above papers when using this code
#############################################################################################################
# Develop a new object from a given object x above a given loglikelihood
##################################################################
function [x xLL] = bsa_explore (data,x, omega_min, omega_max, logLhoodMin,dim,interval,nbackg)
  step=0.1;
  m=20;
  accept=0;
  reject=0;
  xLL=bsa_loglikelihood (data,x,interval,nbackg,dim);
  for i=1:m
#    tryx=x+step*(2.0*rand()-1.0);
    tryx=normrnd(x,step);
# this ensures that the new value is within the specified limits
    tryx=omega_min + mod((tryx-omega_min),(omega_max-omega_min));
    #tryx-=floor(tryx); # this needs adjusting to make sure that the
		       # trial values is within the specified range
    #disp(newx);
    tryxLL=bsa_loglikelihood (data,tryx,interval,nbackg,dim);
    if (tryxLL > logLhoodMin)
      x=tryx;
      xLL=tryxLL;
      #disp(x);
      #disp(newxLL);
      accept++;
    else
      reject++;
    endif
    if (accept>reject)
      step *= exp(1.0/accept);
    endif
    if (accept<reject)
      step /= exp(1.0/reject);
    endif
    # d isp(step);
  endfor
#  if (xLL <= logLhoodMin)
#    disp ("explore failed");
#  endif
endfunction