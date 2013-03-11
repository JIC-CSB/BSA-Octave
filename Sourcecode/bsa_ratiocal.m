#############################################################################################################
# Part of the Bayesian Spectrum Analysis (BSA) package in Octave, by Emma
# Granqvist & Richard J. Morris, June 2011
# This code is only proof-of-principle, it is released freely and without any warranty
# For more information see Granqvist et al., BMC Systems Biology 2011,
# and Bretthorst, Springer-Verlag 1988
# Please cite the above papers when using this code
#############################################################################################################
# bsa_ratiocal, a part of modelratio calculation
# Calculates terms in model comparion section, to be used in the ratio
########################################################################
function [GL1,GL2,GL31,GL6] = bsa_ratiocal (data,interval,start,stop,nsamples,nbackg,resultsfile)
  omega_range = stop-start;
  ndata = size(data,1);
  tpoints = linspace(interval,ndata*interval,ndata);
  [omega,logp] = bsa_post(data,interval,start,stop,nsamples,nbackg,0); # Run Channel
  [omega,logp,p,noise,power,signaltonoise] = \
      bsa_post(data,interval,start,stop,nsamples,nbackg,max(logp)); # Run again,
				# with max(logp) as upper bound
  [normp] = bsa_normalise(data,interval,omega,p,nsamples,nbackg,max(logp),omega_range); # Normalise p distribution
  logpmax = 0;
    for (j=1:nsamples)
    if logp(j) > logpmax
      logpmax = logp(j);
      omega1 = omega(j);  
    endif 
  endfor
# Model comparion part starts
  r = 1; 
  meanosq = sumsq(omega1);
  fpoints=bsa_samplepoints(tpoints,omega1,nbackg);
  ndata=rows(fpoints);
  nfunc=columns(fpoints);
  [meandsq1,meanhsq1] = bsa_prob_short (transpose(data),fpoints,max(logp));
  GL1 = gamma(nfunc/2)* (((nfunc* meanhsq1)/2)^(-nfunc/2));	
  GL2 = gamma(r/2)*(((r*meanosq)/2)^(-r/2));
  GL31 = (abs(ndata-nfunc-r))/2;
  GL4 = ndata* meandsq1;
  GL5 = nfunc* meanhsq1;
  GL32 = (GL4-GL5)/2;
  fact = (nfunc+r-ndata)/2;
  GL6 = (GL32)^fact;
  if (GL6 == Inf)
    GL6 = 10^300;#preventing numerical error
  endif
endfunction
