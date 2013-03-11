#############################################################################################################
# Part of the Bayesian Spectrum Analysis (BSA) package in Octave, by Emma
# Granqvist & Richard J. Morris, June 2011
# This code is only proof-of-principle, it is released freely and without any warranty
# For more information see Granqvist et al., BMC Systems Biology 2011,
# and Bretthorst, Springer-Verlag 1988
# Please cite the above papers when using this code
#############################################################################################################
# bsa_post_norm, calculates the posterior probability density,
# normalised, and prints results to 'resultsfile'
# Input(9):data (1) with some name (2), sampling from start (3) to stop (4), number of samples nsamples (5), nbackg (6) = Legendre polynomials expansion order, interval (8) time interval between sample points, resultsfile (9) name of file where string of results is printed
# Output(3):vector of normalised posterior ditribution, vector of omega, 
# and value of omega value with highest probability (omega_max)
##################################################################################################
function [normp2,omega2,omega_max] = bsa_post_norm (data,name,start,stop,nsamples,nbackg,interval,resultsfile)
  omega_range = stop-start;
  omegaN = omega_range*0.05;
  tpoints = linspace(interval,size(data,1)*interval,size(data,1));
  [omega,logp] = bsa_post(data,interval,start,stop,nsamples,nbackg,0); # Run Channel
  [omega,logp,p,noise,power,signaltonoise] = bsa_post(data,interval,start,stop,nsamples,nbackg,max(logp)); 
  maxp = 0;
  for (j=1:nsamples)
    if (p(j) > maxp)
      maxp = p(j);
      guess = omega(j);
  endif 
  endfor
  [normp] = bsa_normalise_output(data,name,interval,omega,p,nsamples,nbackg,max(logp),resultsfile,omega_range); # Normalise 
  [maxlogp,omega1]=bsa_amoeba(1000,guess,omegaN*2,tpoints,nbackg,data);# Search for maximum with amoeba
  [omega2,logp,p,noise,power,signaltonoise]= bsa_post(data,interval,(omega1-omegaN),(omega1+omegaN),nsamples,nbackg,maxlogp);# Channel again, local
  [normp2,omega_max] = bsa_normalise_output(data,name,interval,omega2,p,nsamples,nbackg,max(logp),resultsfile,omega_range); # Normalise 
endfunction
