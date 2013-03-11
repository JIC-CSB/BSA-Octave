#############################################################################################################
# Part of the Bayesian Spectrum Analysis (BSA) package in Octave, by Emma
# Granqvist & Richard J. Morris, June 2011
# This code is only proof-of-principle, it is released freely and without any warranty
# For more information see Granqvist et al., BMC Systems Biology 2011,
# and Bretthorst, Springer-Verlag 1988
# Please cite the above papers when using this code
#############################################################################################################
# bsa_post, calculates the posterior, before normalisation
# Input(7): data (1) with some interval (2), sampling from start (3) to stop (4), number of samples nsamples (5), nbackg (6) = Legendre polynomials expansion order, normalised logp value normp (7), needed for p, noise and power results to be calculated
# Output(6): vectors of omega, logp, p, noise, power spectrum and signal-noise ratio
# logpmax = peak value, frequency = corresponding key frequency
#####################################################################################################
function [omega,logp,p,noise,power,signaltonoise] = bsa_post (data,interval,start,stop,nsamples,nbackg,normp)
  omega = [];
  logp = [];
  p = [];
  noise = [];
  power = [];
  logpmax = 0;
  omega1 = 0;
  frequency = 0;
  tpoints = linspace(interval,size(data,1)*interval,size(data,1));
  for (i=1:nsamples)
    omega(i)=start+i*(stop-start)/nsamples;
    fpoints=bsa_samplepoints(tpoints,omega(i),nbackg);
    [a,b,c,d,f] = bsa_prob (transpose(data),fpoints,normp);
    logp(i)=a;
    p(i)=b;
    noise(i)=c;
    power(i)=d;
    signaltonoise(i)=f;
  endfor
endfunction
