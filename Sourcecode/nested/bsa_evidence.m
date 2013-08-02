#############################################################################################################
# Part of the Bayesian Spectrum Analysis (BSA) package in Octave, by Emma
# Granqvist & Richard J. Morris, June 2011
# This code is only proof-of-principle, it is released freely and without any warranty
# For more information see Granqvist et al., BMC Systems Biology 2011,
# and Bretthorst, Springer-Verlag 1988
# The Nested Sampling routine, following closely Skilling 2004, Sivia & Skilling 2005
# Please cite the above papers when using this code
#############################################################################################################
# bsa_evidence, computes the evidence and average values of omega
# Input(7):data, omega_min and omega_max are upper and lower bound for
# omega sampling, dim = dimension of search, nsamples = number of
# samples over omega interval, nbackg = number of background functions,
# interval = time (s) between sample points in the data, nposts = number
# of sample points from posterior for calculating evidence
# Output(6): logZ = evidence, logZerror, evidence error, momega = mean
# omega value, stomega = standard deviation of omega, mperiod =  mean
# period, stperiod = standard deviation of omega
##################################################
function [logZ logZerror momega stomega mperiod stperiod posts] = bsa_evidence (data,omega_min,omega_max,dim,nsamples,nbackg,interval,nposts)
  samples=bsa_prior(omega_min, omega_max, dim, nsamples);
  samplesLLHoods = bsa_getLLHoods(data,samples,interval,nbackg,dim);
  samplesLWeights = zeros(nsamples,1);
  posts = zeros(nposts,1);
  postsLWeights = zeros(nposts,1);
  H=0.0;
  logZ=-realmax();
  logwidth=log(1.0-exp(-1.0/nsamples));
  for i=1:nposts 
    worst=1;
    for j=2:nsamples
      if (samplesLLHoods(j) < samplesLLHoods(worst))
	worst=j;
      endif
    endfor
    samplesLWeights(worst)=logwidth+samplesLLHoods(worst);
    logZnew = bsa_logplus(logZ,samplesLWeights(worst));
    H=exp(samplesLWeights(worst)-logZnew)*samplesLLHoods(worst)+exp(logZ-logZnew)*(H+logZ)-logZnew;
    logZ=logZnew;
# save these samples for later analysis (average parameters, etc.)
    posts(i)=samples(worst);
    postsLWeights(i)=samplesLWeights(worst);
    do 
      newi=1+floor(nsamples*rand());
      #newi=randint(1,1,[1,nsamples]);
    until (newi != worst);
    LLHoodmin=samplesLLHoods(worst);
    samples(worst)=samples(newi);
    [samples(worst) \
     samplesLLHoods(worst)]=bsa_explore(data,samples(worst),omega_min, \
				    omega_max, LLHoodmin, dim,interval,nbackg);
    #disp(samplesLLHoods(worst));
    logwidth-=1.0/nsamples;
  endfor
  logZerror=sqrt(H/nsamples);
  [momega stomega mperiod stperiod] = bsa_meanomega(posts,postsLWeights,logZ);
endfunction
