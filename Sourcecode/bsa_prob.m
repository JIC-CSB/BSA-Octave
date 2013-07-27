#############################################################################################################
# Part of the Bayesian Spectrum Analysis (BSA) package in Octave, by Emma
# Granqvist & Richard J. Morris, June 2011
# This code is only proof-of-principle, it is released freely and without any warranty
# For more information see Granqvist et al., BMC Systems Biology 2011,
# and Bretthorst, Springer-Verlag 1988
# Please cite the above papers when using this code
#############################################################################################################
# bsa_prob, calculates the "Student T-distribution" Eq. 3.17 in Bretthorst, 1988
# INPUT(3):
# data = time series (transposed)
# fvalues = sampled points from function (samplef)
# maxlogST = max value of logST, if computed earlier, otherwise 0

# ndata = number of samples from samplef ('N')
# nfunc = number of extensions Legendre ('m')
# orthofvalues = orthonormalized matrix of sampled fvalues
# h = projection of the data onto orthonormal model functions, Eq. 3.13 & 4.3. in Bretthorst, 1988
# meanhsq = Mean-squared of h, Eq. 3.15 in Bretthorst, 1988
# meandsq = Mean-squared of data value observed
# factor = bracketed part of Student t-dist 

# OUTPUT(5):
# logST = log(Student t-dist)
# ST = Student t-dist
# sigma = estimated noise, Eq. 4.6, 4.7 in Bretthorst, 1988
# spden = power spectral density, Eq. 4.15 in Bretthorst, 1988
# signaltonoise = signaltonoiseratio Eq. 4.8 in Bretthorst, 1988
#####################################################################################################
function [logST,ST,sigma,spden,signaltonoise] = bsa_prob (data,fvalues,maxlogST)
  ndata=rows(fvalues);
  nfunc=columns(fvalues);
  [orthofvalues,evalues] = bsa_orthonormalize(fvalues);
  h = data*orthofvalues;
  meanhsq = sumsq(h)/nfunc;
  meandsq = sumsq(data)/ndata;
  factor = 1.0 - nfunc*meanhsq/ndata/meandsq;
    if (abs(factor) < 1.0e-14)
      factor = 1.0e-14; #preventing numerical error
    endif 
  logST = log (factor) * (nfunc-ndata)/2.0;
  logdiff = logST - maxlogST;
  ST = 0.0;
  if (maxlogST != 0.0) 
    ST=exp(logdiff);
  endif 
  sigma = sqrt(ndata/(ndata-nfunc-2) * (meandsq - nfunc*meanhsq/ndata));
  spden = nfunc * meanhsq * ST;
  signaltonoise = sqrt(nfunc/ndata*(1+meanhsq/sigma/sigma));
endfunction
