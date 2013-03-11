#############################################################################################################
# Part of the Bayesian Spectrum Analysis (BSA) package in Octave, by Emma
# Granqvist & Richard J. Morris, June 2011
# This code is only proof-of-principle, it is released freely and without any warranty
# For more information see Granqvist et al., BMC Systems Biology 2011,
# and Bretthorst, Springer-Verlag 1988
# The Nested Sampling routine, following closely Skilling 2004, Sivia & Skilling 2005
# Please cite the above papers when using this code
#############################################################################################################
# compute average results from posterior samples
############################################################
function [momega stomega mperiod stperiod] = bsa_meanomega (posts,postsLWeights,logZ)
  momega=0.0;
  stomega=0.0;
  mperiod=0.0;
  stperiod=0.0;
  period = [];
## Into period, in s
  for (j=1:length(posts))
    period(j) = (2*pi)/posts(j);
  endfor
## For omega
  for i=1:length(posts)
    weight=exp(postsLWeights(i)-logZ);
    momega+= weight*posts(i);
    stomega+=weight*posts(i)*posts(i);
  endfor
  stomega=sqrt(stomega-momega*momega);
## For periods
  for k=1:length(posts)
    weight=exp(postsLWeights(k)-logZ);
    mperiod+= weight*period(k);
    stperiod+=weight*period(k)*period(k);
  endfor
  stperiod=sqrt(stperiod-mperiod*mperiod);
endfunction