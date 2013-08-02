#############################################################################################################
# Part of the Bayesian Spectrum Analysis (BSA) package in Octave, by Emma
# Granqvist & Richard J. Morris, June 2011
# This code is only proof-of-principle, it is released freely and without any warranty
# For more information see Granqvist et al., BMC Systems Biology 2011,
# and Bretthorst, Springer-Verlag 1988
# Please cite the above papers when using this code
#############################################################################################################
# bsa_nest, wrapper function to compute evidence and mean value of omega
# INPUT(9):
# data = vecor of data points (as single column)
# name = name (as string)
# start - stop = interval of omega to search initially
# nsamples = number of samples to obtain over the chosen interval
# nbackg = background function maximum expansion order
# interval = time interval between data points (in seconds)
# resultsfile = File to write results into, appended for each time series
# nposts: number of samples drawn from posterior to compute evidence
# OUTPUT(4): logZ = evidence, logZerror, evidence error, momega = mean
# omega value, stomega = standard deviation of omega
#####################################################################################################
function [logZ logZerror momega stomega posts] = bsa_nest(data,name,start,stop,nsamples,nbackg,interval,resultsfile,nposts)
  fid = fopen(resultsfile,"a");
  c = ["name ","evidence ","evidence_error ","mean_omega ","st_omega ","mean_period ","st_period "];
  fprintf(fid,"%s\n",c);
  fclose(fid);
  [logZ logZerror momega stomega mperiod stperiod posts] = bsa_evidence(data,start,stop,1,nsamples,nbackg,interval,nposts);
  fid = fopen(resultsfile,"a");
  res = [name," " ,num2str(logZ)," ",num2str(logZerror)," " ,num2str(momega)," ",num2str(stomega)," ",num2str(mperiod)," ",num2str(stperiod)];
  fprintf(fid,"%s\n",res);
  fclose(fid);
endfunction
