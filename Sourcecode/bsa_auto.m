#############################################################################################################
# Part of the Bayesian Spectrum Analysis (BSA) package in Octave, by Emma
# Granqvist & Richard J. Morris, June 2011
# This code is only proof-of-principle, it is released freely and without any warranty
# For more information see Granqvist et al., BMC Systems Biology 2011,
# and Bretthorst, Springer-Verlag 1988
# Please cite the above papers when using this code
#############################################################################################################
# bsa_auto, wrapper function that stops adding background functions when modelratio > 1 
# INPUT(8):
# data = vecor of data points (as single column)
# name = name (as string)
# start - stop = interval of omega to search initially
# nsamples = number of samples to obtain over the chosen interval
# nbackg = background function maximum expansion order
# interval = time interval between data points (in seconds)
# resultsfile = File to write results into, appended for each time series
# OUTPUT(3):
# vector of normp and omega, and value of omega_max, the omega with max probability
#####################################################################################################
function [normp,omega,omega_max] = bsa_auto(data,name,start,stop,nsamples,nbackg,interval,resultsfile)
  fid = fopen(resultsfile,"a");
  c = ["name ","probability ","omega ","stdevO ","period ","stdevP "];
  fprintf(fid,"%s\n",c);
  fclose(fid);
  [normp,omega,omega_max] = bsa_modelratio_auto(data,name,start,stop,nsamples,nbackg,interval,resultsfile);
endfunction
