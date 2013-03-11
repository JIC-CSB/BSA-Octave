#############################################################################################################
# Part of the Bayesian Spectrum Analysis (BSA) package in Octave, by Emma
# Granqvist & Richard J. Morris, June 2011
# This code is only proof-of-principle, it is released freely and without any warranty
# For more information see Granqvist et al., BMC Systems Biology 2011,
# and Bretthorst, Springer-Verlag 1988
# Please cite the above papers when using this code
#############################################################################################################
# bsa_normalise_output, normalise probability density function, and
# writes results to resultsfile
#####################################################################################################
function [norm_p,omega_max] = bsa_normalise_output(data,name,interval,omega,p,nsamples,nbackg,logp_max,resultsfile,omega_range)
  omegaN = omega_range/1000;
  dx = 0;
  meanP = 0;
  meanO = 0;
  varP = 0;
  varO = 0;
  stO = 0;
  stP = 0;
  pmax = 0;
  omega_max = 0;
  period_max = 0;
#Start
  for (j=1:size(p)(2))
    period(j) = (2*pi)/omega(j);
    if p(j) > pmax
      pmax = p(j);
      omega_max = omega(j);
      period_max = period(j);
    endif 
  endfor
  for (i=1:size(p)(2))
    if ((p(i) >  0.001) && ((omega_max - omega(i)) > omegaN))
      omegaN = omega_max - omega(i);
    endif 
  endfor
  [omega_new,logp_new,p_new] = bsa_post (data,interval,omega_max-omegaN,omega_max+omegaN,nsamples,nbackg,0);
  [omega_new,logp_new,p_new] = bsa_post (data,interval,omega_max-omegaN,omega_max+omegaN,nsamples,nbackg,max(logp_new));

  dx = (omegaN*2)/nsamples;

  norm_p = [];
  sum1 = sum(p_new*dx);
  for (i=1:size(p)(2))
    norm_p(i) = (p(i)/sum1);
    normp_new(i) = (p_new(i)/sum1);
  endfor
  for (k=1:size(p_new)(2))
    period(k) = 1/(2*pi*omega_new(k));
    varO = varO + dx*normp_new(k)*omega_new(k)^2;
    meanO = meanO + dx*normp_new(k)*omega_new(k);
    varP = varP + dx*normp_new(k)* period(k)^2; 
    meanP = meanP + dx*normp_new(k)*period(k);
  endfor
  varO = varO - (meanO)^2;
  stP = sqrt(varP);
  stO = sqrt(varO);
# Write Results
  results = [max(norm_p),omega_max,stO,period_max,stP];
  results = num2str(results);
  results = [name," ",results];
  fid = fopen(resultsfile,"a");
  fprintf(fid,"%s\n",results);
  fclose(fid); 
endfunction