#############################################################################################################
# Part of the Bayesian Spectrum Analysis (BSA) package in Octave, by Emma
# Granqvist & Richard J. Morris, June 2011
# This code is only proof-of-principle, it is released freely and without any warranty
# For more information see Granqvist et al., BMC Systems Biology 2011,
# and Bretthorst, Springer-Verlag 1988
# Please cite the above papers when using this code
#############################################################################################################
# bsa_modelratio_auto, part of model development,
# calculating model ratio and stops adding background functions when ratio > 1
####################################################################################################
function [normp,omega,omega_max] = bsa_modelratio_auto(data,name,start,stop,nsamples,nbackg,interval,fid)
  for (i=0:nbackg)
    if (i==nbackg)
      [normp,omega,omega_max] = bsa_post_norm(data,name,start,stop,nsamples,i,interval,fid);
    else
      [GL1_1,GL2_1,GL31_1,GL4_1] = bsa_ratiocal (data,interval,start,stop,nsamples,i,fid);
      [GL1,GL2,GL31,GL4] = bsa_ratiocal (data,interval,start,stop,nsamples,i+1,fid);
      if (GL31_1 > GL31)
	n = GL31;
	GL3 = gamma(0.5)/beta(n,0.5);
      else
	n = GL31_1;
	GL3 = beta(n,0.5)/gamma(0.5);
      endif
      ratio = 0;
      prob = gamma(GL31_1);
      like1 = GL1_1*GL2_1*prob*GL4_1;
      ratio = ((GL1_1*GL2_1*GL4_1)/(GL1*GL2*GL4))*GL3;
      disp(['Ratio, ',num2str(i),' and ',num2str(i+1),' : ',num2str(ratio)]);
      if (ratio > 1) 
	[normp,omega,omega_max] = bsa_post_norm(data,name,start,stop,nsamples,i,interval,fid);
	break;
      else
	GL1_1 = GL1;
	GL2_1 = GL2;
	GL31_1 = GL31;
	GL4_1 = GL4;
      endif
    endif
  endfor
endfunction


