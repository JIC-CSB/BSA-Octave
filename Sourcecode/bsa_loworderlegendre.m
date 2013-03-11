#############################################################################################################
# Part of the Bayesian Spectrum Analysis (BSA) package in Octave, by Emma
# Granqvist & Richard J. Morris, June 2011
# This code is only proof-of-principle, it is released freely and without any warranty
# For more information see Granqvist et al., BMC Systems Biology 2011,
# and Bretthorst, Springer-Verlag 1988
# Please cite the above papers when using this code
#############################################################################################################
# bsa_loworderlegendre, calculating the legendre polynomials as
# background functions, works up to order 10
# Input(2)
# n = number of extentions
# x = vector of timepoints, scaled
# Output(1):
# leg = function out
#####################################################################################################
function leg = bsa_loworderlegendre(n,x)
  leg = -1;
  if (n==0) 
    leg = 1;
  elseif (n==1) 
    leg = x;
  elseif (n==2) 
    leg = 0.5*(3*x^2-1);
  elseif (n==3) 
    leg = 0.5*(5*x^3-3*x);
  elseif (n==4) 
    leg = 1/8*(35*x^4-30*x^2+3);
  elseif (n==5) 
    leg = 1/8*(63*x^5-70*x^3+15*x);
  elseif (n==6)
    leg = 1/16*(231*x^6-315*x^4+105*x^2-5);
  elseif (n==7)
    leg = (1/16)*(429*x^7-693*x^5+315*x^3-35*x);
  elseif (n==8)
    leg = (1/128)*(6435*x^8-12012*x^6+6930*x^4-1260*x^2+35);
  elseif (n==9)
    leg = (1/128)*(12155*x^9-25740*x^7+18018*x^5-4620*x^3+315*x);
  elseif (n==10)
    leg = (1/256)*(46189*x^10-109395*x^8+90090*x^6-30030*x^4+3465*x^2-63);
  else
    printf ("nbackg only for n<=10! \n");
  endif
endfunction 