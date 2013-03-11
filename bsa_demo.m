#############################################################################################################
# Part of the Bayesian Spectrum Analysis (BSA) package in Octave, by Emma
# Granqvist & Richard J. Morris, June 2011
# This code is only proof-of-principle, it is released freely and without any warranty
# For more information see Granqvist et al., BMC Systems Biology 2011,
# and Bretthorst, Springer-Verlag 1988
# Please cite the above papers when using this code
#############################################################################################################
# DEMO script
######################################################################################################
disp(['Welcome to the BSA demo!']);
disp(['Here we will go through four basic examples.']);
disp(['Press enter to start with the first one!']); pause
disp(' ');
##############################  Example 1
disp(['Case 1: A time series with omega = 0.5']);
dpoints = [];
tpoints = linspace(1,200,200);
for (i=1:200)
  dpoints(i) = sin(0.5*i);
endfor
figure; # Plot time series
plot(tpoints,dpoints);
title("Example time series","fontsize",16);
xlabel("t","fontsize",14);
ylabel("d(t)","fontsize",14);
set(gca(),"xlim",[1, 200],"ylim",[-1.2, 1.2]);
disp(['To start BSA, we enter the following on the command line:']); 
disp(['[normp,omega,maxomega] = bsa_auto (dpoints,"Test1",0,1,100,0,1,"BSA-DEMO");']); 
disp(['where dpoints is the data, as a single column,']); 
disp(['Test1 is the name of the time series,']);
disp(['0 is the lower limit of omega and 1 is the upper,']);
disp(['100 is the number of sample points between the lower and upper omega,']);
disp(['0 is the upper limit of background functions,']);
disp(['1 is the interval between sampling points in the data (in seconds),']);
disp(['and "BSA-DEMO" is the name of the output file.']);
disp(['Press enter to start BSA']); pause
fprintf (stderr, "Working... This may take a few minutes.\n"); 
[normp,omega,maxomega] = bsa_auto (dpoints',"Test1",0,1,100,0,1,"BSA-DEMO");
disp(['Done. The output is the estimated omega (maxomega = ' ,num2str(maxomega),'), and vectors of omega and the posterior probability.']); 
disp(['Press enter to plot omega and posterior density']); pause
plot(omega,normp,"1","linewidth",1);# Plot results
title("BSA result","fontsize",16);
xlabel('\omega',"fontsize",14);
ylabel("PDF","fontsize",14);
set(gca(),"xlim",[0, 1]);
disp(['Press enter to go on to the next example']); pause
##############################  Example 2
disp(['Case 2: A time series with high noise levels, omega = 0.5']); 
dpoints = [];
tpoints = linspace(1,200,200);
for (i=1:200)
  dpoints(i) = sin(0.5*(0+i)+0.1*(0.5-normrnd(0.0,1.0))) + 0.7*normrnd(0.0,1.0);
endfor
figure; # Plot time series
plot(tpoints,dpoints);
title("Example time series 2","fontsize",16);
xlabel("t","fontsize",14);
ylabel("d(t)","fontsize",14);
set(gca(),"xlim",[1, 200],"ylim",[-2.2, 2.2]);
disp(['Press enter to start BSA, with the same settings as before']); pause
fprintf (stderr, "Working... This may take a few minutes.\n");
[normp,omega,omega_max] = bsa_auto (dpoints',"Test2",0,1,100,0,1,"BSA-DEMO");
disp(['Done, estimated omega: ',num2str(omega_max)]);
disp(['Press enter to plot omega and posterior density']); pause
plot(omega,normp,"1","linewidth",1);# Plot results
title("BSA result 2","fontsize",16);
xlabel('\omega',"fontsize",14);
ylabel("PDF","fontsize",14);
set(gca(),"xlim",[0, 1]);
disp(['Press enter to go on to the next example']); pause
############################## Example 3
disp(['Case3: A time series with background trends, omega = 0.5']); 
dpoints = [];
tpoints = linspace(1,200,200);
for (i=1:200)
  dpoints(i) = sin(0.5*i)- (i^2*0.005)  + 0.1*normrnd(0.0,1.0);
endfor
figure; # Plot time series
plot(tpoints,dpoints);
title("Example time series 3","fontsize",16);
xlabel("t","fontsize",14);
ylabel("d(t)","fontsize",14);
disp(['To start BSA, we enter the following on the command line:']); 
disp(['[normp,omega,maxomega] = bsa_auto (dpoints,"Test3",0,1,100,4,1,"BSA-DEMO");']); 
disp(['where dpoints is the data, as a single column,']); 
disp(['Test3 is the name of the time series,']);
disp(['0 is the lower limit of omega and 1 is the upper,']);
disp(['100 is the number of sample points between the lower and upper omega,']);
disp(['4 is the upper limit of background functions,']);
disp(['1 is the interval between sampling points in the data (in seconds),']);
disp(['and "BSA-DEMO" is the name of the output file.']); 
disp(['Press enter to start BSA']); pause
fprintf (stderr, "Working... This may take a few minutes.\n");
[normp,omega,omega_max] = bsa_auto (dpoints',"Test3",0,1,100,4,1,"BSA-DEMO");
disp(['Done, estimated omega: ',num2str(omega_max)]);
disp(['Press enter to plot omega and posterior density']); pause
plot(omega,normp,"1","linewidth",1);# Plot results
title("BSA result 3","fontsize",16);
xlabel('\omega',"fontsize",14);
ylabel("PDF","fontsize",14);4
set(gca(),"xlim",[0, 1]);
disp(['Press enter to go on to the next example']); pause
##############################  Example 4
disp(['Case 4: A time series with two frequencies, omega = 0.5 and 0.3']); 
dpoints = [];
tpoints = linspace(1,200,200);
for (i=1:200)
  dpoints(i) = sin(0.5*i) + sin(0.3*i)   + 0.1*normrnd(0.0,1.0);
endfor
figure; # Plot time series
plot(tpoints,dpoints);
title("Example time series 4","fontsize",16);
xlabel("t","fontsize",14);
ylabel("d(t)","fontsize",14);
disp(['Now we will calculate the evidence, since there seems to be several frequencies:']); 
disp(['To start BSA, we enter the following on the command line:']); 
disp(['[logZ logZerror momega stomega] = bsa_nest (dpoints,"Test4",0,1,100,0,1,"BSA-DEMO",500);']); 
disp(['where dpoints is the data, as a single column,']); 
disp(['Test4 is the name of the time series,']);
disp(['0 is the lower limit of omega and 1 is the upper,']);
disp(['100 is the number of sample points between the lower and upper omega,']);
disp(['0 is the upper limit of background functions,']);
disp(['1 is the interval between sampling points in the data (in seconds),']);
disp(['"BSA-DEMO" is the name of the output file.']);
disp(['and 500 is the number of samples drawn from the posterior to calculate the evidence.']); 
disp(['Press enter to start BSA']); pause
fprintf (stderr, "Working... This may take a few minutes.\n");
[logZ logZerror momega stomega] = bsa_nest(dpoints',"Test3",0,1,100,0,1,"BSA-DEMO",500);
disp(['Done, estimated omega: ',num2str(momega)]);
disp(['Estimated evidence: ',num2str(logZ)]);
disp(['Demo finished!']);
disp(['Results have been written to file BSA-DEMO.']);