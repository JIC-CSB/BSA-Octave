
Bayesian Spectrum Analysis (BSA) package in Octave
==================================================
by Emma Granqvist & Richard J. Morris, June 2011.
Maintained by Matthew Hartley.

This code is only proof-of-principle, it is released freely and without any warranty. For more information see Granqvist et al., BMC Systems Biology 2011, and Bretthorst, Springer-Verlag 1988. Please cite the above papers when using this code.


To set up the BSA package, start Octave and type: 

   bsa_setup

at the command line. For a demo of test cases, type: 

   bsa_demo


Key Functions
-------------

* bsa_auto - Wrapper function for the automated model comparison, and stops adding background functions when model ratio > 1

* bsa_expand - Wrapper function for the model comparison but does not stop when ratio > 1, keeps going to the chosen expansion order of background functions

* bsa_nest - Wrapper function for Nested Sampling routine, to calculate the evidence and mean omega. Use this for multiple frequencies present in the data

* bsa_post_norm - Function that returns vector of omega and normalised posterior probability distribution

* bsa_prob - Function that calculates the Student-t distribution



Output
------
A results file is automatically created from the wrapping functions.
It will have two lines per case for bsa_expand and bsa_auto, one of which is the initial sampling over the specified interval, with the chosen number of sampling plots. The second line is the result of sampling locally around the peak of maximum probability.
For bsa_nest, the results file will have one line, with the estimated evidence and omega values.
