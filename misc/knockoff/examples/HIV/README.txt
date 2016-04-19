HIV Drug Resistance Demo
------------------------

Requirements: MATLAB R2013b or later

This demo illustrates the use of the knockoff filter on a real (i.e., not
simulated) experiment, including all the pre- and post-processing steps.
The data set consists of drug resistance measurements and patient genotype
information for several classes of HIV drugs. See Section 4 of the knockoff
filter paper for a more detailed description of the data set.

The main entry point for this example is the script 'main.m'. It fetches the
data from the web, pre-processes it, run the knockoff analysis, performs some
minor post-processing, and plots the results.
