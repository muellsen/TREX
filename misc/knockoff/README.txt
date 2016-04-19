Knockoff Filter in MATLAB
=========================

Version: 0.2
Date: 02/04/2015
Authors: Rina Foygel Barber, Emmanuel Candes, Evan Patterson
License: GPL-3

This package is an implementation of the knockoff filter in MATLAB.

Requirements
------------

- MATLAB 2012 or later
- Statistical toolbox
- CVX (optional; required for creating SDP knockoffs)

  Warning: Your CVX release must be from April 2014 or later. Earlier
  versions of CVX contain a serious bug affecting the SDP knockoffs.

Installation
------------

To install the package, unzip this archive anywhere on your machine, say to
'<path>', then add the lines

    addpath('<path>/knockoff_matlab')
    addpath('<path>/knockoff_matlab/glmnet_matlab')

to your MATLAB startup file. (If you already have glmnet_matlab installed, 
the last line is not necessary.)

To test your installation, you can run one of the demo files in the 'examples'
subdirectory. Alternatively, if you have MATLAB 2013a or newer, you can execute
the test suite by typing

    runtests('knockoff.tests')

in the Command Window. All the tests should pass.

Documentation
-------------

For an overview of the functions in this package, run

    help knockoff

For a list of included knockoff statistics, run

    help knockoff.stats

Besides the documentation associated with individual functions 
(accessible via the 'help' function), the main source of documentation
is the collection of examples in the 'examples' subdirectory.

The file 'FirstExamples.m' is a good place to start.
