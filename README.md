Tuning-free sparse linear regression with the TREX
=========

This is the TREX MATLAB package for sparse tuning-free linear regression.



## Package structure 
The TREX package contains the following files and folders

- examples/ (different scenarios), also includes figure creation for [2].
- install_trex.m (script to add package to MATLAB path)
- misc/ (different files, the barweb plotting and the knockoff filter (after download))
- solvers/ (Schmidt's PSG code, the ECOS solver (after download), and SCS solver)
- trex/ (TREX files (both single and multi-thread versions), TREX knockoff filter)

The proximal solvers are fully integrated and do not rely on external software.

## Dependencies

### External solvers
Two solvers have been tested to solve the TREX problem in SOCP form.

- ecos: Conic solver for the cTREX
The software can be downloaded [here](https://github.com/embotech/ecos). The MATLAB interface can be found [here](https://github.com/embotech/ecos-matlab).

- SCS: SOCP solver for the cTREX 
SCS can be downloaded [here](https://github.com/cvxgrp/scs). Place the different files in the solver/ folder.

### Knockoff filter
Knockoff filtering with the TREX requires the MATLAB knockoff filter package 
by Barber-Foygel and Candes. The software can be downloaded at:

http://web.stanford.edu/~candes/Knockoffs/package_matlab.html

Please place it in the misc/ folder


## FAQ

USEFUL FOR SCS CODE COMPILATION

Workaround for compiling code using mex for Mac OS X with MATLAB 2015a and Xcode 7+
http://www.mathworks.com/matlabcentral/answers/246507-why-can-t-mex-find-a-supported-compiler-in-matlab-r2015b-after-i-upgraded-to-xcode-7-0


## References 

The code builds on results from the following papers:

* [1] J Lederer, CL Müller, [Don't fall for tuning parameters: tuning-free variable selection in high dimensions with the TREX](https://www.aaai.org/ocs/index.php/AAAI/AAAI15/paper/viewPaper/9359), Twenty-Ninth AAAI Conference on Artificial Intelligence, 2015
* [2] J Bien, I Gaynanova, J Lederer, CL Müller, [Non-convex global minimization and false discovery rate control for the TREX](https://amstat.tandfonline.com/doi/abs/10.1080/10618600.2017.1341414#.XnPaQi2ZNgc), Journal of Computational and Graphical Statistics 27 (1), 23-33, 2018
* [3] PL Combettes, CL Müller, Perspective functions: [Proximal calculus and applications in high-dimensional statistics](https://www.sciencedirect.com/science/article/pii/S0022247X16308071), Journal of Mathematical Analysis and Applications 457 (2), 1283-1306, 2018
* [4] J Bien, I Gaynanova, J Lederer, CL Müller, [Prediction error bounds for linear regression with the TREX](https://link.springer.com/article/10.1007/s11749-018-0584-4), Test 28 (2), 451-474, 2019

Maintainer:
* Christian L. Müller, Center for Computational Mathematics, Flatiron Institute, Simons Foundation (cmueller@flatironinstitute.org)

