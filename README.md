<img src="https://i.imgur.com/Ei8KgYG.png" alt="TREX" height="150" align="right"/>

Tuning-free sparse linear regression with the TREX 
=========

This is a resource page for the TREX which allows sparse tuning-free linear regression. Currently, the TREX is available 
as MATLAB package.

## Background

The forward model is assumed to be the standard linear model: 

<img src="https://latex.codecogs.com/gif.latex?y&space;=X\beta&space;+\sigma&space;\epsilon&space;" align="middle"/> 

Here, X is a known design matrix and y is a known continuous response vector. The vector &beta; comprises the unknown coefficients and &sigma; an unknown scale.

The package handles the solution of the TREX objective. 

<img src="https://latex.codecogs.com/gif.latex?\hat{\beta_\text{TREX}} = \arg\min_{\beta\in\mathbb&space;R^p}\left\{\frac{\|Y-X\beta\|_2^2}{c\|X^\top(Y-X\beta)\|_\infty}&plus;\|\beta\|_1\right\}." align="middle"/>

The constant c is typically set to c=1/2, thus requiring no tuning parameter (as compared to the Lasso). However, the objective is non-convex and comprises 2p minima.  

Several different algorithmic strategies are available to solve the objective. A proximal gradient descent for an approximate solution has been introduced in [[1]](#references). Via appropriate reformulation and decomposition, the TREX can be solved exactly by solving 2p Second-order Cone Programs (SOCPs) [[2]](#references). Alternatively, the convex subproblems can be reformulated to be amenable to the Douglas-Rachford proximal splitting algorithm [[3]](#references).

The package includes all of the above algorithmic strategies in one framework.

## Package structure 
The TREX package contains the following files and folders

- examples/ (different scenarios), also includes figure creation for [2].
- solvers/ (Schmidt's PSG code, the ecos solver (after additional download), and SCS solver)
- trex/ (TREX solvers (both single and multi-thread versions), TREX knockoff filter)
- misc/ (additional files including barweb plotting and the knockoff filter (after download))

## Dependencies

### External solvers
Two solvers have been tested to solve the TREX problem in SOCP form.

- ecos: Conic solver for the cTREX
The software can be downloaded [here](https://github.com/embotech/ecos). The MATLAB interface can be found [here](https://github.com/embotech/ecos-matlab).

- SCS: SOCP solver for the cTREX 
SCS can be downloaded [here](https://github.com/cvxgrp/scs). 

The solver packages should be compiled and placed in the solver/ folder.

### Knockoff filter
Knockoff filtering with the TREX requires the MATLAB knockoff filter package 
by Barber-Foygel and Candes. The software can be downloaded [here](https://github.com/msesia/knockoff-filter).

The deprecated link during initial development was [here](http://web.stanford.edu/~candes/Knockoffs/package_matlab.html)

Please place it in the misc/ folder

### Other solvers
The proximal solvers from [[1]](#references) and [[3]](#references) are fully integrated and do not rely on external software.

## Basic example
To include the package in your MATLAB environment, type first

```MATLAB
install_trex
````

## References 

The code builds on results from the following papers:

* [1] J Lederer, CL Müller, [Don't fall for tuning parameters: tuning-free variable selection in high dimensions with the TREX](https://www.aaai.org/ocs/index.php/AAAI/AAAI15/paper/viewPaper/9359), Twenty-Ninth AAAI Conference on Artificial Intelligence, 2015
* [2] J Bien, I Gaynanova, J Lederer, CL Müller, [Non-convex global minimization and false discovery rate control for the TREX](https://amstat.tandfonline.com/doi/abs/10.1080/10618600.2017.1341414#.XnPaQi2ZNgc), Journal of Computational and Graphical Statistics 27 (1), 23-33, 2018
* [3] PL Combettes, CL Müller, Perspective functions: [Proximal calculus and applications in high-dimensional statistics](https://www.sciencedirect.com/science/article/pii/S0022247X16308071), Journal of Mathematical Analysis and Applications 457 (2), 1283-1306, 2018
* [4] J Bien, I Gaynanova, J Lederer, CL Müller, [Prediction error bounds for linear regression with the TREX](https://link.springer.com/article/10.1007/s11749-018-0584-4), Test 28 (2), 451-474, 2019

Maintainer:
* Christian L. Müller, Center for Computational Mathematics, Flatiron Institute, Simons Foundation (cmueller@flatironinstitute.org)

## Known issues

Workaround for compiling SCS code using mex for Mac OS X with MATLAB 2015a and Xcode 7+
http://www.mathworks.com/matlabcentral/answers/246507-why-can-t-mex-find-a-supported-compiler-in-matlab-r2015b-after-i-upgraded-to-xcode-7-0

