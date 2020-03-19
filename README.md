PerspeCtive M-estimation (PCM) package
=========

This is the PCM MATLAB package for perspective M-estimation.
The package introduces an optimization model for maximum likelihood-type estimation (M-estimation)
that generalizes a large class of known statistical models, including Huber’s concomitant M-estimation model,
the scaled Lasso, &nu;-Support Vector Machine Regression, and penalized estimation with structured sparsity.
The model, termed perspective M-estimation, leverages the observation that a wide class of
convex M-estimators with concomitant scale as well as structured norms are instances of perspective functions.

The code builds on results from the following papers:

* [1] P. L. Combettes and C. L. Müller, [Perspective functions: Proximal calculus and applications in high-dimensional statistics](https://www.sciencedirect.com/science/article/pii/S0022247X16308071), J. Math. Anal. Appl., vol. 457, no. 2, pp. 1283–1306, 2018.
* [2] P. L. Combettes and C. L. Müller, [Perspective M-estimation via proximal decomposition](https://arxiv.org/abs/1805.06098), Electronic Journal of Statistics, 2020, [Journal version](https://projecteuclid.org/euclid.ejs/1578452535) 
* [3] P. L. Combettes and C. L. Müller, [Regression models for compositional data: General log-contrast formulations, proximal optimization, and microbiome data applications](https://arxiv.org/abs/1903.01050), arXiv, 2019.

Developer:
* Christian L. Müller, Center for Computational Mathematics, Flatiron Institute, Simons Foundation (cmueller@flatironinstitute.org)

## Installation ##

The PCM package is self-contained. No external software needed. However, for testing the code base we rely
on the [cvx package](http://cvxr.com/cvx/).

After downloading the PCM package, use

```MATLAB
% This will add the folders to your MATLAB path
add_pcm
```
to add all subfolders to your MATLAB path. 

## Package structure ##

The PCM package comprises the following folders that contain different
functions and scripts: 

- The examples/ folder contains several test cases about the different modes of usage.
Please refer to the README.md in the folder for further information.

- The prox/ folder implements projection and proximity operators for several perspective functions and
standard regularization functions and set indicators.

- The solvers/ folder implements a generalized Douglas-Rachford scheme for perspective M-estimations. The function pcmC2.m
is the current standard solver. 

- The sqrtlasso-solver/ folder implements coordinate descent solvers for the SQRT-Lasso and the scaled Lasso that solves a specific variant of the perspective M-estimation with the square-loss. 

- The misc/ folder comprises several helper routines and functions for data transformation and analysis.


## Log-contrast models for compositional data - with microbiome applications ##

In examples/LogContrastModels/ we provide all numerical examples used in [Regression models for compositional data: General log-contrast formulations, proximal optimization, and microbiome data applications](https://arxiv.org/abs/1903.01050).

There, we consider the special but important case of estimating a linear log-contrast model for compositional covariates X 
where each of the n rows comprises p-dimensional compositions (or relative abundances) and n continuous outcome variables 
Y that can also contain outliers o (in form of a (sparse) mean shift). The generative model thus reads: 

<a href="https://www.codecogs.com/eqnedit.php?latex=Y=\log(X)\beta&space;&plus;&space;o&space;&plus;&space;\sigma&space;\epsilon&space;\qquad&space;\text{s.t.}\qquad&space;C^T&space;\beta&space;=&space;0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?Y=\log(X)\beta&space;&plus;&space;o&space;&plus;&space;\sigma&space;\epsilon&space;\qquad&space;\text{s.t.}\qquad&space;C^T&space;\beta&space;=&space;0" title="Y=\log(X)\beta + o + \sigma \epsilon \qquad \text{s.t.}\qquad C^T \beta = 0" /></a>

The folder comprises code and data for reproducing the numerical experiments in [[3]](https://arxiv.org/abs/1903.01050). 

