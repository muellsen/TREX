% This demo illustrates the basic and advanced usage of this package on a
% synthetic data set.

%% Synthetic problem parameters

n = 600;          % Number of data points
p = 200;          % Number of variables
k = 30;           % Number of variables with nonzero coefficients
amplitude = 3.5;  % Magnitude of nonzero coefficients
sigma = 1;        % Noise level
q = 0.20;         % Target false discovery rate (FDR)

rng(45678);       % Random seed

%% Synthetic problem construction

X = randn(n,p) / sqrt(n);
S0 = randsample(p,k);
beta = zeros(p,1);
beta(S0) = amplitude;
sampleY = @() X*beta + sigma .* randn(n,1);

trueDiscoveries = @(S) sum(beta(S) > 0);
FDP = @(S) sum(beta(S) == 0) / max(1, length(S));
printSummary = @(S) fprintf(...
    ['%d true discoveries\n' ...
     'FDP = %2.2f%% (target FDR = %2.f%%)\n'], ...
    trueDiscoveries(S), 100*FDP(S), 100*q);

%% Running the knockoff filter

% Here we call the knockoff filter with all the default settings. We will
% explore some variations below.

y = sampleY();
S = knockoff.filter(X, y, q);
printSummary(S);

%% Using a different method for creating knockoff variables

% By default, equi-correlated knockoff variables are created. It is also
% possible to create optimized knockoff variables that are the solution to
% a semi-definite programming (SDP) problem.

S = knockoff.filter(X, y, q, 'Knockoffs', 'SDP');
printSummary(S);

%% Using a different test statistic

% By default, a test statistic based on the lasso is used. Here we use
% a different statistic based on forward selection.

S = knockoff.filter(X, y, q, 'Statistic', @knockoff.stats.forwardSelection);
printSummary(S);

%% Using a custom test statistic

% It is also possible to define your own test statistic. To illustrate
% this, we implement a very simple statistic from the knockoff paper.

myKnockoffStatistic = @(X, X_ko, y) ...
    abs(X' * y) - abs(X_ko' * y);

S = knockoff.filter(X, y, q, 'Statistic', myKnockoffStatistic);
printSummary(S);

% As another example, we show how to change the number of lambda values
% used to approximate the lasso path in the default test statistic
% (cf. the documentation for knockoff.stats.lassoSignedMax).

myLassoStatistc = @(X, X_ko, y) ...
    knockoff.stats.lassoSignedMax(X, X_ko, y, 10*p);

S = knockoff.filter(X, y, q, 'Statistic', myLassoStatistc);
printSummary(S);

%% Running the knockoff filter steps manually

% The main function 'knockoff.filter' is a wrapper around simpler functions
% that create knockoffs, compute test statistics, and perform variable
% selection. When more control is necessary, these functions may be
% called directly. We demonstrate this below in reproducing the plot of
% Figure 1.

X_ko = knockoff.create(X);
[W,Z] = knockoff.stats.lassoSignedMax(X, X_ko, y);
t = knockoff.threshold(W, q);

fig = figure();
hold on
set(fig, 'DefaultTextInterpreter', 'latex');
gscatter(Z(1:p), Z(p+1:2*p), ismember(1:p, S0), 'kr');
plot([t t 0], [0 t t], 'k');
hold off

xlabel('Value of $\lambda$ when $X_j$ enters model');
ylabel('Value of $\lambda$ when $\tilde X_j$ enters model');
limits = [0 ceil(max(Z))];
xlim(limits); ylim(limits);
title('Knockoff Filter with Lasso Statistic');
legend('Null feature', 'Non-null feature');
line = refline(1,0);
set(line, 'LineStyle', ':', 'Color', 'black');
