% Test script that analyzes the designed TREX knockoff statistic

% This demo illustrates the basic and advanced usage of the knockoff package 
% including the TREX on a synthetic data set.

%% Synthetic problem parameters

n = 101;          % Number of data points
p = 50;           % Number of variables
k = 20;           % Number of variables with nonzero coefficients
amplitude = 3.5;  % Magnitude of nonzero coefficients
sigma = 1;        % Noise level
q = 0.20;         % Target false discovery rate (FDR)

%rng(456789);       % Random seed

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

y = sampleY();

%% Running the knockoff filter steps manually

% The main function 'knockoff.filter' is a wrapper around simpler functions
% that create knockoffs, compute test statistics, and perform variable
% selection. When more control is necessary, these functions may be
% called directly. We demonstrate this below in reproducing the plot of
% Figure 1.

X_ko = knockoff.create(X);

tic
[Wf,Zf] = knockoff.stats.forwardSelection(X, X_ko, y);
tf = knockoff.threshold(Wf, q);
timeFS = toc

tic
[Wt,Zt,trexTrace,trexopts] = trexKnockoffStatistic(X,X_ko,y);
tt = knockoff.threshold(Wt, q);
timeTREX = toc

tic
[Wl,Zl] = knockoff.stats.lassoSignedMax(X, X_ko, y);
tl = knockoff.threshold(Wl, q);
timeLASSO = toc


%% Plot TREX path
cpath = trexopts.cpath;

figure;
plot(cpath,trexTrace')
grid on
xlabel('c path')
title('Trace plots of the coefficients fit by qTREX')
%set(gca,'XTick',1:length(cpath))
%set(gca,'XTickLabel',cpath)


 