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

%% Running the knockoff filter

% Here we call the knockoff filter with all the default settings. We will
% explore some variations below.

y = sampleY();
S_EQILASSO = knockoff.filter(X, y, q);
disp('Standard LASSO with equi-correlated knockoffs')
printSummary(S_EQILASSO);
disp(' ')


%% Using a different method for creating knockoff variables

% By default, equi-correlated knockoff variables are created. It is also
% possible to create optimized knockoff variables that are the solution to
% a semi-definite programming (SDP) problem.

S_SDP = knockoff.filter(X, y, q, 'Knockoffs', 'SDP');
disp('Standard LASSO with SDP knockoffs')
printSummary(S_SDP);
disp(' ')

%% Using a different test statistic

% By default, a test statistic based on the lasso is used. Here we use
% a different statistic based on forward selection.

S_FORWARD = knockoff.filter(X, y, q, 'Statistic', @knockoff.stats.forwardSelection);
disp('Stepwise Forward Selection')
printSummary(S_FORWARD);
disp(' ')

%% Using a custom test statistic

% It is also possible to define your own test statistic. To illustrate
% this, we implement a very simple statistic from the knockoff paper.

myKnockoffStatistic = @(X, X_ko, y) ...
    abs(X' * y) - abs(X_ko' * y);

S_SS = knockoff.filter(X, y, q, 'Statistic', myKnockoffStatistic);
disp('Simple screening statistics')
printSummary(S_SS);
disp(' ')

%% Using the TREX test statistic

% It is also possible to define your own test statistic. To illustrate
% this, we implement a very simple statistic from the knockoff paper.

S_TREX = knockoff.filter(X, y, q, 'Statistic', @trexDifference);
disp('TREX-diff statistics with equi-correlated knockoffs')
printSummary(S_TREX);
disp(' ')

S_TREXSDP = knockoff.filter(X, y, q, 'Knockoffs', 'SDP','Statistic', @trexDifference);
disp('TREX-diff statistics with SDP knockoffs')
printSummary(S_TREXSDP);
disp(' ')

S_TREX = knockoff.filter(X, y, q, 'Statistic', @trexSignedMax);
disp('TREX-SignMax statistics with equi-correlated knockoffs')
printSummary(S_TREX);
disp(' ')

S_TREXSDP = knockoff.filter(X, y, q, 'Knockoffs', 'SDP','Statistic', @trexSignedMax);
disp('TREX-SignMax statistics with SDP knockoffs')
printSummary(S_TREXSDP);
disp(' ')

%% Using the convex TREX
S_TREXfun = knockoff.filter(X, y, q, 'Statistic', @trexKnockoffStatisticC);
disp('TREX-fun statistics with equi-correlated knockoffs')
printSummary(S_TREXfun);
disp(' ')

% As another example, we show how to change the number of lambda values
% used to approximate the lasso path in the default test statistic
% (cf. the documentation for knockoff.stats.lassoSignedMax).


myLassoStatistic = @(X, X_ko, y) ...
    knockoff.stats.lassoSignedMax(X, X_ko, y, 20*p);

S = knockoff.filter(X, y, q, 'Statistic', myLassoStatistic);
disp('LASSO statistics with different path')
printSummary(S);
disp(' ')

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

fig = figure();
hold on
set(fig, 'DefaultTextInterpreter', 'latex');
gscatter(Zf(1:p), Zf(p+1:2*p), ismember(1:p, S0), 'kr');
plot([tf tf 0], [0 tf tf], 'k');
hold off

xlabel('Value of \lambda when X_j enters model');
ylabel('Value of \lambda when \tilde X_j enters model');
limits = [0 ceil(max(Zf))];
xlim(limits); ylim(limits);
title('Knockoff Filter with Forward Selection Statistic');
legend('Null feature', 'Non-null feature');
line = refline(1,0);
set(line, 'LineStyle', ':', 'Color', 'black');


tic
[Wtf,Ztf] = trexKnockoffStatisticC(X,X_ko,y);
tt = knockoff.threshold(Wtf, q);
timeTREX = toc

tic
[Wt,Zt] = trexDifference(X,X_ko,y);
tt = knockoff.threshold(Wt, q);
timeTREX = toc

toc
fig = figure();
hold on
set(fig, 'DefaultTextInterpreter', 'latex');
gscatter(Zt(1:p), Zt(p+1:2*p), ismember(1:p, S0), 'kr');
plot([tt tt 0], [0 tt tt], 'k');
hold off
xlabel('Value of c when X_j enters model');
ylabel('Ranking when \tilde X_j enters model');
limits = [0 ceil(max(Zt))];
xlim(limits); ylim(limits);
title('Knockoff Filter with TREX Statistic');
legend('Null feature', 'Non-null feature');
line = refline(1,0);
set(line, 'LineStyle', ':', 'Color', 'black');


tic
[Wl,Zl] = knockoff.stats.lassoSignedMax(X, X_ko, y);
timeLASSO = toc

tl = knockoff.threshold(Wl, q);

fig = figure();
hold on
set(fig, 'DefaultTextInterpreter', 'latex');
gscatter(Zl(1:p), Zl(p+1:2*p), ismember(1:p, S0), 'kr');
plot([tl tl 0], [0 tl tl], 'k');
hold off

xlabel('Value of \lambda when X_j enters model');
ylabel('Value of \lambda when \tilde X_j enters model');
limits = [0 ceil(max(Zl))];
xlim(limits); ylim(limits);
title('Knockoff Filter with LASSO Statistic');
legend('Null feature', 'Non-null feature');
line = refline(1,0);
set(line, 'LineStyle', ':', 'Color', 'black');

