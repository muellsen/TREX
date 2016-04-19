% This demo plots the target FDR vs the average FDP when the knockoff
% filter is applied to a synthetic data set.

%% Synthetic problem parameters

n = 300;          % Number of data points
p = 100;          % Number of variables
k = 15;           % Number of variables with nonzero coefficients
amplitude = 3.5;  % Magnitude of nonzero coefficients
sigma = 1;        % Noise level

rng(45678);       % Random seed

%% Synthetic problem construction

X = randn(n,p) / sqrt(n);
S0 = randsample(p,k);
beta = zeros(p,1);
beta(S0) = amplitude;
sampleY = @() X*beta + sigma .* randn(n,1);

%% Plot target FDR vs achieved FDR

FDP = @(S) sum(beta(S) == 0) / max(1, length(S));

ntrials = 20;
q = 0.05:0.05:0.5;
fdp = zeros(length(q), ntrials);

for i = 1:length(q)
    for j = 1:ntrials
        y = sampleY();
        S = knockoff.filter(X, y, q(i));
        fdp(i,j) = FDP(S);
    end
end

plot(q, mean(fdp,2));
xlabel('Target FDR'), ylabel('Average FDP'), title('False Discovery Rate');
xlim([0 max(q)]), ylim([0 inf]);
line = refline(1,0);
set(line, 'LineStyle', ':', 'Color', 'black');