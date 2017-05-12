% This demo plots the target FDR vs the median FDP when the knockoff
% filter (LASSO and TREX) is applied to a synthetic data set.

%% Synthetic problem parameters

n = 101;          % Number of data points
p = 50;          % Number of variables
k = 15;           % Number of variables with nonzero coefficients
amplitude = 3.5;  % Magnitude of nonzero coefficients
sigma = 1;        % Noise level

% rng(45678);       % Random seed

%% Synthetic problem construction

X = randn(n,p) / sqrt(n);
S0 = randsample(p,k);
beta = zeros(p,1);
beta(S0) = amplitude;
sampleY = @() X*beta + sigma .* randn(n,1);

%% Plot target FDR vs achieved FDR

FDP = @(S) sum(beta(S) == 0) / max(1, length(S));

ntrials = 20;
q = 0.05:0.1:0.55;
fdp = zeros(length(q), ntrials);
fdp_trex = zeros(length(q), ntrials);

for i = 1:length(q)
    disp(['Current target FDR: ',num2str(q(i))])
    for j = 1:ntrials
        y = sampleY();

        S = knockoff.filter(X, y, q(i));
        S_TREX = knockoff.filter(X, y, q(i),'Statistic', @trexKnockoffStatistic);

        fdp(i,j) = FDP(S);
        fdp_trex(i,j) = FDP(S_TREX);
    end
end

figure;
plot(q, median(fdp,2),'LineWidth',5);
hold on
plot(q, median(fdp_trex,2),'LineWidth',5);
grid on
xlabel('Target FDP'), ylabel('Median FDP'), title('False Discovery Rate (');
legend('LASSO','TREX','Location','NorthWest')
xlim([0 max(q)]), ylim([0 inf]);
line = refline(1,0);
set(line, 'LineStyle', ':', 'Color', 'black');
