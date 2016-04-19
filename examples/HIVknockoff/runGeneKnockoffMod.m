function [discoveries,fdp,discoveries_bh,fdp_bh,nData,pData] = ...
    runGeneKnockoffMod(geneData, drugData, tsmData,statMethod,targetFDP)
%% Select variables

q = targetFDP; % Target FDP.

% Perform variable selection separately for each drug.
nDrugs = size(drugData,2);
nData = zeros(nDrugs,1);
pData = zeros(nDrugs,1);
S = cell(nDrugs, 1);
S_bh = cell(nDrugs, 1);

for i = 1:nDrugs
    
    % Remove rows with missing drug data. This varies from drug to drug.
    y = drugData{:,i};
    missing = isnan(y);
    y = log(y(~missing)); % Log-transform response data.
    X = geneData(~missing,:);
    
    % Remove predictors that appear less than `X_thresh` times.
    X_thresh = 3;
    X = X(:,sum(table2array(X)) >= X_thresh);
    
    % Remove duplicated predictors. Without this step, the design matrix
    % will be rank-deficient, which is bad for reproducibility.
    X(:,sum((corr(table2array(X))-1).^2 < 1e-8) > 1) = [];
    nData(i) = size(X,1);
    pData(i) = size(X,2);
    
    % Run the knockoff filter.
    S{i} = knockoff.filter(X, y, q, 'Knockoffs', 'equi','Statistic',statMethod);
    
    % Compare to Benjamini-Hochberg (on least-squares coefficients).
    lm = fitlm(table2array(X), y, 'Intercept', false);
    pvals_LS = lm.Coefficients.pValue;
    khat_bh = max([0; find(sort(pvals_LS) <= q*(1:size(X,2))'/size(X,2))]);
    if (khat_bh == 0)
        S_bh{i} = [];
    else
        S_bh{i} = X.Properties.VariableNames(...
            pvals_LS <= q*khat_bh/size(X,2));
    end
    
    % Display progress
    disp(['Percentage of drugs checked: ',num2str(i./nDrugs*100)])
    
end

%% Compute FDP

discoveries = zeros(1, nDrugs);
falseDiscoveries = zeros(1, nDrugs);
discoveries_bh = zeros(1, nDrugs);
falseDiscoveries_bh = zeros(1, nDrugs);

for i = 1:nDrugs
    % Compare gene position only (not mutation type).
    positions = cellfun(@(cell) str2double(cell{1}), ...
        regexp(S{i},'P(\d+)','tokens'));
    positions = unique(positions); % Remove possible duplicates.
    discoveries(i) = length(positions);
    falseDiscoveries(i) = length(setdiff(positions, tsmData));
    
    % Repeat for BH.
    positions = cellfun(@(cell) str2double(cell{1}), ...
        regexp(S_bh{i},'P(\d+)','tokens'));
    positions = unique(positions); % Remove possible duplicates.
    discoveries_bh(i) = length(positions);
    falseDiscoveries_bh(i) = length(setdiff(positions, tsmData));
end

fdp = falseDiscoveries ./ max(1, discoveries);
fdp_bh = falseDiscoveries_bh ./ max(1, discoveries_bh);

end