% Script that creates the Wainwright example used in
%
% Prediction Error Bounds for Linear Regression With the TREX
% Jacob Bien, Irina Gaynanova, Johannes Lederer, Christian Müller

% Number of replicates
nReps = 51;

% Store generated data for post analysis
XCell = cell(nReps,1);
YCell = cell(nReps,1);
noiseVecCell = cell(nReps,1);

% Store TREX solutions
betaTREXCell = cell(nReps,1);
betaTREXMatCell = cell(nReps,1);
funTREXCell = cell(nReps,1);
optsTREXCell = cell(nReps,1);

% Store Lasso solutions
betaLassoCell = cell(nReps,1);
infoLassoMat = cell(nReps,1);

for i=1:nReps
    
    disp(['Repetition ',num2str(i),' out of ',num2str(nReps)])
    
    %% Simulation setup (following Wainwright, Sharp oracle... 2009)
    
    % Dimension = number of columns of X
    p = 32;
    disp(['Dimension of predictors: ',num2str(p)])
    
    % Number of non-zero entries of the regression vector
    nnzs=ceil(0.40*p^0.75);
    disp(['Number of nnzs: ',num2str(nnzs)])
    
    % Sample sizes
    alpha = 2;
    n = round(alpha*nnzs*log(p-nnzs));
    disp(['Number of data points: ',num2str(n)])
    
    % Generate leading non-zero entries of size one
    const = 2;
    firstEntries = const*[(-1).^(1:nnzs)]';
    
    % True beta vector
    betaTrue = [firstEntries;zeros(p-nnzs,1)]; % only nnz nonzero coefficients
    
    % Correlation kappa
    kappa = 0;
    
    % Generate covariance matrix
    covMat = kappa*ones(p,p);
    covMat(1:p+1:p^2) = 1;
    cholCov = chol(covMat);
    
    %  Noise scalar
    sig = 0.5;
    
    % Generate data
    X0 = (cholCov'*randn(p,n))';
    
    % Normalize X to length sqrt(n)
    normX0 = repmat(sqrt(sum(X0.^2)),n,1);
    X = sqrt(n)*X0./normX0;
    
    normX = repmat(sqrt(sum(X.^2)),n,1);
    
    % Gaussian noise
    noiseVec = randn(n,1);
    
    % Response with sigma * standardized noise vector
    % Here: do not use the normalization for numerical reasons
    Y = X*betaTrue + sig*noiseVec;
    %X=X0;
    
    % Save generated data and noise
    XCell{i}=X;
    YCell{i}=Y;
    noiseVecCell{i} = noiseVec;
    
    % Compute Lasso first
    % Do not allow an intercept
    % Fix the number of lambda values
    
    nRegParamsInit = 50;
    lamRatio = 2e-3
    tic;[betaLasso,infoLasso] = lasso(X,Y,'Standardize',0,'CV',5,'NumLambda',nRegParamsInit,'LambdaRatio',lamRatio); solveLasso = toc
    
    betaLassoCell{i} = betaLasso;
    infoLassoMat{i} = infoLasso;
    
    % Upper bound on Lasso path
    u_max = max(abs(X'*Y));
    
    % Percentage of lambda in terms of u_max
    relLassoPath = (infoLasso.Lambda*n)./u_max;
    
    % Lasso may reduce number of lambdas internally
    nRegParams = length(relLassoPath);
    
    % Find the upper bound for the TREX
    [c_max,cVec,u_hatVec,u_max] = findTREXUB(X,Y)
    
    % Corresponding TREX path (same number of c's relative to u_max)
    cpath = c_max.*relLassoPath;
    nCParams = length(cpath);
    
    %% Solution of TREX for approximately valid path (or c=1/2)
    
    % Convex TREX options for ECOS solver
    clear trexEcosopts
    trexEcosopts.cpath = cpath;         % Constant c
    trexEcosopts.verbose = 0;           % No diagnostic output
    trexEcosopts.feastol = 1e-12;       % Tolerances in the solver
    trexEcosopts.reltol = 1e-12;
    trexEcosopts.abstol = 1e-12;
    
    % Multi-thread version
    tic;[betaTREX,betaTREXMat,funVecTREX] = trex_ecosp(X,Y,trexEcosopts); solveCTREX = toc
    
    % Threshold numerical solutions
    betaTREX(abs(betaTREX(:))<1e-6)=0;
    
    betaTREXCell{i} = betaTREX;
    betaTREXMatCell{i} = betaTREXMat;
    funTREXCell{i} = funVecTREX;
    optsTREXCell{i} = trexEcosopts;
    
end

% Save workspace
saveString = ['Wainwright_p',num2str(p),'_alpha',num2str(alpha),'_norm_kap_',num2str(kappa),'_const',num2str(const),'_nReps',num2str(nReps)];
save(saveString);

