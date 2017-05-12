% Compares c-trex and q-trex in the setting when p>n

% Run install_trex
install_trex

% Add glmnet_matlab in case using random restarts from lasso path
% addpath glmnet_matlab

% Dimension = number of columns of X
pVec = [100,500]; %%%% Add one more dimension
numP=length(pVec);
disp(['Dimension of predictors: ',num2str(pVec)])

% Noise
noiseType = {'norm'};
noiseInd = 1;

% Sample size is fixed at n<p
n = 50;
numSamples = n;
disp(['Number of data points: ',num2str(n)])

% Correlation vector
kappaVec = [0,0.3,0.6,0.9];
numKappa = length(kappaVec);

% Noise level
sigVec = [0.1,0.5,3];
numSig = length(sigVec);

% Number of repetitions
numRep = 21;

runTimeMat = zeros(numP,numKappa,numSig,numRep,2);
estMat = zeros(numP,numKappa,numSig,numRep,2);
predMat = zeros(numP,numKappa,numSig,numRep,2);
tpMat = zeros(numP,numKappa,numSig,numRep,2);
fpMat = zeros(numP,numKappa,numSig,numRep,2);
fnMat = zeros(numP,numKappa,numSig,numRep,2);

betaCellQ = cell(numP,numKappa,numSig,numRep);
betaCellC = cell(numP,numKappa,numSig,numRep);

funCellQ = cell(numP,numKappa,numSig,numRep);
funCellC = cell(numP,numKappa,numSig,numRep);

XCell = cell(numP,numKappa,numSig,numRep);
YCell = cell(numP,numKappa,numSig,numRep);

% Number of non-zero entries of the regression vector
nnzs=5;

% Generate leading non-zero entries of size one
firstEntries = ones(nnzs,1);

% True beta vector
betaTrueVec = cell(numP);

for pind=1:numP
    betaTrueVec{pind}=[firstEntries;zeros(pVec(pind)-nnzs,1)]; % only nnz nonzero coefficients
end
    

% Loop over values of p

for pind=1:numP

    p= pVec(pind);

    % True beta vector
    betaTrue = betaTrueVec{pind}; % only nnz nonzero coefficients

    % Loop over correlation vector kappa
    for k = 1:numKappa
    
        kappa = kappaVec(k);
    
        % Generate covariance matrix
        covMat = kappa*ones(p,p);
        covMat(1:p+1:p^2) = 1;
        cholCov = chol(covMat);
    
        % loop over noise vector
        for s = 1:numSig

            sig = sigVec(s);
        
            % loop over number of repetitions
            for r = 1:numRep
            
                disp(['Rep: ',num2str(r),' with kappa=',num2str(kappa),' sigma=',num2str(sig)])
            
                % Generate data in rth repetition
            
                X = (cholCov'*randn(p,n))';
            
                % Normalize X to length sqrt(n)
                normX = repmat(sqrt(sum(X.^2)),n,1);
                X = sqrt(n)*X./normX;
            
                % Gaussian noise
                if noiseInd==1
                     noiseVec = randn(n,1);
                     % T-distribution noise
                elseif noiseInd==2
                     noiseVec = trnd(4,n,1)./sqrt(2);
                     % Chi-2 noise
                elseif noiseInd==3
                     noiseVec =  exprnd(1,n,1)-1;
                else
                     error('No corrent noise model specified')
                end
            
                % Response with sigma * standardized noise vector
                Y = X*betaTrue + sig*noiseVec;
                
                
                
            
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % TREX regression
            
                % Number of bootstrap runs
                L = 0;
            
                
            
                %%%%%%%%%%%%%%%%%%%
                % Lasso
                nl=20; % number of lasso restarts for trex, 1st one is zlways zero
                opts_lasso=[];
                opts_lasso.nlambda=nl; 
                fit=glmnet(X,Y,[],opts_lasso);
                betalasso=fit.beta;
                
                % Number of random restarts
                nTREXRep = 40; %first zero, then 19 lasso, then 20 random
            
                % Non-convex TREX options
                trexopts.boots = L;
                trexopts.rep = nTREXRep;
                trexopts.verbose = 0;
                % HERE YOU CAN ALSO PROVIDE A MATRIX OF START VECTORS. IF YOU DO
                % NOT, the TREX WILL GENERATE INTERNALLY nTREXRep START VECTORS
                % WITH THE ZERO ONE ALWAYS INCLUDED
                trexopts.beta0 = betalasso; % 1st lasso is always zero
                
                trexopts.q = 40;             % Order of the norm to approximate max norm
                trexopts.normConst = 0.5;    % Constant c
            
                % Single thread version
                %tic;[betaQ,betaMatQ,funVecQ] = trex(X,Y,trexopts); solveQTREX = toc
            
                % Multi-thread version
                tic;[betaQ,betaMatQ,funVecQ] = trexp(X,Y,trexopts); solveQTREX = toc
            
                % beta vector entries are only approximately zero due to the solver; set small values to 0
                epsPssg = 1e-5;
                betaQ(abs(betaQ)<epsPssg) = 0;
                betaMatQ(betaMatQ(:)<epsPssg) = 0;
            
                % Convex TREX options
                trexEcosopts.c = 0.5;        % Constant c
                trexEcosopts.verbose = 0;    % No diagnostic output
                inopts.feastol = 1e-6;       % Tolerances in the solver
                inopts.reltol = 1e-6;
                inopts.abstol = 1e-6;
            
                % Single thread version
                %tic;[betaC,betaMatC,funVecC] = trex_ecos(X,Y,trexEcosopts); solveCTREX = toc
            
                % Multi-thread version
                tic;[betaC,betaMatC,funVecC] = trex_ecosp(X,Y,trexEcosopts); solveCTREX = toc
            
                % beta vector entries are only approximately zero due to the solver; set small values to 0
                epsEcos = 1e-5; % Irina: potentially make this larger
                betaC(abs(betaC)<epsEcos) = 0;
                betaMatC(betaMatC(:)<epsEcos) = 0;
            
                % Estimation of beta
                estTrexQ = sqrt(sum((betaQ-betaTrue).^2));
                estTrexC = sqrt(sum((betaC-betaTrue).^2));
            
                % Prediction of Y
                predTrexQ = sum((X*betaTrue-X*betaQ).^2)/n;
                predTrexC = sum((X*betaTrue-X*betaC).^2)/n;
            
                % Recovery of beta entries (just non-zero pattern)
                tpTrexQ = length(intersect(find(betaQ~=0),find(betaTrue~=0)));
                fpTrexQ = length(intersect(find(betaQ~=0),find(betaTrue==0)));
                fnTrexQ = length(intersect(find(betaQ==0),find(betaTrue~=0)));
            
                tpTrexC = length(intersect(find(betaC~=0),find(betaTrue~=0)));
                fpTrexC = length(intersect(find(betaC~=0),find(betaTrue==0)));
                fnTrexC = length(intersect(find(betaC==0),find(betaTrue~=0)));
            
            
                % Save all in matrix
                runTimeMat(pind,k,s,r,1) = solveQTREX;
                runTimeMat(pind,k,s,r,2) = solveCTREX;
            
                estMat(pind,k,s,r,1) = estTrexQ;
                estMat(pind,k,s,r,2) = estTrexC;
            
                predMat(pind,k,s,r,1) = predTrexQ;
                predMat(pind,k,s,r,2) = predTrexC;
            
                tpMat(pind,k,s,r,1) = tpTrexQ;
                tpMat(pind,k,s,r,2) = tpTrexC;
            
                fpMat(pind,k,s,r,1) = fpTrexQ;
                fpMat(pind,k,s,r,2) = fpTrexC;
            
                fnMat(pind,k,s,r,1) = fnTrexQ;
                fnMat(pind,k,s,r,2) = fnTrexC;
            
                betaCellQ{pind,k,s,r} = betaMatQ;
                betaCellC{pind,k,s,r} = betaMatC;
            
                funCellQ{pind,k,s,r} = funVecQ;
                funCellC{pind,k,s,r} = funVecC;
            
                XCell{pind,k,s,r} = X;
                YCell{pind,k,s,r} = Y;
            
                end
            end
    end
end

save(['SmallNExample_',date,'_p',num2str(pVec),'_n',num2str(n),'_nRep',num2str(numRep)])


