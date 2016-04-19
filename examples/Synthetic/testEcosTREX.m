% TEMPLATE for a script that compares statistical performance of the
% cTREX (Ecos) and qTREX(PSS) versions

% Dimension = number of columns of X
p = 100;
disp(['Dimension of predictors: ',num2str(p)])

% Noise
noiseType = {'norm'};
noiseInd = 1;

% Sample sizes
n = 50;
numSamples = n;
disp(['Number of data points: ',num2str(n)])

% Correlation vector
kappaVec = [0,0.9];
numKappa = length(kappaVec);

% Noise level
sigVec = [0.1,0.5,3];
numSig = length(sigVec);

% Number of repetitions
numRep = 21;

runTimeMat = zeros(numKappa,numSig,numRep,2);
estMat = zeros(numKappa,numSig,numRep,2);
predMat = zeros(numKappa,numSig,numRep,2);
tpMat = zeros(numKappa,numSig,numRep,2);
fpMat = zeros(numKappa,numSig,numRep,2);
fnMat = zeros(numKappa,numSig,numRep,2);

betaCellQ = cell(numKappa,numSig,numRep);
betaCellC = cell(numKappa,numSig,numRep);

funCellQ = cell(numKappa,numSig,numRep);
funCellC = cell(numKappa,numSig,numRep);

XCell = cell(numKappa,numSig,numRep);
YCell = cell(numKappa,numSig,numRep);

% Number of non-zero entries of the regression vector
nnzs=5;

% Generate leading non-zero entries of size one
firstEntries = ones(nnzs,1);

% True beta vector
betaTrue = [firstEntries;zeros(p-nnzs,1)]; % only nnz nonzero coefficients

% Regularization parameter c for TREX
cVal = 0.5;

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
            
            disp([' '])
            disp(['########'])
            disp(['Rep: ',num2str(r),' with kappa=',num2str(kappa),' sigma=',num2str(sig)])
            disp(['########'])            
            
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
            
            % Number of random restarts
            nTREXRep = 1;
            
            % Non-convex TREX options
            trexopts.boots = L;
            trexopts.rep = nTREXRep;
            trexopts.verbose = 0;
            % HERE YOU CAN ALSO PROVIDE A MATRIX OF START VECTORS. IF YOU DO
            % NOT, the TREX WILL GENERATE INTERNALLY nTREXRep START VECTORS
            % WITH THE ZERO ONE ALWAYS INCLUDED
            trexopts.beta0 = zeros(p,1); % Include the zero vector as start point
            
            trexopts.q = 40;              % Order of the norm to approximate max norm
            trexopts.normConst = cVal;    % Constant c
            
            % Single thread version
            tic;[betaQ,betaMatQ,funVecQ] = trex(X,Y,trexopts); solveQTREX = toc;
            
            % Multi-thread version
            %tic;[betaQ,betaMatQ,funVecQ] = trexp(X,Y,trexopts); solveQTREX = toc;
            
            % beta vector entries are only approximately zero due to the solver; set small values to 0
            epsPssg = 1e-2;
            betaQ(abs(betaQ)<epsPssg) = 0;
            betaMatQ(betaMatQ(:)<epsPssg) = 0;            
            
            % Convex TREX options for ECOS solver
            trexEcosopts.cpath = cVal;          % Constant c
            trexEcosopts.verbose = 0;           % No diagnostic output
            trexEcosopts.feastol = 1e-6;        % Tolerances in the solver
            trexEcosopts.reltol = 1e-6;
            trexEcosopts.abstol = 1e-6;
            
            % Single thread version
            tic;[betaC,betaMatC,funVecC] = trex_ecos(X,Y,trexEcosopts); solveCTREX = toc;
            
            % Multi-thread version
            %tic;[betaC,betaMatC,funVecC] = trex_ecosp(X,Y,trexEcosopts); solveCTREX = toc;
            
            % beta vector entries are only approximately zero due to the solver; set small values to 0
            epsEcos = 1e-2;
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
            runTimeMat(k,s,r,1) = solveQTREX;
            runTimeMat(k,s,r,2) = solveCTREX;
            
            estMat(k,s,r,1) = estTrexQ;
            estMat(k,s,r,2) = estTrexC;
            
            predMat(k,s,r,1) = predTrexQ;
            predMat(k,s,r,2) = predTrexC;
            
            
            tpMat(k,s,r,1) = tpTrexQ;
            tpMat(k,s,r,2) = tpTrexC;
            
            fpMat(k,s,r,1) = fpTrexQ;
            fpMat(k,s,r,2) = fpTrexC;
            
            fnMat(k,s,r,1) = fnTrexQ;
            fnMat(k,s,r,2) = fnTrexC;
            
            betaCellQ{k,s,r} = betaMatQ;
            betaCellC{k,s,r} = betaMatC;
            
            funCellQ{k,s,r} = funVecQ;
            funCellC{k,s,r} = funVecC;
            
            XCell{k,s,r} = X;
            YCell{k,s,r} = Y;
            
            % Display output
            disp(['Non-zero betas from qTREX: ',mat2str(betaQ(betaQ~=0))])
            disp(['Non-zero betas from cTREX: ',mat2str(betaC(betaC~=0))])
            
            disp(['qTREX TP/FP/FN: ',num2str(tpTrexQ),' ',num2str(fpTrexQ),' ',num2str(fnTrexQ)]);
            disp(['cTREX TP/FP/FN: ',num2str(tpTrexC),' ',num2str(fpTrexC),' ',num2str(fnTrexC)]);
            
            disp(['RunTime qTREX, cTREX: ',num2str(solveQTREX),' ',num2str(solveCTREX)]);
            
            disp(' ')

        end
    end
end


%save(['ICMLExampleEcos_',date,'_p',num2str(p),'_n',num2str(n),'_nRep',num2str(numRep)])


