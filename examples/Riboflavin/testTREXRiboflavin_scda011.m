% Example Riboflavin dataset from
%
% Bühlmann, P., Kalisch, M., & Meier, L. (2014).
% High-Dimensional Statistics with a View Toward Applications in Biology.
% Annual Review of Statistics and Its Application, 1(1), 255?278.
% doi:10.1146/annurev-statistics-022513-115545

% Initialize global variables
global X
global Y
global normConst
global qNorm

qNorm = 100;
normConst = 0.5;

% Test TREX and B-TREX on the riboflavin data set

% Load relevant data
load riboflavinData.mat

% B. subtilis gene names
geneNames = riboTags;

% Unnormalized response log-production rates
% of riboflavin from n=71 B. subtilis strains
Yunorm = riboData(1,:)';                                                                                                                                   
                                                                                                                                                           
% Center response variable                                                                                                                                 
Y = Yunorm-mean(Yunorm);                                                                                                                                   
                                                                                                                                                           
% Gene expression data for p=4088 genes across all n=71 strains                                                                                            
Xunorm = riboData(2:end,:)';

% Dimension of the problem
[n,p] = size(Xunorm);

% Standardize data
X0 = Xunorm-repmat(mean(Xunorm),n,1);

% Normalize X to length sqrt(n)
normX = repmat(sqrt(sum(X0.^2)),n,1);

X = sqrt(n)*X0./normX;

% Run the cross-validated LASSO as reference
%tic; [B Fitinfo] = lasso(X,Y,'CV',10,'LambdaRatio',1e-3); solveLASSOCV = toc

% Choose lambda
%lamMSE = Fitinfo.IndexMinMSE;
%betaLASSO_MSE = B(:,lamMSE);

%lam1SE = Fitinfo.Index1SE;
%betaLASSO_1SE = B(:,lam1SE);

% Run the standard q-TREX

% Number of bootstrap runs
L = 0;

% TREX options
% trexopts.boots = L;
% trexopts.rep = p/2;
% trexopts.verbose = 0;
% trexopts.MaxIter = p/2;
% trexopts.progTol = 1e-9;
% trexopts.optTol = 1e-9;
% trexopts.q = qNorm;                % Order of the norm to approximate max norm
% trexopts.normConst = normConst;    % Normalizing constant
% trexopts.beta0 = zeros(p,1);       % Include the zero vector as start point
% 
% tic;[betaTREX,betaTREXMat,funTREXMat] = trexp(X,Y,trexopts); solveQTREX = toc
% 
% % Filter numerical zeros
% epsPssg = 1e-4;
% betaTREX(abs(betaTREX)<epsPssg) = 0;
% betaTREXMat(betaTREXMat(:)<epsPssg) = 0;

clear trexEcosopts;
trexEcosopts.cpath = normConst;          % Constant c
trexEcosopts.verbose = 0;          % No diagnostic output
trexEcosopts.feastol = 1e-9;       % Tolerances in the solver
trexEcosopts.reltol = 1e-9;
trexEcosopts.abstol = 1e-9;

tic;[betaTREXC,betaTREXMatC,funTREXMatC,outEcos] = trex_ecosp(X,Y,trexEcosopts); solveCTREX = toc

% Filter numerical zeros
%epsEcos = 1e-12;
%betaTREXC(abs(betaTREXC)<epsEcos) = 0;
%betaTREXMatC(betaTREXMatC(:)<epsEcos) = 0;