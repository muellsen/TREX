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

qNorm = 40;
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
tic; [B Fitinfo] = lasso(X,Y,'CV',10,'LambdaRatio',1e-3); solveLASSOCV = toc

% Choose lambda
lamMSE = Fitinfo.IndexMinMSE;
betaLASSO_MSE = B(:,lamMSE);

lam1SE = Fitinfo.Index1SE;
betaLASSO_1SE = B(:,lam1SE);

% Run the standard q-TREX

% Number of bootstrap runs
L = 0;

% TREX options
trexopts.boots = L;
trexopts.rep = 11;
trexopts.verbose = 0;
trexopts.q = qNorm;             % Order of the norm to approximate max norm
trexopts.normConst = normConst;    % Normalizing constant
trexopts.beta0 = zeros(p,1); % Include the zero vector as start point
tic;[betaTREX,betaTREXMat,funTREXMat] = trexp(X,Y,trexopts); solveQTREX = toc

% Filter numerical zeros
epsPssg = 1e-4;
betaTREX(abs(betaTREX)<epsPssg) = 0;
betaTREXMat(betaTREXMat(:)<epsPssg) = 0;


trexEcosopts.cpath = 0.5;              % Constant c
trexEcosopts.verbose = 0;          % No diagnostic output
trexEcosopts.feastol = 1e-6;       % Tolerances in the solver
trexEcosopts.reltol = 1e-6;
trexEcosopts.abstol = 1e-6;

tic;[betaTREXC,betaTREXMatC,funTREXMatC] = trex_ecosp(X,Y,trexEcosopts); solveCTREX = toc

% Filter numerical zeros
epsECOS = 1e-6;
betaTREXC(abs(betaTREXC)<epsECOS) = 0;
betaTREXMatC(betaTREXMatC(:)<epsECOS) = 0;

% Support estimation and least-squares refitting
betaInds = (betaTREX~=0);
Xrf = X(:,betaInds);
temp = (Xrf'*Xrf)\Xrf'*Y;
betaRFTREX=zeros(p,1);
betaRFTREX(find(betaTREX~=0)) = temp;

% betaBootInds = (betaBoot~=0);
% Xrf = X(:,betaBootInds);
% temp = (Xrf'*Xrf)\Xrf'*Y;
% betaBootRFTREX=zeros(p,1);
% betaBootRFTREX(find(betaBoot~=0)) = temp;

Y_true = Y;
[Y_sorted, sortedInds] = sort(Y_true);

% Estimated response
Y_LassoMSE = X*betaLASSO_MSE+Fitinfo.Intercept(lamMSE);
Y_Lasso1SE = X*betaLASSO_1SE+Fitinfo.Intercept(lam1SE);
Y_TREX = X*betaTREX;
Y_TREXRF = X*betaRFTREX;

figure;
plot(Y_sorted,'k','LineWidth',5);
hold on
grid on

% Plot the LASSO solutions
plot(Y_Lasso1SE(sortedInds),'LineWidth',3);
plot(Y_LassoMSE(sortedInds),'LineWidth',3);

% Plot the TREX solutions
plot(Y_TREX(sortedInds),'LineWidth',3)
plot(Y_TREXRF(sortedInds),'LineWidth',3)
legend({' Response',' LASSO-CV-MSE',' LASSO-CV-1SE',' TREX',' TREX-RF'},'Location','NW','FontSize',10)
xlabel('Strains','FontSize',14)
ylabel('log(production rate)','FontSize',14)
title('Normalized riboflavin production rate','FontSize',14)
set(gca,'LineWidth',2)

% saveas(gcf,'Response_all','png')

fTrex = objTREX(betaTREX) + sum(abs(betaTREX));
fTrexRF = objTREX(betaRFTREX) + sum(abs(betaRFTREX));

f_MSE = objTREX(betaLASSO_MSE) + sum(abs(betaLASSO_MSE));
f_1SE = objTREX(betaLASSO_1SE) + sum(abs(betaLASSO_1SE));

figure;
% Plot the LASSO solutions
stem(betaLASSO_1SE,'LineWidth',5,'MarkerSize',5);
hold on
grid on
stem(betaLASSO_MSE,'LineWidth',5,'MarkerSize',5);
% Plot the TREX solutions
stem(betaTREX,'LineWidth',5,'MarkerSize',5)
stem(betaRFTREX,'LineWidth',5,'MarkerSize',5)
xlabel('Genes','FontSize',14)
ylabel('\beta values','FontSize',14)
legend({'LASSO-CV-MSE',' LASSO-CV-1SE',' TREX',' TREX-RF'},'Location','NW')
% saveas(gcf,'beta_solutions','png')

% Show all non-zero TREX solutions and corresponding gene names
geneSubsetInds = sum(betaTREXMat,2)>0;

figure;
imagesc(betaTREXMat(geneSubsetInds,:));
set(gca,'YTick',1:length(geneSubsetInds));
set(gca,'YTickLabel',riboTags(geneSubsetInds));

