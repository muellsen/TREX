% Analysis of TREX solutions on the Riboflavin (vitamin B) production test case from
% Buehlmann et al. 2014
% The TREX files are available at https://tinyurl.com/kcz5jz3
% Download the .mat files and save them in this folder

% clear all
close all

global X
global Y
global normConst
global qNorm

% Load data
% riboData
% riboTags

if ~exist('riboData')    
    load riboflavinData.mat
end

% Response
Yunorm = riboData(1,:)';
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

geneNames = riboTags;
for i=1:length(riboTags)
   temp = strsplit(riboTags{i},'_');
   geneNames{i} = temp{1};
end


% Load cTREX (SCS/ECOS) and qTREX solutions on riboflavin data
try
    
load riboTestCaseData_scsOnly.mat
load riboTestCaseData_ecosOnly.mat
load riboTestCaseData_qTREXOnly.mat

catch

    error('Download mat files data from https://tinyurl.com/kcz5jz3 and store in directory')

end


% Analyzing the non-convex optimization landscape as explored by c-TREX (with SCS) and
% q-TREX (with PSG) 

% Number of random restarts
nRepPSG = size(betaTREXMatQ,2);

% Number of 2p local solutions
nRepSCS = 2*p;

% Number of 2p local solutions
nRepECOS = 2*p;


% epsSCS from loaded data
epsSCS = 1e-5;
betaTREXSCS(abs(betaTREXSCS)<epsSCS) = 0;
betaTREXMatSCS(abs(betaTREXMatSCS(:))<epsSCS) = 0;

% epsEcos from loaded data
epsEcos = 1e-10;
%epsEcos = 1e-6;
betaTREXEcos(abs(betaTREXEcos)<epsEcos) = 0;
betaTREXMatEcos(abs(betaTREXMatEcos(:))<epsEcos) = 0;

% epsPsg set here
epsPsg = 1e-10;
%epsPsg = 1e-6;
betaTREXQ(abs(betaTREXQ)<epsPsg) = 0;
betaTREXMatQ(abs(betaTREXMatQ(:))<epsPsg) = 0;

% Recompute the acutal function value of the thresholded solutions

normConst = trexopts.normConst;
qNorm = trexopts.q;

funTREXMatTh = funTREXMatQ;

for i=1:nRepPSG
    
    temp = objTREX_SA(betaTREXMatQ(:,i),X,Y,normConst);
    funTREXMatTh(i) = temp;
    %temp = objTREX_wog(betaTREXMat(:,i))
    %funTREXMatTh(i) = temp;
    
end

[~,minQInd] = min(funTREXMatTh);


funTREXMatEcosTh = funTREXMatEcos;
for i=1:nRepECOS
    
    funTREXMatEcosTh(i) = objTREX_SA(betaTREXMatEcos(:,i),X,Y,normConst);
    
end

[~,minEcosInd] = min(funTREXMatEcosTh);

% Redefine new beta solution based on true ranking
betaTREXEcosOrig = betaTREXEcos;
betaTREXQOrig = betaTREXQ;

betaTREXEcos = betaTREXMatEcos(:,minEcosInd);
betaTREXQ = betaTREXMatQ(:,minQInd);


% Plot thresholded function value vs. original function value
figure;
plot(funTREXMatQ,funTREXMatTh,'.','MarkerSize',30);
grid on
set(gca,'FontSize',20)
ylabel('Exact TREX function value')
xlabel('qNorm TREX Function value')
title('Exact vs. qNorm function value')
hold on
hline = refline(1,0);
hline.LineStyle=('--');

figure;
plot(funTREXMatEcos,funTREXMatEcosTh,'.','MarkerSize',30);
grid on
set(gca,'FontSize',20)
ylabel('Exact TREX function value')
xlabel('cTREX Function value')
title('Exact vs. cTREX function value')
hold on
hline = refline(1,0);
hline.LineStyle=('--');

% Plot thresholded function value histogram

figure;
histogram(funTREXMatTh,'Normalization','probability','BinLimits',[2 4]);
hold on;
histogram(funTREXMatEcosTh,'Normalization','probability','BinLimits',[2 4]);
title('Thresholded exact TREX function value distribution')
legend('qTREX (PSG)','cTREX(Ecos)')
grid on
set(gca,'FontSize',20)

% Plot original function value histogram
figure;
histogram(funTREXMatQ,'Normalization','probability','BinLimits',[2 4]);
hold on;
histogram(funTREXMatEcos,'Normalization','probability','BinLimits',[2 4]);
title('Numerical TREX function value distributions')
legend('qTREX (PSG)','cTREX(Ecos)')
set(gca,'FontSize',20)
grid on


% qTREX solution vectors sorted by function value

[sortedFunQTREX,sortedQTREXInds] = sort(funTREXMatTh);

figure;
imagesc(betaTREXMatQ(:,sortedQTREXInds)~=0);
colorbar
grid on
set(gca,'FontSize',20)
title('qTREX solutions sorted by function value')


% cTREX solution vectors sorted by function value

[sortedFunTREXEcos,sortedTREXEcosInds] = sort(funTREXMatEcosTh);

figure;
imagesc(betaTREXMatSCS(:,sortedTREXEcosInds)~=0);
colorbar
set(gca,'FontSize',20)
grid on
title('cTREX solutions sorted by function value')


% Sorted function value trace
figure;
loglog(sortedFunQTREX,'LineWidth',10)
hold on
loglog(sortedFunTREXEcos,'LineWidth',10)
title('Exact TREX function values')
legend(' qTREX (PSG)',' cTREX(Ecos)')
set(gca,'FontSize',20)
grid on

% Sparsity of the solutions sorted by function value
figure;
loglog(sum(betaTREXMatQ(:,sortedQTREXInds)~=0));
hold on;
loglog(sum(betaTREXMatEcos(:,sortedTREXEcosInds)~=0));
grid on
legend(' qTREX (PSG)',' cTREX(Ecos)')
set(gca,'FontSize',20)
ylabel('Sparsity of solution')
title('Sparsity of solutions sorted by function value')
ylim([0 100])

% Sparsity of the solutions vs  function value
figure;
plot(sortedFunQTREX,sum(betaTREXMatQ(:,sortedQTREXInds)~=0),'.','MarkerSize',30);
grid on
set(gca,'FontSize',20)
ylabel('Sparsity of solution')
xlabel('TREX Function value')
title('Sparsity vs. function value of qTREX')
ylim([0 100])
xlim([2.2 3])

figure;
plot(sortedFunTREXEcos,sum(betaTREXMatEcos(:,sortedTREXEcosInds)~=0),'.','MarkerSize',30);
grid on
set(gca,'FontSize',20)
ylabel('Sparsity of solution')
xlabel('TREX Function value')
title('Sparsity vs. function value of cTREX')
ylim([0 100])
xlim([2.2 3])


% Hitting probability for q-TREX

epsPVec = [0.0001:0.0001:5];
probVec = zeros(length(epsPVec),1); 

mincTREX = min(funTREXMatEcosTh)

i=1;
for epsP=epsPVec
    probVec(i) = sum(funTREXMatTh<mincTREX+epsP)./p;
    i=i+1; 
end

figure;
plot(epsPVec,probVec,'LineWidth',5)
grid on
box on
xlim([0 0.1])
xlabel('Difference to P^*','FontSize',30)
ylabel('Hitting probabilty','FontSize',30)
set(gca,'FontSize',30)


% TREX indices

trexSCSInds = find(betaTREXSCS~=0)
trexEcosInds = find(betaTREXEcos~=0)
trexQInds = find(betaTREXQ~=0)


% Top K indices according to function value
topK = 30;%length(trexEcosInds);
topKInds = round(sortedTREXEcosInds(1:topK)/2)
[sortedTopKInds,~] = sort(topKInds);

% Gene names selected
geneNamesTREXTopK = geneNames(sortedTopKInds)
geneNamesTREX = geneNames(sort(trexEcosInds))
geneNamesTREXQ = geneNames(sort(trexQInds))


% Refitting and plot of the prediction results

Xrf = X(:,topKInds);
temp = (Xrf'*Xrf)\Xrf'*Y;
betaTopTREXRF=zeros(p,1);
betaTopTREXRF(topKInds) = temp;

Xrf = X(:,trexEcosInds);
temp = (Xrf'*Xrf)\Xrf'*Y;
betaTREXRF=zeros(p,1);
betaTREXRF(trexEcosInds) = temp;

Xrf = X(:,trexQInds);
temp = (Xrf'*Xrf)\Xrf'*Y;
betaTREXQRF=zeros(p,1);
betaTREXQRF(trexQInds) = temp;


% Visualization of best TREX solutions
figure;
stem(betaTREXRF,'LineWidth',3);
hold on
stem(betaTREXEcos,'--','LineWidth',7);
grid on
stem(betaTREXQRF,'--','LineWidth',5);
stem(betaTREXQ,'LineStyle','--','LineWidth',5);

legend(' cTREX (Ecos)',' cTREX(RF)',' qTREX(PSG)',' qTREX(RF)')
set(gca,'FontSize',20)
ylabel(' beta hat value')
title('Solutions found by cTREX(Ecos) and qTREX')
xlabel('Variable index j')
xlim([0 p+1])

Y_true = Y;
[Y_sorted, sortedInds] = sort(Y_true);

% Estimated response qTREX
Y_TREXQ = X*betaTREXQ;
Y_TREXQRF = X*betaTREXQRF;

% Estimated response cTREX
Y_TREX = X*betaTREXEcos;
Y_TREXRF = X*betaTREXRF;

Y_TOPTREXRF = X*betaTopTREXRF;

figure;
plot(Y_sorted,'k','LineWidth',5);
hold on
grid on

% Plot the qTREX solutions
plot(Y_TREXQ(sortedInds),'LineWidth',2);
plot(Y_TREXQRF(sortedInds),'LineWidth',5)

% Plot the cTREX solutions
plot(Y_TREX(sortedInds),'LineWidth',2)
plot(Y_TREXRF(sortedInds),'LineWidth',5)
legend({' Response',' cTREX',' cTREX-RF',' qTREX',' qTREX-RF'},'Location','NW','FontSize',10)
xlabel('Strains','FontSize',14)
ylabel('log(production rate)','FontSize',14)
title('Normalized riboflavin production rate','FontSize',14)
set(gca,'LineWidth',2)
set(gca,'FontSize',20)


% Prediction error for all methods
PredTrexQ = sum((Y_true-Y_TREXQ).^2)
PredTrexQRF = sum((Y_true-Y_TREXQRF).^2)

PredTrex = sum((Y_true-Y_TREX).^2)
PredTrexRF = sum((Y_true-Y_TREXRF).^2)
PredTopTrexRF = sum((Y_true-Y_TOPTREXRF).^2)

AICTrexQRF = n + n*log(pi)+n*log(PredTrexQRF/n) + 2*(nnz(betaTREXQRF)+1)
corrTREXQRF = (corr(Y_true,Y_TREXQRF)).^2
MAETrexQRF = mean(abs(Y_true-Y_TREXQRF))

AICTrexRF = n + n*log(pi)+n*log(PredTrexRF/n) + 2*(nnz(betaTREXRF)+1)
corrTREXRF = (corr(Y_true,Y_TREXRF)).^2
MAETrexRF = mean(abs(Y_true-Y_TREXRF))


% Correlation matrices of the different subsets
figure; 
imagesc(abs(cov(X(:,sortedTopKInds)))>0.5)
set(gca,'YTick',1:1:length(topKInds))
set(gca,'YTickLabel',geneNamesTREXTopK)
set(gca,'FontSize',15)
colorbar
title(['Correlation matrix among the top',num2str(length(topKInds)),' variables (by function value)'])

figure; 
imagesc(abs(cov(X(:,sort(trexEcosInds))))>0.5)
set(gca,'YTick',1:1:length(trexEcosInds))
set(gca,'YTickLabel',geneNamesTREX)
set(gca,'FontSize',15)
colorbar
title(['Correlation matrix among the ', num2str(length(trexEcosInds)),' cTREX variables'])

figure; 
imagesc(abs(cov(X(:,sort(trexQInds))))>0.5)
set(gca,'YTick',1:1:length(trexQInds))
set(gca,'YTickLabel',geneNamesTREXQ)
set(gca,'FontSize',15)
colorbar
title(['Correlation matrix among the ', num2str(length(trexQInds)),' qTREX variables'])


% figure; 
% imagesc(inv(cov(X(:,sortedTopKInds))))
% set(gca,'YTick',1:1:length(topKInds))
% set(gca,'YTickLabel',topGeneNamesTREX)
% set(gca,'FontSize',15)
% colorbar
% title('Inverse correlation matrix among the top 25 variables (by function value)')
% 
% figure; 
% imagesc(abs(inv(cov(X(:,sortedTopKInds))))>20)
% set(gca,'YTick',1:1:length(topKInds))
% set(gca,'YTickLabel',topGeneNamesTREX)
% set(gca,'FontSize',15)
% colorbar
% title('Graph induced by inverse matrix among the top 25 variables (by function value)')








