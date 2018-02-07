% Plot figures for illustration of prediction error of LASSO and TREX
% used in
%
% Prediction Error Bounds for Linear Regression With the TREX
% Jacob Bien, Irina Gaynanova, Johannes Lederer, Christian Müller

% Select dimension
p=64;

if p==64
    load('Wainwright_p64_alpha2_norm_kap_0_const2_nReps51')
elseif p==128
    load('Wainwright_p128_alpha2_norm_kap_0_const2_nReps51')
elseif p==256
    load('Wainwright_p256_alpha2_norm_kap_0_const2_nReps51')
end

% Number of repetitions
nReps = 51;

% Flag if uhat should also be plotted
plot_uhatFlag = 0;

% Various matrices to collect
nRegParams = nRegParamsInit;
nCParams = nRegParamsInit;

predMatTREX = zeros(nReps,nCParams);
predBICTREX = zeros(nReps,1);
predMatLasso = zeros(nReps,nRegParams);
predMSELasso = zeros(nReps,1);
pred1SELasso = zeros(nReps,1);
predBICLasso = zeros(nReps,1);

lamMSELasso = zeros(nReps,1);
lam1SELasso = zeros(nReps,1);
lamBICLasso = zeros(nReps,1);


cminVecTREX = zeros(nReps,1);
cmaxVecTREX = zeros(nReps,1);
minBICVecTREX = zeros(nReps,1);

minPredLasso = zeros(nReps,1);
minPredTREX = zeros(nReps,1);

% Ratio of uhat(TREX) for small regularization and |X'eps|_inf
uhatRatio = zeros(nReps,1);

lminLasso = zeros(nReps,1);

xepsVec = zeros(nReps,1);
xyVec = zeros(nReps,1);

%% Computing prediction errror and other statistics across all repetitions
for i=1:nReps
    
    X = XCell{i};
    Y = YCell{i};
    
    
    noiseVec = noiseVecCell{i};
    betaLasso = betaLassoCell{i};
    infoLasso = infoLassoMat{i};
    betaTREX = betaTREXCell{i};
    betaTREXMat = betaTREXMatCell{i};
    
    % Element-wise median over all 2p TREX solutions
    betaTREXM = squeeze(median(betaTREXMat,2));
    
    nRegParams2 = size(betaLasso,2);
    nCParams2 = size(betaTREX,2);
    
    pred_trex = 1/n*sum((repmat(X*betaTrue,1,nCParams2)-X*betaTREX).^2);
    uhat_trex = max(abs(X'*(repmat(Y,1,nCParams2)-X*betaTREX)));
    
    
    [~,minBICT] = min(n*log(sum((repmat(Y,1,size(betaTREX,2))-X*betaTREX).^2,1)/n) + sum(betaTREX~=0).*log(n));
    minBICVecTREX(i) = optsTREXCell{i}.cpath(minBICT);
    predBICTREX(i) = pred_trex(minBICT);
    predMatTREX(i,nCParams-nCParams2+1:nCParams) = pred_trex;
    
    % Prediction error of the median TREX solution (not plotted here)
    pred_trexm = 1/n*sum((repmat(X*betaTrue,1,nCParams2)-X*betaTREXM).^2);
    
    pred_lasso = 1/n*sum((repmat(X*betaTrue,1,nRegParams2)-X*betaLasso).^2);
    
    uhat_lasso = max(abs(X'*(repmat(Y,1,nRegParams2)-X*betaLasso)));
    
    
    [~,minBICL] = min(n*log(sum((repmat(Y,1,size(betaLasso,2))-X*betaLasso).^2,1)/n) + sum(betaLasso~=0).*log(n));
    
    lamBICLasso(i) = infoLasso.Lambda(minBICL);
    lamMSELasso(i) = infoLasso.LambdaMinMSE;
    lam1SELasso(i) = infoLasso.Lambda1SE;
    
    predMatLasso(i,nRegParams-nRegParams2+1:nRegParams) = pred_lasso;
    
    predMSELasso(i) = pred_lasso(infoLasso.IndexMinMSE);
    pred1SELasso(i) = pred_lasso(infoLasso.Index1SE);
    predBICLasso(i) = pred_lasso(minBICL);
    
    cminVecTREX(i) = findTREXLB(X,Y,betaTrue,noiseVec,sig);
    cmaxVecTREX(i) = optsTREXCell{i}.cpath(end);
    
    xepsVec(i) = max(abs(X'.*sig*noiseVec));
    xyVec(i) = max(abs(X'*Y));
    lminLasso(i) = xepsVec(i)/(infoLasso.Lambda(end)*n);
    
    relLassoPath = n*infoLasso.Lambda./(xyVec(i));
    
    rel_cminVecTREX(i) = cminVecTREX(i)./cmaxVecTREX(i);
    
    valIndsL = find(relLassoPath>=lminLasso(i));
    invalIndsL = find(relLassoPath<lminLasso(i));
    
    valIndsT = find(relLassoPath>=rel_cminVecTREX(i));
    invalIndsT = find(relLassoPath<rel_cminVecTREX(i));
    
    hFig = figure(p);
    
    linestyles = {'--';'--';'--';'--'};
    colors = {[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250];[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880]};
    
    alphavalues = [.2 .2 .2 .2];
    
    k=1;
    loglog(relLassoPath(valIndsL),pred_lasso(valIndsL),'Color',colors{k},'LineWidth',2)
    hold on
    k=2;
    loglog(relLassoPath(valIndsT),pred_trex(valIndsT),'Color',colors{k},'LineWidth',2)
    % This would be the median TREX...
    %loglog(relLassoPath(valIndsT),pred_trexm(valIndsT),'Color',colors{k},'LineWidth',2)
    %patchline(relLassoPath,pred_trexm,'linestyle',linestyles{k},'edgecolor',colors{k},'linewidth',1,'edgealpha',alphavalues(k))
    grid on
    plot(n*lamMSELasso(i)./xyVec(i),predMSELasso(i),'.','MarkerSize',40,'Color',colors{3})
    plot(n*lam1SELasso(i)./xyVec(i),pred1SELasso(i),'.','MarkerSize',40,'Color',colors{4})
    plot(n*lamBICLasso(i)./xyVec(i),predBICLasso(i),'.','MarkerSize',40,'Color',colors{5})
    
    if i==1
        title(['p = ',num2str(p)])
        legend(' LASSO path',' TREX path',' LASSO MSE',' LASSO 1SE',' LASSO BIC', 'Location','NW')
    end
    
    k=1;
    patchline(relLassoPath,pred_lasso,'linestyle',linestyles{k},'edgecolor',colors{k},'linewidth',1,'edgealpha',alphavalues(k))
    k=2;
    patchline(relLassoPath,pred_trex,'linestyle',linestyles{k},'edgecolor',colors{k},'linewidth',1,'edgealpha',alphavalues(k))
    
    % Standard scaling for all dimensions
    xlim([2e-3 1])
    ylim([3e-2 2e2])
    
    xlabel('Relative regularization path \rho')
    ylabel('Prediction error')
    
    
    set(gca,'FontSize',30)
    drawnow
    box on
    set(hFig, 'Position', [100 100 1500 1000])
    
    minPredLasso(i) = min(pred_lasso(valIndsL));
    minPredTREX(i) = min(pred_trex(valIndsT));
    
    
    if plot_uhatFlag
        
        hFig2 = figure(p+10);
        
        k=1;
        loglog(relLassoPath(valIndsL),uhat_lasso(valIndsL),'Color',colors{k},'LineWidth',2)
        hold on
        k=2;
        loglog(relLassoPath(valIndsT),uhat_trex(valIndsT),'Color',colors{k},'LineWidth',2)
        
        if i==1
            loglog(relLassoPath,xepsVec(i)*ones(length(relLassoPath),1),'k--','LineWidth',2)
        else
            patchline(relLassoPath,xepsVec(i)*ones(length(relLassoPath),1),'linestyle',linestyles{4},'edgecolor',[0 0 0],'linewidth',1,'edgealpha',alphavalues(4))
        end
        
        grid on
        
        if i==1
            legend(' LASSO path',' TREX path',' \|X^T \eps\|', 'Location','NW')
        end
        
        k=1;
        patchline(relLassoPath,uhat_lasso,'linestyle',linestyles{k},'edgecolor',colors{k},'linewidth',1,'edgealpha',alphavalues(k))
        k=2;
        patchline(relLassoPath,uhat_trex,'linestyle',linestyles{k},'edgecolor',colors{k},'linewidth',1,'edgealpha',alphavalues(k))
        
        xlim([2e-3 1])
        xlabel('Regularization scale \rho')
        ylabel('\hat u')
        
        
        set(gca,'FontSize',30)
        drawnow
        box on
        set(hFig, 'Position', [100 100 1500 1000])
    end
    uhatRatio(i) = uhat_trex(1)/xepsVec(i);
    
end

disp('Average minimally achievable prediction error')
minAvgLassoPredError = mean(minPredLasso)

minAvgTREXPredError = mean(minPredTREX)
