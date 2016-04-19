% Analysis script for runs generated in testEcosTREX.m or testSCSTREX.m
 
% Assume that all data from previous test script are in the workspace

% Here for completeness the settings from a typical testTREX script
% % Dimension = number of columns of X
% p = 100;
% disp(['Dimension of predictors: ',num2str(p)])
%
% % Noise
% noiseType = {'norm'};
% noiseInd = 1;
%
% % Sample sizes
% n = 50;
% numSamples = n;
% disp(['Number of data points: ',num2str(n)])
%
% % Correlation vector
% kappaVec = [0,0.9];
% numKappa = length(kappaVec);
%
% % Noise level
% sigVec = [0.1,0.5,3];
% numSig = length(sigVec);
%
% % Number of repetitions
% numRep = 21;

% Collect data
meanRunTime = squeeze(mean(runTimeMat,3));
meanEstMat = squeeze(mean(estMat,3));
meanPredMat = squeeze(mean(predMat,3));
meanTpMat = squeeze(mean(tpMat,3));
meanFpMat = squeeze(mean(fpMat,3));
meanFnMat = squeeze(mean(fnMat,3));

allMeanData = zeros(numKappa,numSig,2,6);
allMeanData(:,:,:,1) = meanRunTime;
allMeanData(:,:,:,2) = meanEstMat;
allMeanData(:,:,:,3) = meanPredMat;
allMeanData(:,:,:,4) = meanTpMat;
allMeanData(:,:,:,5) = meanFpMat;
allMeanData(:,:,:,6) = meanFnMat;

stdRunTime = squeeze(std(runTimeMat,[],3));
stdEstMat = squeeze(std(estMat,[],3));
stdPredMat = squeeze(std(predMat,[],3));
stdTpMat = squeeze(std(tpMat,[],3));
stdFpMat = squeeze(std(fpMat,[],3));
stdFnMat = squeeze(std(fnMat,[],3));

allStdData = zeros(numKappa,numSig,2,6);
allStdData(:,:,:,1) = stdRunTime;
allStdData(:,:,:,2) = stdEstMat;
allStdData(:,:,:,3) = stdPredMat;
allStdData(:,:,:,4) = stdTpMat;
allStdData(:,:,:,5) = stdFpMat;
allStdData(:,:,:,6) = stdFnMat;


kapNames ={'0','09'};

for k=1:numKappa
    
    kapNam = kapNames{k};
    
    % RUNTIME
    nameVec = ['Runtime','_kap_',kapNam];
    yVec = ['Run time '];
    
    yMax = max(max(max(allMeanData(k,:,:,1))));
    
    figName = nameVec;
    titleName = [yVec,'for \kappa=',num2str(kappaVec(k))];
    
    currMeanBar = squeeze(allMeanData(k,:,:,1));
    currStdBar = squeeze(allStdData(k,:,:,1));
    
    figure;
    barweb(currMeanBar,currStdBar)
    % Adapt y scale to have enough space for legend
    ylim([0 yMax+0.4*yMax])
    xlabel('\sigma','FontSize',22)
    ylabel([yVec,' /[s]'],'FontSize',22)
    title (titleName,'FontSize',22)
    grid on
    box on
    set(gca,'XTick',[1:numSig])
    set(gca,'XTickLabel',sigVec)
    set(gca,'FontSize',22)
    xlabel('\sigma','FontSize',22)
    set(gca,'FontSize',22)
    legend({' q-TREX',' c-TREX'},'Location','NorthEast', 'FontSize',20)
    
    saveas(gcf,[figName,'_p',num2str(p),'_n',num2str(n)],'png')
end

nameVec = {'TruePositives','FalsePositives','FalseNegatives'};
yVec = {'True positives','False positives','False negatives'};

for k=1:numKappa
    
    kapNam = kapNames{k};
    
    for i=1:length(nameVec)
        
        yMax = max(max(max(allMeanData(k,:,:,i+3))));
        
        figName = [nameVec{i},'_kap_',kapNam];
        titleName = [yVec{i},' for  \kappa=',num2str(kappaVec(k))];
        
        figure;
        barweb(squeeze(allMeanData(k,:,:,i+3)),squeeze(allStdData(k,:,:,i+3)))
        
        % Adapt y scale to have enough space for legend
        ylim([0 yMax+0.4*yMax])
        
        ylabel(yVec{i},'FontSize',22)
        title (titleName,'FontSize',22)
        grid on
        box on
        set(gca,'XTick',[1:numSig])
        set(gca,'XTickLabel',sigVec)
        set(gca,'FontSize',22)
        xlabel('\sigma','FontSize',22)
        set(gca,'FontSize',22)
        legend({' q-TREX',' c-TREX'},'Location','NorthEast', 'FontSize',20)
        
        saveas(gcf,[figName,'_p',num2str(p),'_n',num2str(n)],'png')
        
    end
end

% Hamming distance for q-TREX and c-TREX
supportCell = cell(numSig,numKappa,2);
hammingMat = zeros(numSig,numKappa,2,numRep);

for k=1:numKappa
    
    for s=1:numSig
        
        tempMat = zeros(p,2,numRep);
        
        for r=1:numRep
            % Extract best solution out of all restarts
            [~,minIndQ] = min(funCellQ{k,s,r});
            tempQ = betaCellQ{k,s,r};
            tempMat(:,1,r) = tempQ(:,minIndQ);
            
            % Extract best solution across all 2p solutions
            [~,minIndC] = min(funCellC{k,s,r});
            tempC = betaCellC{k,s,r};
            tempMat(:,2,r) = tempC(:,minIndC);
        end
        
        supportCell{k,s,1} = tempMat(:,1,:)~=0;
        supportCell{k,s,2} = tempMat(:,2,:)~=0;
        
        hammingMat(k,s,1,:) = sum((squeeze(supportCell{k,s,1})-repmat(betaTrue,1,numRep))~=0);
        hammingMat(k,s,2,:) = sum((squeeze(supportCell{k,s,2})-repmat(betaTrue,1,numRep))~=0);
        
        
    end
end

% Plot the Hamming distance

for k=1:numKappa
    
    kapNam = kapNames{k};
    
    figName = ['Hamming_kap_',kapNam];
    
    figure;
    aveHamming = squeeze(mean(hammingMat(k,:,:,:),4));
    stdHamming = squeeze(std(hammingMat(k,:,:,:),[],4));
    
    barweb(1e-6+aveHamming,1e-6+stdHamming)
    grid on
    box on
    set(gca,'XTick',[1:numSig])
    set(gca,'XTickLabel',sigVec)
    set(gca,'FontSize',22)
    xlabel('\sigma','FontSize',22)
    ylabel('Hamming distance','FontSize',22)
    title(['\kappa = ',num2str(kappaVec(k))])
    %xlim([0.05 3.1])
    %ylim([0 50])
    legend({' q-TREX',' c-TREX'},'Location','NorthEast', 'FontSize',20)
    
    saveas(gcf,[figName,'_p',num2str(p),'_n',num2str(n)],'png')
    
end