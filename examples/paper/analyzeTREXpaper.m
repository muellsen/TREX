% Analysis script for runs generated in smallnTREXsimulations.m and
% largenTREXsimulations.m

% Assume that all data from one of these test scripts are in the workspace

% Here for completeness the settings from the n<p case script
% % Dimension = number of columns of X
% p = [100,500];
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
% kappaVec = [0,0.3,0.6,0.9];
% numKappa = length(kappaVec);
%
% % Noise level
% sigVec = [0.1,0.5,3];
% numSig = length(sigVec);
%
% % Number of repetitions
% numRep = 21;

% Collect data (4 is the position of numRep)
meanRunTime = squeeze(mean(runTimeMat,4));
meanEstMat = squeeze(mean(estMat,4));
meanPredMat = squeeze(mean(predMat,4));
meanTpMat = squeeze(mean(tpMat,4));
meanFpMat = squeeze(mean(fpMat,4));
meanFnMat = squeeze(mean(fnMat,4));

allMeanData = zeros(numP,numKappa,numSig,2,6);
allMeanData(:,:,:,:,1) = meanRunTime;
allMeanData(:,:,:,:,2) = meanEstMat;
allMeanData(:,:,:,:,3) = meanPredMat;
allMeanData(:,:,:,:,4) = meanTpMat;
allMeanData(:,:,:,:,5) = meanFpMat;
allMeanData(:,:,:,:,6) = meanFnMat;

stdRunTime = squeeze(std(runTimeMat,[],4));
stdEstMat = squeeze(std(estMat,[],4));
stdPredMat = squeeze(std(predMat,[],4));
stdTpMat = squeeze(std(tpMat,[],4));
stdFpMat = squeeze(std(fpMat,[],4));
stdFnMat = squeeze(std(fnMat,[],4));

allStdData = zeros(numP,numKappa,numSig,2,6);
allStdData(:,:,:,:,1) = stdRunTime;
allStdData(:,:,:,:,2) = stdEstMat;
allStdData(:,:,:,:,3) = stdPredMat;
allStdData(:,:,:,:,4) = stdTpMat;
allStdData(:,:,:,:,5) = stdFpMat;
allStdData(:,:,:,:,6) = stdFnMat;


kapNames ={'0','03','06','09'};
sigNames={'01','05','3'};

for pind=1:numP

    p=pVec(pind);

for k=1:numKappa
    
    kapNam = kapNames{k};
    
    % RUNTIME
    nameVec = ['Runtime','_kap_',kapNam];
    yVec = ['Run time '];
    
    yMax = max(max(max(allMeanData(pind,k,:,:,1))));
    
    figName = nameVec;
    titleName = [yVec,'for p=', num2str(p),' and \kappa=',num2str(kappaVec(k))];
    currMeanBar = squeeze(allMeanData(pind,k,:,:,1));
    currStdBar = squeeze(allStdData(pind,k,:,:,1));
    
    figure;
    barweb(currMeanBar,currStdBar)
    % Adapt y scale to have enough space for legend
    ylim([0 yMax+0.4*yMax])
    xlabel('\sigma','FontSize',22)
    ylabel([yVec,'[s]'],'FontSize',22)
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

%%%%%%%%%%%%%%%%%%%%%%%
% Estimation error
for pind=1:numP

    p=pVec(pind);

for k=1:numKappa
    
    kapNam = kapNames{k};
    
    % Estimation
    nameVec = ['Estimation','_kap_',kapNam];
    yVec = ['Estimation error '];
    
    yMax = max(max(max(allMeanData(pind,k,:,:,2))));
    
    figName = nameVec;
    titleName = [yVec,'for p=', num2str(p),' and \kappa=',num2str(kappaVec(k))];
    currMeanBar = squeeze(allMeanData(pind,k,:,:,2));
    currStdBar = squeeze(allStdData(pind,k,:,:,2));
    
    figure;
    barweb(currMeanBar,currStdBar)
    % Adapt y scale to have enough space for legend
    ylim([0 yMax+0.4*yMax])
    xlabel('\sigma','FontSize',22)
    ylabel(yVec,'FontSize',22)
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


%%%%%%%%% Computational Comparison
%%% funCellQ for beta=0 and funCellC

eps=1e-4; %tolerance for the difference in function values
agreeP=zeros(numP,numKappa,numSig);
objQmat=zeros(numP,numKappa,numSig,numRep); %objective values from q-trex
objCmat=objQmat;% objective values from c-trex


%%%%%%%%%%%%%%%%%%%%%
% only zero vector
for pind=1:numP

p=pVec(pind);
for k=1:numKappa
    
    for s=1:numSig
        
        tempMat = zeros(p,2,numRep);
        
        for r=1:numRep
            % Extract solution corresponding to zero initial point
            objQ = funCellQ{pind,k,s,r}(1);
            objQmat(pind,k,s,r)=objQ;
            
            % Extract best solution across all 2p solutions
            [objC,~] = min(funCellC{pind,k,s,r});
            objCmat(pind,k,s,r)=objC;
            
            % Compare objectives
            if (abs(objQ-objC)<eps)
                agreeP(pind,k,s)=agreeP(pind,k,s)+1;
            end         
        end        
    end
end
end


% plot # of times out of numRep that the difference in objectives is less
% than eps

for pind=1:numP
p=pVec(pind);
for k=1:numKappa
    
    kapNam = kapNames{k};
    figName = ['ObjectiveDiff_kap_',kapNam];
    
    figure;
    
    bar(squeeze(agreeP(pind,k,:))/numRep)
    grid on
    box on
    set(gca,'XTick',[1:numSig])
    set(gca,'XTickLabel',sigVec)
    set(gca,'FontSize',22)
    xlabel('\sigma','FontSize',22)
    ylabel('Fraction of agreement','FontSize',22)
    title(['p = ',num2str(p),' and \kappa = ',num2str(kappaVec(k))])
    ylim([0 1.05])
    
    saveas(gcf,[figName,'_p',num2str(p),'_n',num2str(n)],'png')
end
end

%%% funCellQ for beta[1:nTREXrep] and funCellC
% only random restarts

eps=1e-4;
nRand=nTREXRep-(nl-1);
agreePnrep=zeros(numP,numKappa,numSig,nRand);


for pind=1:numP

p=pVec(pind);
for k=1:numKappa
    
    for s=1:numSig
        
        %tempMat = zeros(p,2,numRep);
        
        for r=1:numRep
            % Extract best solution across all 2p solutions
            [objC,~] = min(funCellC{pind,k,s,r});
            
            % Extract best solution from first j non lasso restarts
            for j=1:nRand
                if j==1
                    % Extract best q-tres solution from first j restarts
                    [objQ,~] = min(funCellQ{pind,k,s,r}(1));
                else
                    [objQ,~] = min(funCellQ{pind,k,s,r}([1,(nl+1):(nl+(j-1))]));
                end

                % Compare objectives
                if (abs(objQ-objC)<eps)
                    agreePnrep(pind,k,s,j)=agreePnrep(pind,k,s,j)+1;
                end  
            end
        end        
    end
end
end



% plot # of times out of numRep that the difference in objectives is less
% than eps
lstyle={'-','--',':'}; % different line styles for different sigmas

for pind=1:numP
p=pVec(pind);
for k=1:numKappa
    
    kapNam = kapNames{k};
    figName = ['ObjectiveDiff_NrepRand_kap_',kapNam];

    figure;
    ylim([-0.05 1.05])
    xlim([1 nRand])
        for s=1:numSig
            hold on;
            plot(1:nRand,squeeze(agreePnrep(pind,k,s,:))/numRep,'LineWidth',5,'LineStyle',lstyle{s});
        end
    grid on
    box on
    legend({'\sigma=0.1','\sigma=0.5','\sigma=3'},'Location','SouthEast', 'FontSize',20)
    set(gca,'FontSize',22)
    xlabel('Number of random restarts','FontSize',22)
    ylabel('Fraction of agreement','FontSize',22)
    title(['n = ',num2str(n),', p = ',num2str(p),' and \kappa = ',num2str(kappaVec(k))])
    
    saveas(gcf,[figName,'_p',num2str(p),'_n',num2str(n)],'png')
end
end

%%% funCellQ for beta[1:nTREXrep] and funCellC
% only LASSO staring points, the results are similar and are not used in the paper

eps=1e-4;
agreePnrepLasso=zeros(numP,numKappa,numSig,nl);


for pind=1:numP

p=pVec(pind);
for k=1:numKappa
    
    for s=1:numSig
        
        
        for r=1:numRep
            % Extract best solution across all 2p solutions
            [objC,~] = min(funCellC{pind,k,s,r});
            
            % Extract best solution from first j non lasso restarts
            for j=1:nl
                % Extract best q-tres solution from first j restarts
                [objQ,~] = min(funCellQ{pind,k,s,r}(1:j));
                    
                % Compare objectives
                if (abs(objQ-objC)<eps)
                    agreePnrepLasso(pind,k,s,j)=agreePnrepLasso(pind,k,s,j)+1;
                end  
            end
        end        
    end
end
end



% plot # of times out of numRep that the difference in objectives is less
% than eps

for pind=1:numP
p=pVec(pind);
for k=1:numKappa
    
    kapNam = kapNames{k};
    figName = ['ObjectiveDiff_NrepLasso_kap_',kapNam];

    figure;
    ylim([-0.05 1.05])
    xlim([1 nl])
        for s=1:numSig
            hold on;

            plot(1:nl,squeeze(agreePnrepLasso(pind,k,s,:))/numRep);
        end
    grid on
    box on
    legend({'\sigma=0.1','\sigma=0.3','\sigma=3'},'Location','SouthEast', 'FontSize',20)
    set(gca,'FontSize',22)
    xlabel('Number of random restarts','FontSize',22)
    ylabel('Fraction of agreement','FontSize',22)
    title(['p = ',num2str(p),' and \kappa = ',num2str(kappaVec(k))])
    
    saveas(gcf,[figName,'_p',num2str(p),'_n',num2str(n)],'png')
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Histogram of 2p values

% number of restarts for q-TREX
nl=20;
nRand=nTREXRep-(nl-1);
for pind=1:numP

    p=pVec(pind);
    for k=1:numKappa
    
        for s=1:numSig
        
            % Extract 2p solutions from c-TREX
            r=1; % fix 1 model instance
            obj2p = squeeze(funCellC{pind,k,s,r});
            
            % Extract solutions to random starts q-TREX
            objQnRand=squeeze(funCellQ{pind,k,s,r}([1,(nl+1):end]));
            % Only take those that are within the borders of obj2p
            objQnRand=objQnRand(objQnRand<max(obj2p)+1);
            nRand=length(objQnRand);
            % Create histogram 
            [sortRand,ind]=sort(objQnRand);
            yaxisval=0.5*ones(nRand,1);
            epsv=1e-3;
            for l=2:nRand
                if abs(sortRand(l)-sortRand(l-1))<epsv
                    sortRand(l-1)=sortRand(l);
                    yaxisval(l)=yaxisval(l-1)+0.5;
                end
            end
                  
           
            figName = ['2pvaluesHistQhistCapped_kap_',kapNames{k}];
            figure;
            histogram(obj2p,20)
            % Put q-trex solutions as red dots to see where they go
            hold on
            plot(sortRand,yaxisval,'dr','MarkerSize',10)
            drawnow
            set(gca,'FontSize',18)
            title(['p = ',num2str(p),', \kappa = ',num2str(kappaVec(k)),', \sigma = ',num2str(sigVec(s))])
            xlabel('Function value','FontSize',18)
    
            saveas(gcf,[figName,'_p',num2str(p),'_n',num2str(n),'_sig',sigNames{s}],'png')
        end
    end
end

