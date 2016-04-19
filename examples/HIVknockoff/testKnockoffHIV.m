% This script runs the knockoff analysis for three specified drug class
% and displays a plot showing the number of true and false discoveries and
% saves the data in a mat file

% Clear workspace
clear all

%% Load the data
% Possible drug types are 'NRTI', 'NNRTI', and 'PI'.
drugClassCell = {'NRTI','PI','NNRTI'};

% Loop over all drug classes
for d = 1:length(drugClassCell)
    
    drugClass = drugClassCell{d}
    
    [geneData,drugData] = readGeneData(drugClass);
    tsmData = readTSMData(drugClass);
    
    %% Run the knockoff analysis

    % Target FDP
    targetFDP = 0.20;

    % TREX (qTREX or cTREX)
    % trexKnockoffStatisticC returns the W^f knockoff statistic
    % trexKnockoffStatistic returns the W^\phi knockoff statistic
	
    statMethod = @trexKnockoffStatisticC;
    tic
    [discoveriesTREX, fdpTREX, discoveriesTREX_bh, fdpTREX_bh, nDataTREX, pDataTREX] = runGeneKnockoffMod(geneData, drugData, tsmData,statMethod,targetFDP);
    eval(['timeTREX_',drugClass,' = toc;']);
    
    % Original call
    %[discoveriesTREX, fdpTREX, discoveries_bh, fdp_bh, nData, pData] ...
    %    = runGeneKnockoff(geneData, drugData, tsmData,statMethod);
    
    % LASSO
    statMethod = @knockoff.stats.lassoSignedMax;
    
    tic
    [discoveriesLASSO, fdpLASSO, discoveries_bh, fdp_bh, nData, pData] ...
        = runGeneKnockoffMod(geneData, drugData, tsmData,statMethod,targetFDP);
    eval(['timeLASSO_',drugClass,' = toc;']);
    
    %% Plot the results
    
    drugNames = drugData.Properties.VariableNames;
    nDrugs = size(drugData,2);
    nTSM = length(tsmData);
    
    plotDataTREX = [discoveriesTREX .* (1-fdpTREX); discoveriesTREX .* fdpTREX]';
    plotDataLASSO = [discoveriesLASSO .* (1-fdpLASSO); discoveriesLASSO .* fdpLASSO]';
    
    plotData_bh = [discoveries_bh .* (1-fdp_bh); discoveries_bh .* fdp_bh]';
    
    h = figure;
    set(h, 'Position', [ 440   378   600   40+180*ceil(nDrugs/3)]);
    for i=1:nDrugs,
        subplot(ceil(nDrugs/3),3,i)
        R=[plotDataTREX(i,:);plotDataLASSO(i,:);plotData_bh(i,:)]';
        barplot=bar(R','stacked');
        set(barplot(1),'FaceColor',[.15 0 .5])
        set(barplot(2),'FaceColor',[.8 .3 .1])
        
        if(i==1)
            if(strcmp(drugClass,'PI'))
                ylabel('# HIV-1 PROT positions selected');
                legend({'In TSM list','Not in TSM list'},'FontSize',5);
                
            else
                ylabel('# HIV-1 RT positions selected');
                legend('In TSM list','Not in TSM list');
                
            end
        end
        ymax=max([7+nTSM;1+discoveriesTREX(:);1+discoveries_bh(:)]);
        %axis([0 3 0 ymax])
        title(['Resistance to ',char(drugNames(i))])
        hold on
        plot(xlim,[nTSM nTSM],'--k')
        hold off
        set(gca,'XTick',1:3,'XTickLabel',{'TREX' 'LASSO' 'BHq'});
        text(0,-ymax/6,strcat('(Data set size: n=',num2str(nData(i)),...
            ', p=',num2str(pData(i)),')'));
        axis square;
        if(strcmp(drugClass,'PI'))
            set(gca,'FontSize',7)
        else
            set(gca,'FontSize',9)
        end
        
    end
    
    %% Save the results
    saveas(gcf,['SummaryStats_',drugClass],'png')
    saveas(gcf,['SummaryStats_',drugClass],'fig')

    save(['workspace_',date,'_',drugClass])
end
