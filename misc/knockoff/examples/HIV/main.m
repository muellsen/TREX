% This script runs the knockoff analysis for a specified drug class
% and displays a plot showing the number of true and false discoveries.

%% Load the data

% Possible drug types are 'NRTI', 'NNRTI', and 'PI'.
drugClass = 'PI';

[geneData,drugData] = readGeneData(drugClass);
tsmData = readTSMData(drugClass);

%% Run the knockoff analysis

[discoveries, fdp, discoveries_bh, fdp_bh, nData, pData] ...
    = runGeneKnockoff(geneData, drugData, tsmData);

%% Plot the results

drugNames = drugData.Properties.VariableNames;
nDrugs = size(drugData,2);
nTSM = length(tsmData);

plotData = [discoveries .* (1-fdp); discoveries .* fdp]';
plotData_bh = [discoveries_bh .* (1-fdp_bh); discoveries_bh .* fdp_bh]';

h = figure;
set(h, 'Position', [ 440   378   600   40+180*ceil(nDrugs/3)]);
for i=1:nDrugs,
    subplot(ceil(nDrugs/3),3,i)
    R=[plotData(i,:);plotData_bh(i,:)]';
    barplot=bar(R','stacked');
    set(barplot(1),'FaceColor',[.15 0 .5])
    set(barplot(2),'FaceColor',[.8 .3 .1])
    if(i==1)
        legend('In TSM list','Not in TSM list');
        if(strcmp(drugClass,'PI'))
            ylabel('# HIV-1 protease positions selected');
        else
            ylabel('# HIV-1 RT positions selected');
        end
    end
    ymax=max([7+nTSM;1+discoveries(:);1+discoveries_bh(:)]);
    axis([0 3 0 ymax])
    title(['Resistance to ',char(drugNames(i))])
    hold on
    plot(xlim,[nTSM nTSM],'--k')
    hold off
    set(gca,'XTick',1:2,'XTickLabel',{'Knockoff' 'BHq'});
    text(0,-ymax/6,strcat('(Data set size: n=',num2str(nData(i)),...
        ', p=',num2str(pData(i)),')'));
    axis square;
end
