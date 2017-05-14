% Generate Figure 6 (Riboflavin data analysis with TREX) in the manuscipt
% as well as a number of additional Figures (including Figures in the Appendix)
% Data needs to be loaded and pre-processed in analyzeTREXRiboData

% Plot thresholded function value histogram
% Magnified by hand in the viewer for creating inset
figure;
histogram(funTREXMatTh,'Normalization','probability','BinLimits',[2.15 4],'NumBins',80);
hold on;
histogram(funTREXMatEcosTh,'Normalization','probability','BinLimits',[2.15 4],'NumBins',80);
title('Thresholded exact TREX function value distribution')
legend({'q-TREX','c-TREX(ECOS)'},'FontSize',30)
grid on
set(gca,'FontSize',30)

% Sparsity of the solutions vs  function value
figure;
plot(sortedFunQTREX,sum(betaTREXMatQ(:,sortedQTREXInds)~=0),'.','MarkerSize',60);
grid on
set(gca,'FontSize',40)
ylabel('Sparsity of solution')
xlabel('TREX Function value')
title('Sparsity vs. function value of qTREX')
ylim([0 100])
xlim([2.2 3])
hold on
plot(sortedFunTREXEcos,sum(betaTREXMatEcos(:,sortedTREXEcosInds)~=0),'.','MarkerSize',60);
grid on
set(gca,'FontSize',40)
ylabel('Sparsity of solution')
xlabel('TREX Function value')
title('Sparsity vs. function value of q-TREX (blue) and c-TREX(red)')
ylim([10 100])
xlim([2.2 3])
set(gca,'FontSize',40)
hold on
plot(sortedFunQTREX(1),sum(betaTREXMatQ(:,sortedQTREXInds(1))~=0),'b.','MarkerSize',120);
plot(sortedFunTREXEcos(1),sum(betaTREXMatEcos(:,sortedTREXEcosInds(1))~=0),'r.','MarkerSize',120);

% Sparsity of the solutions sorted by function value
figure;
loglog(sum(betaTREXMatQ(:,sortedQTREXInds)~=0),'LineWidth',5);
hold on;
loglog(sum(betaTREXMatEcos(:,sortedTREXEcosInds)~=0),'LineWidth',5);
grid on
legend({' q-TREX    ',' c-TREX(ECOS)'},'FontSize',30)
set(gca,'FontSize',30)
ylabel('Sparsity of solution')
title('Sparsity of solutions sorted by function value')
ylim([0 100])

% Correlation among selected predictors for c-TREX and q-TREX
figure; 
imagesc(corr(X(:,sort(trexEcosInds))))
set(gca,'YTick',1:1:length(trexEcosInds))
set(gca,'YTickLabel',geneNamesTREX)
set(gca,'FontSize',15)
colorbar
title(['Correlation matrix among the ', num2str(length(trexEcosInds)),' c-TREX variables'])

figure; 
imagesc(corr(X(:,sort(trexQInds))))
set(gca,'YTick',1:1:length(trexQInds))
set(gca,'YTickLabel',geneNamesTREXQ)
set(gca,'FontSize',15)
colorbar
title(['Correlation matrix among the ', num2str(length(trexQInds)),' q-TREX variables'])


% Plot the TREX refitted solutions vs measured log-production rates
figure;
plot(Y_sorted,Y_TREXQRF(sortedInds),'.','MarkerSize',70)
hold on
grid on
plot(Y_sorted,Y_TREXRF(sortedInds),'.','MarkerSize',70)
legend({' q-TREX (LS-RF)',' c-TREX (LS-RF)'},'Location','NW','FontSize',40)
line([-3 2],[-3 2],'LineStyle','--','Color','k','LineWidth',2)
ylabel('Fitted log production rates','FontSize',40)
xlabel('Measured log production rate','FontSize',40)
title('Riboflavin data','FontSize',40)
set(gca,'LineWidth',2)
set(gca,'FontSize',40)
xlim([-3 2])
ylim([-3 2])



