% The script generates Figure 5 from the
% manuscript based on the information provided by
% smallnTREXsimulations.m and
% largenTREXsimulations.m. 


% Large n setting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('LargeNExample_09-Apr-2016_p100_n500_nRep21.mat')

kapNames ={'0','03','06','09'};
sigNames={'01','05','3'};

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
            title(['n = ',num2str(n),', p = ',num2str(p),', \kappa = ',num2str(kappaVec(k)),', \sigma = ',num2str(sigVec(s))])
            xlabel('Function value','FontSize',18)
    
            saveas(gcf,[figName,'_p',num2str(p),'_n',num2str(n),'_sig',sigNames{s}],'png')
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Small n setting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('SmallNExample_26-Jan-2016_p100 _500_n50_nRep21.mat')

kapNames ={'0','03','06','09'};
sigNames={'01','05','3'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Histogram of 2p values

% number of restarts for q-TREX
numP = 2;
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
            title(['n = ',num2str(n),', p = ',num2str(p),', \kappa = ',num2str(kappaVec(k)),', \sigma = ',num2str(sigVec(s))])
            xlabel('Function value','FontSize',18)
    
            saveas(gcf,[figName,'_p',num2str(p),'_n',num2str(n),'_sig',sigNames{s}],'png')
        end
    end
end
