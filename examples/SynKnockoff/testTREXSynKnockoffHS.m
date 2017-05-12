% This script repeats and extends the experiment from Barber-Foygel and Candes page 11

% Script that tests target FDR vs the average FDP and corresponding power
% when the knockoff filter (TREX, LASSO) is applied to a synthetic data set.

rng(12345) % Random seed

%% Synthetic problem parameters

%pVec = [50,100,200,500];
pVec = 100;

% Number of variables
for p=pVec
    p
    nVec = p+1;%[p+1,p+11,p+26,p+51];
    
    k = 30; % Number of variables with nonzero coefficientss
    
    rhoVal = 0.3;
    covMat = rhoVal*ones(p,p);
    covMat(1:p+1:p^2) = 1;
    cholMat = chol(covMat);
    
    % Number of data points
    for n=nVec
        n
        
        %% Introduce heteroscedastic noise
        for sigma1 = 1;%[0.5, 0.7, 1]
            
            sigma2 = sqrt(2-sigma1^2);
            
            for amplitude = 3.5;%[1.5,2.5,3.5];  % Magnitude of nonzero coefficients
                
                if amplitude==1.5 && sigma1==0.5 && n==101
                    continue
                end
                
                
                %% Synthetic problem construction
                
                X = (cholMat*randn(p,n))';
                % Center and normalize the data
                X = X-repmat(mean(X),n,1);
                X = X./repmat(std(X,[],1),n,1);
                
                S0 = 1:k;
                beta = zeros(p,1);
                beta(S0) = amplitude*[(-1).^(S0)]';
                
                sampleY = @() X*beta + randSigmaHS(n,sigma1,sigma2) .* randn(n,1);
                
                %% Save achieved FDR and #TP
                
                FDP = @(S) sum(beta(S) == 0) / max(1, length(S));
                TPP = @(S) sum(beta(S) ~= 0) ;%/ length(S);
                
                ntrials = 51;
                qVec = [0.05:0.05:0.95];
                nQ = length(qVec);
                
                fdp_TREXq = zeros(nQ, ntrials);
                fdp_TREXf = zeros(nQ, ntrials);
                
                fdp_LASSO = zeros(nQ, ntrials);
                
                tpp_TREXq = zeros(nQ, ntrials);
                tpp_TREXf = zeros(nQ, ntrials);
                tpp_LASSO = zeros(nQ, ntrials);
                
                for j = 1:ntrials
                    
                    y = sampleY();
                    
                    [~,W_TREXq] = knockoff.filter(X, y, qVec(1),'Statistic',@trexSignedMax);
                    %X_ko = knockoff.create(X);
                    %[W_TREX,Z_TREX,trexTrace,trexopts] = trexSignedMax(X,X_ko,y)
                    %t = knockoff.threshold(W_TREX, q);
                    
                    [~,W_TREXf] = knockoff.filter(X, y, qVec(1),'Statistic',@trexKnockoffStatisticC);
                    
                    [~,W_LASSO] = knockoff.filter(X, y, qVec(1));
                    
                    
                    for i=1:length(qVec)
                        
                        disp(['Current target FDR: ',num2str(qVec(i))])
                        
                        S_TREXq = knockoff.selectVars(W_TREXq, qVec(i));
                        S_TREXf = knockoff.selectVars(W_TREXf, qVec(i));
                        S_LASSO = knockoff.selectVars(W_LASSO, qVec(i));
                        
                        if ~isempty(S_TREXq)
                            
                            fdp_TREXq(i,j) = FDP(S_TREXq);
                            tpp_TREXq(i,j) = TPP(S_TREXq);
                        end
                        
                        if ~isempty(S_TREXf)
                            
                            fdp_TREXf(i,j) = FDP(S_TREXf);
                            tpp_TREXf(i,j) = TPP(S_TREXf);
                        end
                        
                        if ~isempty(S_LASSO)
                            fdp_LASSO(i,j) = FDP(S_LASSO);
                            tpp_LASSO(i,j) = TPP(S_LASSO);
                        end
                        
                    end
                    
                    figure(23);
                    plot(qVec,sum(fdp_TREXq,2)./j,'LineWidth',3)
                    hold on
                    plot(qVec,sum(fdp_TREXf,2)./j,'LineWidth',3)
                    plot(qVec,sum(fdp_LASSO,2)./j,'LineWidth',3)
                    plot([0,1],[0,1],'k--')
                    grid on
                    legend({'W^{\Phi}','W^f','W'},'Location','NorthWest')
                    title(['FDP in dimension p: ',num2str(p),'Sample size n: ',num2str(n)])
                    
                    hold off
                    drawnow
                    
                    figure(42);
                    plot(qVec,sum(tpp_TREXq,2)./j,'LineWidth',3)
                    hold on
                    plot(qVec,sum(tpp_TREXf,2)./j,'LineWidth',3)
                    plot(qVec,sum(tpp_LASSO,2)./j,'LineWidth',3)
                    grid on
                    legend({'W^{\Phi}','W^f','W'},'Location','NorthWest')
                    title(['TPP in dimension p: ',num2str(p),'Sample size n: ',num2str(n)])
                    
                    hold off
                    drawnow
                    
                    save(['workspace_FDP_TPP_p',num2str(p),'_n',num2str(n),'_amp',num2str(10*amplitude),'_sigma',num2str(10*sigma1)]);
                    
                end
            end
        end
    end
end
