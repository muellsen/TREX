function [W,Z,trexTrace,trexopts] = trexSignedMax(X,X_ko,y)
% X:    Original design
% X_ko: Knockoff features
% y:    response

% Combined feature matrix
Xunorm = [X,X_ko];

% Response vector
Y = y;

[n,p] = size(Xunorm);

% Standardize data
X0 = Xunorm-repmat(mean(Xunorm),n,1);

% Normalize X to length sqrt(n)
normX = repmat(sqrt(sum(X0.^2)),n,1);
X = sqrt(n)*X0./normX;

% Run the TREX

% qTREX regularization path
cpath = [1.5:-0.005:0.01];

% Number of bootstrap runs
L = 0;

% TREX options
trexopts.boots = L;
trexopts.rep = 6;
trexopts.verbose = 0;
trexopts.q = 40;                    % Order of the norm to approximate max norm
trexopts.beta0 = zeros(p,1);        % Include the zero vector as start point
trexopts.cpath = cpath;  % Regularization path

trexTrace = trexp(X,Y,trexopts);

% Thresholding (added feature)
epsPSG = 1e-6; % Previous threshold 1e-3,1e-4
trexTrace(abs(trexTrace(:))<epsPSG) = 0;

figure(1);
plot(cpath,sum(trexTrace~=0),'LineWidth',10)
grid on
hold on
xlabel('c regularization path')
ylabel('Sparsity of the solution')
drawnow

% % Check the order in which variables enter the TREX path
% % S is the list of variable indices
% S = [];
% for i=1:length(cpath)
%     % Check non-zero indices
%     currInds = find(trexTrace(:,i)~=0);
%     % Check which ones are new
%     addInds = setdiff(currInds,S)';
%     % Sort them by coefficient value
%     [~,sortedInds] = sort(abs(trexTrace(addInds,i)),'descend');
%     S = [S,addInds(sortedInds)];
% end
% 
% % Fill up the final variables
% currInds = [1:p]';
% addInds = setdiff(currInds,S)';
% S = [S,addInds];
% 
% % Copied from other statistics
% added = S;
% [~,order] = sort(added);
% 
% p = size(X_ko,2);
% Z = 2*p + 1 - order;
% orig = 1:p; ko = (p+1):(2*p);
% W = max(Z(orig), Z(ko)) .* sign(Z(orig) - Z(ko));

% Put the maximum c for each variable when it enters the path
Z = zeros(2*p,1);
for i=1:length(cpath)
    % Check non-zero indices
    currInds = find(trexTrace(:,i)~=0);
    % Check which ones are new
    addInds = setdiff(currInds,find(Z~=0))';
    Z(addInds) = cpath(i);
end

% Original size of the problem
p = size(X_ko,2);
orig = 1:p; ko = (p+1):(2*p);
W = max(Z(orig), Z(ko)) .* sign(Z(orig) - Z(ko));

end

