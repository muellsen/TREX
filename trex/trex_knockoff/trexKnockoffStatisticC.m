function [W,Z,trexTrace,trexSCSopts] = trexKnockoffStatisticC(X,X_ko,y)
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

% Run the convex TREX

cpath = min(p./n,1/2);  % Regularization path on one value

% %% Convex TREX options for SCS solver
% trexSCSopts.verbose = 0;    % No diagnostic output
% trexSCSopts.maxit = 5000;   % Number of iterations in SCS solver
% trexSCSopts.abstol = 1e-4;  % Solution accuracy in SCS solver
% trexSCSopts.cpath = cpath;  % Regularization path
% 
% % Multi-thread version
% [~,trexTrace,trexFunTrace] = trex_scsp_direct(X,Y,trexSCSopts);

%% Convex TREX options for ECOS solver

trexEcosopts.cpath = cpath;          % Constant c
trexEcosopts.verbose = 0;           % No diagnostic output
trexEcosopts.feastol = 1e-6;        % Tolerances in the solver
trexEcosopts.reltol = 1e-6;
trexEcosopts.abstol = 1e-6;

% Single thread version
[~,trexTrace,trexFunTrace] = trex_ecosp(X,Y,trexEcosopts);


minFunTrace = min(trexFunTrace(1:2:end-1),trexFunTrace(2:2:end));
[~,order]=sort(minFunTrace);
p = size(X_ko,2);
Z = 2*p + 1 - order;
Z = 1./minFunTrace;
orig = 1:p; ko = (p+1):(2*p);
W = max(Z(orig), Z(ko)) .* sign(Z(orig) - Z(ko));

end

