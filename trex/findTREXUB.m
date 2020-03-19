function [c_max,cVec,u_hatVec,u_max] = findTREXUB(X,Y)
% Function that determines the upper bound on the TREX c-path
% cf. Assumption... Bien et al, 2018

% General path for non-convex TREX
cpath = [0.0001,0.01:0.01:2];
nPath = length(cpath);

% Try to get a rough idea where the bound is using the fast nonconvex TREX

% Non-convex TREX
clear trexncopts;
trexncopts.cpath = cpath;
trexncopts.rep = 10;
trexncopts.verbose = 0;
tic;[betaTREXnc,betaTREXMatnc,funVecTREXnc] = trexp(X,Y,trexncopts); solveQTREX = toc;

% Upper bound
u_max = max(abs(X'*Y))

% Compute u_hat trace
u_hatnc = max(abs(X'*(repmat(Y,1,nPath)-X*betaTREXnc)))

% Index of first u_hat value that does not satisfy the bound (note: some of high c
% values may lead to near sparse vectors that satisfy the bound again; we
% are only interested in beta solutions that have non-zero support

maxInd = find(u_hatnc>(u_max),1);

% In case all solutions satisfy the bound, start at c=1;
if isempty(maxInd)
   maxInd = (nPath-1)/2; 
end

% Use the maximum as start point of a step search (+1/-1 dependent on the
% value)
firstInd = maxInd

c = cpath(firstInd);

% Convex TREX options for ECOS solver
trexEcosopts.cpath = c;  % Constant c
trexEcosopts.verbose = 0;           % No diagnostic output
trexEcosopts.feastol = 1e-12;        % Tolerances in the solver
trexEcosopts.reltol = 1e-12;
trexEcosopts.abstol = 1e-12;

% Multi-thread version
tic;betaTREX = trex_ecosp(X,Y,trexEcosopts); solveCTREX = toc

u_hat = max(abs(X'*(Y-X*betaTREX)))

if u_hat<u_max
    incStep = +1; % Current c too large
else
    incStep = -1; % Current c too large
end

% Record computed u_hat and corresponding c values
u_hatVec = [];
u_hatVec = [u_hatVec,u_hat];
cVec = [];
cVec = [cVec,c];
c_max = c;

% Indicator whether search was successful
foundC = 0;

currInd = firstInd;

while foundC==0 && currInd<=nPath
    
    currInd = currInd + incStep;
    c = cpath(currInd);

    % Convex TREX options for ECOS solver
    trexEcosopts.cpath = c;  % Constant c
    trexEcosopts.verbose = 0;           % No diagnostic output
    trexEcosopts.feastol = 1e-12;        % Tolerances in the solver
    trexEcosopts.reltol = 1e-12;
    trexEcosopts.abstol = 1e-12;
    
    % Multi-thread version
    tic;betaTREX= trex_ecosp(X,Y,trexEcosopts); solveCTREX = toc;
    
    u_hat = max(abs(X'*(Y-X*betaTREX)));
   

    % Check conditions on whether the bound was identified
    if u_hat < u_max && incStep==-1
        foundC = 1;
        c_max = c;
        cVec = [cVec,c];
        u_hatVec = [u_hatVec,u_hat];

    elseif u_hat > u_max && incStep==1
        foundC = 1;
        c_max = cpath(currInd-1);
    else
        foundC = 0;
        cVec = [cVec,c];
        u_hatVec = [u_hatVec,u_hat];
    end
   
end





