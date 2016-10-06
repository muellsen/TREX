function [betaTREX, betaTREXMat, funTREXMat,out] = trex_prox(Xmat, Yvec, inopts)
% Tuning-free Regression that adapts to the design matrix X (TREX)
% Exact solution using prox solver and enumerate all 2p solutions
% (multi-thread mode across all 2p problems)
% Input:  Response Yvec in R^n
%         Data Xmat in R^nxp (n p-dimensional measurements)
%         structure of inopts with free parameters
% Output: TREX solutions: betaTREX (best),betaTREXMat (2*p x p matrix sorted by dimension
%         and +/- signs),funTREXMat (2*p function values)

[n, p] = size(Xmat);

% Options for prox
defopts.abstol = 1e-6;      % Tolerance
defopts.maxit = p*1e3;        % Number of iteration in Douglas Rachford
defopts.gamma = 10;         % Arbitrary scaling
defopts.dr_mu = 1.9;        % Step size (must be in [0 2])
defopts.cpath = 0.5;        % Regularization parameter
defopts.activeSet = 1:p;    % Active set of variables
defopts.verbose = 0;        % Plot additional output
defopts.plotting = 1;       % Plot the solution trajectory
defopts.outit = 1e1;        % Interval for writing out prox trace

% Setting a single c is also possible
if ~isfield(inopts,'cpath') && isfield(inopts,'c')
    inopts.cpath = inopts.c;
end

% Merge options inopts and defopts
if nargin < 3 || isempty(inopts) % no input options available
    opts = defopts;
else
    opts = getoptions(inopts, defopts);
end

% Active set on which TREX is conditioned on
actInds = opts.activeSet;
nA = length(actInds);

% TREX parameter
cpath = opts.cpath;

% Length of the c path
nC = length(cpath);

% Create storage for 2 nA beta vectors, function values and timings
betaTREXMat = zeros(p, 2*nA,nC);
funTREXMat = zeros(2*nA, nC);
runTimes = zeros(2*nA, nC);

% Build list of active variables (+-)
pActInds = sort([2*(actInds-1)+1,2*actInds]);
maxit = opts.maxit;
abstol = opts.abstol;
beta0 = zeros(p,1);
gamma = opts.gamma;
outit = opts.outit;

nIt = round(maxit/outit);

betaProxTrace = zeros(p,2*nA,nC,nIt);
traceLenMat = zeros(2*nA,nC);
parfor i = 1:2*nA
    
    ii = pActInds(i);
    
    j=round(ii/2);
    
    s = (-1)^i*(-1);
    
    tempRunTime = zeros(nC,1);
    tempBetaMat = zeros(p,nC);
    tempFunMat = zeros(nC,1);
    betaProxTraceTemp = zeros(p,nC,nIt);
    
    for k=1:nC
        
        % Current regularization parameter
        alpha = cpath(k);
        
        x_s = s * Xmat(:, j);
        
        eta = x_s'*Xmat*beta0;
        
        y = Yvec+1e-3*randn(n,1);
        
        tic
        % Prox solver (feed -x_s because of flipped constraint term)
        [beta_j,betaTrace] = dougRach(eta,y,x_s,Yvec,Xmat,beta0,alpha,gamma,opts);
        tempRunTime(k) = toc;
        
        tempBetaMat(:,k) = beta_j;
        tempFunMat(k) = objTREX_SA(beta_j,Xmat,Yvec,alpha);
        % Warmstart
        %beta0 = beta_j;
        traceLen = size(betaTrace,2);
        traceLenMat(i,k) = traceLen;
        betaProxTraceTemp(:,k,1:traceLen) = betaTrace;
    end
    
    betaProxTrace(:,i,:,:) = betaProxTraceTemp;
    
    runTimes(i,:) = tempRunTime;
    betaTREXMat(:, i,:) = tempBetaMat;
    funTREXMat(i,:) = tempFunMat;
    if opts.verbose
        disp([num2str(j),'th variable done...'])
    end
end

% Find best solution for each value of c
[~, ihat] = min(funTREXMat,[],1);

betaTREX=zeros(p,nC);
for k=1:nC
    betaTREX(:,k) = betaTREXMat(:, ihat(k),k);
end

% Put function values into array
temp = funTREXMat;
funTREXMat = inf(2 * p, nC);
funTREXMat(pActInds,:) = temp(1:2*nA,:);

% Put beta solutions into array
temp = betaTREXMat;
betaTREXMat = zeros(p, 2 * p,nC);
betaTREXMat(:,pActInds,:) = temp(:,1:2*nA,:);

% Squeeze to proper dimension for scalar c
if nC==1
    betaTREX = squeeze(betaTREX);
    betaTREXMat = squeeze(betaTREXMat);
    funTREXMat = squeeze(funTREXMat);
end

betaProxTrace = betaProxTrace(:,:,:,1:max(traceLenMat(:)));
betaProxTrace = squeeze(betaProxTrace);

out.opts = opts;
out.runTimes = runTimes;
out.betaProxTrace = betaProxTrace;
out.X = Xmat;
out.Y = Yvec;
end

% ---------------------------------------------------------------
% FUNCTIONS BELOW IMPLEMENT THE DOUGLAS RACHFORD AND PROXIMITY
% ---------------------------------------------------------------
function [betaTREX,betaTrace] = dougRach(eta0,y0,v,Yvec,Xmat,beta0,alpha,gamma,dropts)

maxit = dropts.maxit;
outit = dropts.outit;
nit = round(maxit/outit);
abstol = dropts.abstol;
mu_k = dropts.dr_mu;

[n,p] = size(Xmat);

% n+1 x p matrix
M = [v'*Xmat;Xmat];

% Matrices
Id_n = eye(n+1);
R = M'/(Id_n+M*M');

%         if opts.plotting
%             figure(42)
%             imagesc(R)
%             colorbar
%         end
%
% Initialize DR sequences
x_k = beta0;
y_k = [eta0;y0];

% Step size
K_old = 1;
b_old = beta0;
b_k1 = b_old;
betaTrace = zeros(p,maxit);
outk=1;

for kk=1:maxit
    
    q_k = M*x_k - y_k;
    b_k = x_k - R*q_k;
    c_k = M*b_k;
    z_k = proxST(2*b_k-x_k,alpha*gamma);
    t_k = proxgTREX(2*c_k-y_k,v,Yvec,1,gamma);
    
    %temp_tk = proxgTREX2(2*c_k-y_k,v,Yvec,alpha,gamma)
    
    %diffT = norm(t_k-alpha/2*temp_tk)
    %pause
    
    x_k1 = x_k + mu_k*(z_k-b_k);
    y_k1 = y_k + mu_k*(t_k-c_k);
    
    if norm(b_k1-b_k)<abstol || ((kk>1e1) && norm(y_k1-y_k)<abstol)
        x_k = x_k1;
        y_k = y_k1;
        %K
        %disp('Early convergence')
        break
    end
    
    b_k1 = b_k;
    
    if ~mod(kk,outit)
        betaTrace(:,outk) = b_k;
        outk = outk+1;
    end
    x_k = x_k1;
    y_k = y_k1;
    
    if dropts.plotting
        
        if ~mod(kk,2e1)
            semilogy([K_old,kk],abs([b_old,b_k]),'k-','LineWidth',5)
            grid on
            hold on
            drawnow
            b_old = b_k;
            K_old = kk;
        end
    end
    %mu_k = 0.1+(0.9)*mu_k;
end
%pause
hold off
betaTREX = b_k;
betaTrace = betaTrace(:,1:(outk-1));
end

function betaT = proxST(b,gamma)
betaT = sign(b).*max(abs(b) - gamma,0);
end

% TREX prox
function etaProx_yProx = proxgTREX(eta_y,v,Yvec,alpha,gamma)

eta = eta_y(1);
y = eta_y(2:end);

yDiff = (y-Yvec);
yDiffNorm2 = sum(yDiff.^2);

% Shift from TREX
xTYvec = v'*Yvec;

% Shifted eta
etaPrime = eta-xTYvec;

% Check whether prox calculation is needed
if (4*gamma*etaPrime+alpha*yDiffNorm2)>0
        
    yDiffNorm = sqrt(yDiffNorm2);
    
    % Compute mu
    mu = (4/alpha^2)*sqrt(yDiffNorm2/(gamma.^2) + 32/(27*alpha^2)*((alpha*etaPrime)/(2*gamma)+1)^3);
    
    % Compute root explicitly
    p_s1 = gamma + (alpha*etaPrime)/2;
    
    rt3term = 4*yDiffNorm/(alpha^2*gamma);
    
    p_denom = p_s1 + ((alpha^2*gamma)/8)*(sign(rt3term + mu)*(abs(rt3term + mu))^(1/3) ...
        + sign(rt3term - mu)*(abs(rt3term - mu))^(1/3)).^2;
    %yDiffNorm
    p = yDiff./(p_denom);
    
    % Root checking
    % t = yDiffNorm/p_denom;
    % t^3+4*(alpha*etaPrime+2*gamma)/(alpha^2*gamma)*t - 8*yDiffNorm/(alpha^2*gamma)
    
    % Compute polynomial
    etaProx_yProx = [eta + alpha*gamma*sum((p).^2)/4;y-gamma*p];
    
else
    
    etaProx = xTYvec;
    yProx = Yvec;
    etaProx_yProx = [etaProx;yProx];
    
end
end

function etaProx_yProx = proxgTREX_dep(eta_y,v,Yvec,alpha,gamma)
% Deprecated prox function

eta = eta_y(1);
y = eta_y(2:end);

gammaPrime = 2*gamma/alpha;

% Check whether prox calculation is needed

yDiff = (y-Yvec);
yDiffNorm2 = sum(yDiff.^2);

xTYvec = v'*Yvec;

if (2*gamma*eta+yDiffNorm2)>(2*gamma*xTYvec)
    
    %disp(['Hello with eta=',num2str(eta)])
    yDiffNorm = sqrt(yDiffNorm2);
    
    c1 = (gammaPrime + eta - xTYvec);
    
    % Compute mu
    mu = sqrt(yDiffNorm2./gammaPrime^2 + (8/27*(c1./gammaPrime).^3));
    
    % Compute root explicitly
    p_denom = (c1 + gammaPrime/2*((sign(yDiffNorm./gammaPrime+mu)*(abs(yDiffNorm./gammaPrime+mu)).^(1/3) ...
        + sign(yDiffNorm./gammaPrime-mu)*(abs(yDiffNorm./gammaPrime-mu)).^(1/3))).^2);
    
    p = yDiff./p_denom;
    
    %p_norm = sqrt(sum(p.^2));
    %cubEq = p_norm.^3 + 2*(c1)*p_norm - 2*yDiffNorm./gamma
    %pause
    
    %             pr1 = yDiff/p_denom1;
    %             pr2 = yDiff/p_denom2;
    %             pr3 = yDiff/p_denom3;
    %
    %             %p1 = yDiff/p_denom1;
    %             %p2 = yDiff/p_denom2;
    %
    %             pr1_norm = sqrt(sum(pr1.^2));
    %             pr2_norm = sqrt(sum(pr2.^2));
    %             pr3_norm = sqrt(sum(pr3.^2));
    %
    %             cubEqR1 = pr1_norm.^3 + 2*(c1)*pr1_norm - 2*yDiffNorm;
    %             cubEqR2 = pr2_norm.^3 + 2*(c1)*pr2_norm - 2*yDiffNorm;
    %             cubEqR3 = pr3_norm.^3 + 2*(c1)*pr3_norm - 2*yDiffNorm;
    
    %             % Solve the cubic root using matlab
    %             rootfacs = [1 0 2*(c1) -2*yDiffNorm];
    %             cubroots = roots(rootfacs);
    %
    %             p1 = yDiff./yDiffNorm*cubroots(1);
    %             p2 = yDiff./yDiffNorm*cubroots(2);
    %             p3 = yDiff./yDiffNorm*cubroots(3);
    %
    %             p1_norm = sqrt(sum(p1.^2));
    %             p2_norm = sqrt(sum(p2.^2));
    %             p3_norm = sqrt(sum(p3.^2));
    %
    %             cubEq1 = p1_norm.^3 + 2*(c1)*p1_norm - 2*yDiffNorm;
    %             cubEq2 = p2_norm.^3 + 2*(c1)*p2_norm - 2*yDiffNorm;
    %             cubEq3 = p3_norm.^3 + 2*(c1)*p3_norm - 2*yDiffNorm;
    %
    %             if abs(cubEq1)<1e-3
    %                 p = p1;
    %             elseif abs(cubEq2)<1e-3
    %                 p = p2;
    %             elseif abs(cubEq3)<1e-3
    %                 p = p3;
    %             else
    %                 error('No root is good')
    %             end
    
    %etaProx = eta + gammaPrime*sum((p).^2)/2;y-gammaPrime*p;
    %yProx = y-gammaPrime*p;
    %etaProx_yProx = [etaProx;yProx];
    etaProx_yProx = [eta + gammaPrime*sum((p).^2)/2;y-gammaPrime*p];
    
else
    
    etaProx = xTYvec;
    yProx = Yvec;
    etaProx_yProx = [etaProx;yProx];
    
end
end

function [funVal,L1beta] = objTREX_SA(betaVec,X,Y,c)
% TREX standalone objective function
% funVal: function value at current betaVec
% L1beta: L1 norm of betaVec

% Global data vector X, response Y, regularization parameter c

[n,p]= size(X);

funVal = sum((Y-X*betaVec).^2)/(c*max(abs(X'*(Y-X*betaVec))));

L1beta = sum(abs(betaVec));

funVal = funVal+L1beta;
end

% ---------------------------------------------------------------
% FUNCTIONS BELOW ARE TAKEN FROM Niko Hansen's CMA-ES code
% ---------------------------------------------------------------
function opts=getoptions(inopts, defopts)
% OPTS = GETOPTIONS(INOPTS, DEFOPTS) handles an arbitrary number of
% optional arguments to a function. The given arguments are collected
% in the struct INOPTS.  GETOPTIONS matches INOPTS with a default
% options struct DEFOPTS and returns the merge OPTS.  Empty or missing
% fields in INOPTS invoke the default value.  Fieldnames in INOPTS can
% be abbreviated.
%
% The returned struct OPTS is first assigned to DEFOPTS. Then any
% field value in OPTS is replaced by the respective field value of
% INOPTS if (1) the field unambiguously (case-insensitive) matches
% with the fieldname in INOPTS (cut down to the length of the INOPTS
% fieldname) and (2) the field is not empty.
%


if nargin < 2 || isempty(defopts) % no default options available
    opts=inopts;
    return;
elseif isempty(inopts) % empty inopts invoke default options
    opts = defopts;
    return;
elseif ~isstruct(defopts) % handle a single option value
    if isempty(inopts)
        opts = defopts;
    elseif ~isstruct(inopts)
        opts = inopts;
    else
        error('Input options are a struct, while default options are not');
    end
    return;
elseif ~isstruct(inopts) % no valid input options
    error('The options need to be a struct or empty');
end

opts = defopts; % start from defopts
% if necessary overwrite opts fields by inopts values
defnames = fieldnames(defopts);
idxmatched = []; % indices of defopts that already matched
for name = fieldnames(inopts)'
    name = name{1}; % name of i-th inopts-field
    if isoctave
        for i = 1:size(defnames, 1)
            idx(i) = strncmp(lower(defnames(i)), lower(name), length(name));
        end
    else
        idx = strncmp(lower(defnames), lower(name), length(name));
    end
    if sum(idx) > 1
        error(['option "' name '" is not an unambigous abbreviation. ' ...
            'Use opts=RMFIELD(opts, ''' name, ...
            ''') to remove the field from the struct.']);
    end
    if sum(idx) == 1
        defname  = defnames{find(idx)};
        if ismember(find(idx), idxmatched)
            error(['input options match more than ones with "' ...
                defname '". ' ...
                'Use opts=RMFIELD(opts, ''' name, ...
                ''') to remove the field from the struct.']);
        end
        idxmatched = [idxmatched find(idx)];
        val = getfield(inopts, name);
        % next line can replace previous line from MATLAB version 6.5.0 on and in octave
        % val = inopts.(name);
        if isstruct(val) % valid syntax only from version 6.5.0
            opts = setfield(opts, defname, ...
                getoptions(val, getfield(defopts, defname)));
        elseif isstruct(getfield(defopts, defname))
            % next three lines can replace previous three lines from MATLAB
            % version 6.5.0 on
            %   opts.(defname) = ...
            %      getoptions(val, defopts.(defname));
            % elseif isstruct(defopts.(defname))
            warning(['option "' name '" disregarded (must be struct)']);
        elseif ~isempty(val) % empty value: do nothing, i.e. stick to default
            opts = setfield(opts, defnames{find(idx)}, val);
            % next line can replace previous line from MATLAB version 6.5.0 on
            % opts.(defname) = inopts.(name);
        end
    else
        warning(['option "' name '" disregarded (unknown field name)']);
    end
end
end

% ---------------------------------------------------------------
% ---------------------------------------------------------------
function res = isoctave
% any hack to find out whether we are running octave
s = version;
res = 0;
if exist('fflush', 'builtin') && eval(s(1)) < 7
    res = 1;
end
end
