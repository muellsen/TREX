function [betaTREX, betaTREXMat, funTREXMat,out] = trex_gprox(Xmat, Yvec, inopts)
% Tuning-free Regression that adapts to the design matrix X (TREX)
% Generalized TREX is based on a generalized objective function 
% Exact solution using prox solver and enumerate all 2p solutions
% (multi-thread mode across all 2p problems)
% Input:  Response Yvec in R^n
%         Data Xmat in R^nxp (n p-dimensional measurements)
%         structure of inopts with free parameters
% Output: TREX solutions: betaTREX (best),betaTREXMat (2*p x p matrix sorted by dimension
%         and +/- signs),funTREXMat (2*p function values), additional
%         output in out

[n, p] = size(Xmat);

% Options for prox
defopts.abstol = 1e-6;      % Tolerance
defopts.maxit = p*1e2;      % Number of iteration in Douglas Rachford
defopts.gamma = 10;         % Arbitrary scaling
defopts.dr_mu = 1.9;        % Step size (must be in [0 2])
defopts.cpath = 0.5;        % Regularization parameter
defopts.activeSet = 1:p;    % Active set of variables
defopts.verbose = 0;        % Plot additional output
defopts.plotting = 1;       % Plot the solution trajectory
defopts.qPower = 2;         % qPower: TREX is 2, all others are generalizations
defopts.outit = 1e2;        % Interval for writing out prox trace

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

% TREX exponent
qPower = opts.qPower;

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
beta_j = beta0;
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
        % Prox solver with warm start
        [beta_j,betaTrace] = dougRach(eta,y,x_s,Yvec,Xmat,beta0,alpha,qPower,gamma,opts);
        tempRunTime(k) = toc;
        
        tempBetaMat(:,k) = beta_j;
        tempFunMat(k) = objTREX_SA(beta_j,Xmat,Yvec,alpha,qPower);
        
        % Warmstart (not used right now)
        %beta0 = beta_j;
        
        traceLen = size(betaTrace,2);
        traceLenMat(i,k) = traceLen;
        betaProxTraceTemp(:,k,1:traceLen) = betaTrace;
        if opts.verbose
            disp([num2str(j),'th variable done...'])
        end
    end
    
    betaProxTrace(:,i,:,:) = betaProxTraceTemp;
    
    runTimes(i,:) = tempRunTime;
    betaTREXMat(:, i,:) = tempBetaMat;
    funTREXMat(i,:) = tempFunMat;
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
% FUNCTIONS BELOW IMPLEMENT DOUGLAS RACHFORD AND THE PROXIMITY OPERATOR
% ---------------------------------------------------------------
function [betaTREX,betaTrace] = dougRach(eta0,y0,v,Yvec,Xmat,beta0,alpha,qPower,gamma,dropts)

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
    z_k = proxST(2*b_k-x_k,gamma);
    t_k = proxggTREX(2*c_k-y_k,v,Yvec,alpha,qPower,gamma);
    
    x_k1 = x_k + mu_k*(z_k-b_k);
    y_k1 = y_k + mu_k*(t_k-c_k);
    
    if norm(b_k1-b_k)<abstol || ((kk>1e2) && norm(y_k1-y_k)<abstol)
        x_k = x_k1;
        y_k = y_k1;
        %K
        %disp('Early convergence')
        break
    end
    b_k1 = b_k;
    
    if ~mod(kk,outit)
        betaTrace(:,outk) = b_k;
        ok = outk+1;
    end
    
    x_k = x_k1;
    y_k = y_k1;
    
    %     if dropts.plotting
    %
    %         if ~mod(kk,2e1)
    %             semilogy([K_old,kk],abs([b_old,b_k]),'k-','LineWidth',5)
    %             grid on
    %             hold on
    %             drawnow
    %             b_old = b_k;
    %             K_old = kk;
    %         end
    %     end
    %mu_k = 0.1+(0.9)*mu_k;
end
%pause
%hold off
betaTREX = b_k;
betaTrace = betaTrace(:,1:(outk-1));
end

function betaT = proxST(b,gamma)
betaT = sign(b).*max(abs(b) - gamma,0);
end

% TREX prox
function etaProx_yProx = proxggTREX(eta_y,v,Yvec,alpha,qPower,gamma)

eta = eta_y(1);
y = eta_y(2:end);

xTYvec = v'*Yvec;

% Shifted eta
etaPrime = eta-xTYvec;

% Dual norm
q_s = qPower/(qPower-1);

if (abs(q_s-round(q_s)))>1e-6 && qPower~=3
    error('Dual norm of q must be integer!')
end

% Factor
delta = (alpha*(1-(1/q_s)))^(q_s-1);

% Check whether prox calculation is needed

% First summand of prox check
q_t = q_s*gamma^(q_s-1)*etaPrime;

yDiff = (y-Yvec);
yDiffNorm2 = sum(yDiff.^2);
yDiffNorm = sqrt(yDiffNorm2);

testVal = (q_t + delta*(yDiffNorm.^q_s));

% Check whether prox calculation is needed

if testVal>0
    
    % Solve the q* root using matlab
    if qPower==3
        % Highest root + 1
        rootfacs = zeros(1,5);
        rootfacs(5) = 1;
        rootfacs(3) = q_s.*etaPrime/(gamma*delta);
        rootfacs(2) = q_s/(delta^2); % for q=q_s=2 we already have a component
        rootfacs(1) = -q_s*yDiffNorm./(gamma*delta^2);
    else
        q_s = round(q_s);
        
        % Highest root + 1
        rootfacs = zeros(1,2*q_s);
        rootfacs(2*q_s) = 1;
        rootfacs(q_s) = q_s.*etaPrime/(gamma*delta);
        rootfacs(2) = rootfacs(2)+q_s/(delta^2); % for q=q_s=2 we already have a component
        rootfacs(1) = -q_s*yDiffNorm./(gamma*delta^2);
        
    end
    
    % Reverse root factors (factor for highest degree first,...)
    rootfacs = rootfacs(end:-1:1);
    qroots = roots(rootfacs);
    
    % Only the real root
    t = qroots(imag(qroots)==0);
    
    if sum(t>=0)>1
        error('Too many positive real roots')
    end
    
    % Take the largest real root
    t = max(t);
    
    if qPower==3
        t = sqrt(t);
    end
    
    % p_denom  = gamma+delta*(eta+gamma*delta*(t.^q_s)/q_s)*t.^(q_s-2);
    % p = y/(p_denom);
    
    % Simplified prox computation
    p = (yDiff/yDiffNorm) * t;
    
    etaProx_yProx = [eta + gamma*delta*(t.^q_s)/q_s;y-gamma*p];;
    
else
    
    etaProx_yProx = [xTYvec;Yvec];
    
end
end


function [funVal,L1beta] = objTREX_SA(betaVec,X,Y,c,qPower)
% Generalized TREX stand-alone objective function
% funVal: function value at current betaVec
% L1beta: L1 norm of betaVec

% Global data vector X, response Y, regularization parameter c, exponent q

[n,p]= size(X);

funVal = (sqrt(sum((Y-X*betaVec).^2))^qPower)/(c*max(abs(X'*(Y-X*betaVec))).^(qPower-1));

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
