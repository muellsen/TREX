function [betaTREX, betaTREXMat, funTREXMat,out] = trex_sprox(Xmat, Yvec, inopts)
% Tuning-free Regression that adapts to the design matrix X (TREX)
% Exact solution using prox solver and enumerate all 2p solutions
% (multi-thread mode across all 2p problems) but only returns p solutions
% due to early abort (DR-Sel)
% Input:  Response Yvec in R^n
%         Data Xmat in R^nxp (n p-dimensional measurements)
%         structure of inopts with free parameters
% Output: TREX solutions: betaTREX (best),betaTREXMat (p x p matrix sorted by dimension
%         and the best amonng all +/- signs),funTREXMat (p function
%         values), additional ouput data in out struct

[n, p] = size(Xmat);

% Options for prox
defopts.abstol = 1e-10;     % Tolerance
defopts.maxit = 1e3;        % Number of iteration in Douglas Rachford
defopts.screenit = 5e1;     % Number of iteration in screening
defopts.gamma = 70;         % Arbitrary scaling
defopts.dr_mu = 1.95;       % Step size (must be in [0 2])
defopts.cpath = 0.5;        % Regularization parameter
defopts.activeSet = 1:p;    % Active set of variables
defopts.verbose = 0;        % Plot additional output
defopts.plotting = 1;       % Plot the solution trajectory

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
betaTREXMat = zeros(p, nA,nC);
funTREXMat = zeros(nA, nC);
runTimes = zeros(nA, nC);

% Build list of active variables (+-)
pActInds = sort(actInds);
maxit = opts.maxit;
abstol = opts.abstol;
beta0 = zeros(p,1);
gamma = opts.gamma;

betaProxTrace = zeros(p,2*nA,nC,maxit);
traceLenMat = zeros(nA,nC);
parfor j = 1:nA
    
    tempRunTime = zeros(nC,1);
    tempBetaMat = zeros(p,nC);
    tempFunMat = zeros(nC,1);
    betaProxTraceTemp = zeros(p,nC,maxit);
    
    for k=1:nC
        
        % Current regularization parameter
        alpha = cpath(k);
        
        x_s = Xmat(:, j);
        
        eta = x_s'*Xmat*beta0;
        
        y = Yvec+1e-3*randn(n,1);
        
        tic
        % Prox solver
        [beta_j,betaTrace] = sdougRach(eta,y,x_s,Yvec,Xmat,beta0,alpha,gamma,opts);
        tempRunTime(k) = toc;
        
        tempBetaMat(:,k) = beta_j;
        tempFunMat(k) = objTREX_SA(beta_j,Xmat,Yvec,alpha);
        % Warmstart
        %beta0 = beta_j;
        traceLen = size(betaTrace,2);
        traceLenMat(j,k) = traceLen;
        betaProxTraceTemp(:,k,1:traceLen) = betaTrace;
    end
    
    betaProxTrace(:,j,:,:) = betaProxTraceTemp;
    
    runTimes(j,:) = tempRunTime;
    betaTREXMat(:, j,:) = tempBetaMat;
    funTREXMat(j,:) = tempFunMat;
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
funTREXMat = inf(p, nC);
funTREXMat(pActInds,:) = temp(1:nA,:);

% Put beta solutions into array
temp = betaTREXMat;
betaTREXMat = zeros(p, p,nC);
betaTREXMat(:,pActInds,:) = temp(:,1:nA,:);

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
% FUNCTIONS BELOW IMPLEMENT THE DOUGLAS RACHFORD AND PROXIMITY OPERATORS
% ---------------------------------------------------------------
function [betaTREX,betaTrace] = sdougRach(eta0,y0,v,Yvec,Xmat,beta0,alpha,gamma,dropts)
% Simultaneous Douglas Rachford that checks both half spaces simultaneously

maxit = dropts.maxit;
screenit = dropts.screenit;
abstol = dropts.abstol;
mu_k = dropts.dr_mu;

[n,p] = size(Xmat);

% n+1 x p matrix
M1 = [v'*Xmat;Xmat];
M2 = M1;
M2(1,:) = -M2(1,:);

% Matrices
Id_n = eye(n+1);
R1 = M1'/(Id_n+M1*M1');
R2 = R1;
R2(:,1) = -R2(:,1);

% Initialize DR sequences for both sign v's
x1_k = beta0;
y1_k = [eta0;y0];

x2_k = x1_k;
y2_k = y1_k;
y2_k(1) = -y2_k(1);

v1 = v;
v2 = -v;

betaTrace1 = zeros(p,screenit);
betaTrace2 = zeros(p,screenit);

% Pre-screening phase
for kk=1:screenit
    
    % Update first variable
    q1_k = M1*x1_k - y1_k;
    b1_k = x1_k - R1*q1_k;
    c1_k = M1*b1_k;
    z1_k = proxST(2*b1_k-x1_k,gamma);
    t1_k = proxgTREX(2*c1_k-y1_k,v1,Yvec,alpha,gamma);
    
    % Update
    x1_k = x1_k + mu_k*(z1_k-b1_k);
    y1_k = y1_k + mu_k*(t1_k-c1_k);
    
    betaTrace1(:,kk) = b1_k;
    
    % Update second variable
    q2_k = M2*x2_k - y2_k;
    b2_k = x2_k - R2*q2_k;
    c2_k = M2*b2_k;
    z2_k = proxST(2*b2_k-x2_k,gamma);
    t2_k = proxgTREX(2*c2_k-y2_k,v2,Yvec,alpha,gamma);
    
    % Update state variables
    x2_k = x2_k + mu_k*(z2_k-b2_k);
    y2_k = y2_k + mu_k*(t2_k-c2_k);
    
    betaTrace2(:,kk) = b2_k;
    
end

f1 = objTREX_SA(b1_k,Xmat,Yvec,alpha);
f2 = objTREX_SA(b2_k,Xmat,Yvec,alpha);

if f1<f2
    M = M1;
    R = R1;
    v = v1;
    x_k = x1_k;
    y_k = y1_k;
    betaTraceScreen = betaTrace1;
    
    K_old = screenit;
    b_old = b1_k;
    b_k1 = b1_k;
    
else
    M = M2;
    R = R2;
    v = v2;
    x_k = x2_k;
    y_k = y2_k;
    
    betaTraceScreen = betaTrace2;
    K_old = screenit;
    b_old = b2_k;
    b_k1 = b2_k;
end

maxit = maxit-screenit;

betaTrace = zeros(p,maxit);

for kk=1:maxit
    q_k = M*x_k - y_k;
    b_k = x_k - R*q_k;
    c_k = M*b_k;
    z_k = proxST(2*b_k-x_k,gamma);
    t_k = proxgTREX(2*c_k-y_k,v,Yvec,alpha,gamma);
    
    x_k1 = x_k + mu_k*(z_k-b_k);
    y_k1 = y_k + mu_k*(t_k-c_k);
    
    if norm(b_k1-b_k)<abstol || ((kk>screenit) && norm(y_k1-y_k)<abstol)
        x_k = x_k1;
        y_k = y_k1;
        break
    end
    b_k1 = b_k;
    betaTrace(:,kk) = b_k;
    
    x_k = x_k1;
    y_k = y_k1;
    
    if dropts.plotting
        figure(23);
        if ~mod(kk,2e1)
            figure(23);
            semilogy([K_old,kk],abs([b_old,b_k]),'k-','LineWidth',5)
            grid on
            hold on
            drawnow
            b_old = b_k;
            K_old = kk;
        end
    end
    % Adaptive mu_k
    % mu_k = 0.1+(0.9)*mu_k;
end

hold off
betaTREX = b_k;
betaTrace = betaTrace(:,1:kk);
betaTrace = [betaTraceScreen,betaTrace];
end


function betaT = proxST(b,gamma)
betaT = sign(b).*max(abs(b) - gamma,0);
end

% TREX prox
function etaProx_yProx = proxgTREX(eta_y,v,Yvec,alpha,gamma)

eta = eta_y(1);
y = eta_y(2:end);

gammaPrime = 2*gamma/alpha;

% Check whether prox calculation is needed

yDiff = (y-Yvec);
yDiffNorm2 = sum(yDiff.^2);

xTYvec = v'*Yvec;

if (2*gamma*eta+yDiffNorm2)>(2*gamma*xTYvec)
    
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
