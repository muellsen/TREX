function [betaTREX, betaTREXMat, funTREXMat,out] = trex_scsp(Xmat, Yvec, inopts)
% Tuning-free Regression that adapts to the design matrix Xmat (TREXmat)
% Exact solution using scs sdp solver and enumerate all 2p solutions
% (multi-thread mode to solve 2p problems)
% Input:  Response Yvec in R^n
%         Data Xmat in R^nxp (n p-dimensional measurements)
%         structure of inopts with free parameters
% Output: TREX solutions: betaTREX matrix of best solutions for each c,betaTREXmatMat (2*p x p matrix sorted by dimension
%         and +/- signs),funTREXmatMat (2*p function values)

% Size of the input desgin matrix
[n, p] = size(Xmat);

% Options for c-trex
%defopts.feastol = 1e-12;
%defopts.reltol = 1e-12;
defopts.abstol = 1e-2;  % eps parameter is solution tolerance in SCS solver
defopts.maxit = 5000;
defopts.cpath = 0.5;
defopts.activeSet = 1:p;
defopts.verbose = 0;

% Merge options inopts and defopts
if nargin < 3 || isempty(inopts) % no input options available
    opts = defopts;
else
    opts = getoptions(inopts, defopts);
end

% Set explicit options
maxit = opts.maxit;
verbose = opts.verbose;
cpath = opts.cpath;

% Length of the c path
nC = length(cpath);

% Active set on which TREX is conditioned on
actInds = opts.activeSet;
nA = length(actInds);

v = sparse(n, 1); % initialization
% Data initialization for scs
Ip = speye(p);In = speye(n);
AA = [sparse(n, p + 1), Xmat, In];
bb_scs = Yvec;
GG = [sparse(2 * p,1), [-Ip; -Ip], [Ip; -Ip], sparse(2 * p, n)];
% second-order cone part:
GG = [
    GG;
    [-1, sparse(1, 2 * p), -v'];
    [-1, sparse(1, 2 * p), v'];
    [sparse(n, 2 * p + 1), 2 * In]
    ];
hh_scs = zeros(size(GG, 1), 1);
cc = [ones(p + 1, 1); zeros(p + n, 1)];
%   data must consist of data.A, data.b, data.c, where A,b,c used as above.
A=[AA;GG];
b=[bb_scs;hh_scs];
c2=cc;

% From here need to change for SCS
% Set up K -  a product of cones in the following order:
% free cone, lp cone, second order cone(s), semi-definite cone(s), primal
%   exponential cones, dual exponential cones

%   cone.f, length of free cone (for equality constraints)
K.f=n;
%   cone.l, length of lp cone
K.l=2*p;
%   cone.q, array of SOC lengths
K.q=[n+2]; %??not sure what is meant by the array of lengths, I believe I only have 1 SOC of size n+2 (n from Yvec-Xmatb, 1 from t_0 and 1 for RHS)
%   cone.s, array of SD lengths
K.s=[];
% no exponential cones as far as I know (number of primal and dual exp cones)
K.ep=0;
K.ed=0;

% Set the parameters for SCS
params.verbose = verbose;
params.max_iters = maxit;
params.eps = opts.abstol;

% Create storage for 2 nA beta vectors, function values and timings
betaTREXMat = zeros(p, 2*nA,nC);
funTREXMat = zeros(2*nA, nC);
runTimes = zeros(2*nA, nC);

% Build list of active variables (+-)
pActInds = sort([2*(actInds-1)+1,2*actInds]);

parfor i = 1:2*nA
    
    ii = pActInds(i);
    
    j=round(ii/2);
    
    if verbose && (mod(j, 10) == 0)
        disp(j)
    end
    
    % Set sign of the jth variable
    s = (-1)^i;
    
    %   Remove data struct for parfor
    A=[AA;GG];
    b=[bb_scs;hh_scs];
    c2=cc;
    
    % Dummy initialization for parallel code
    xx = [];
    yy = [];
    ss = [];
    
    % Temporary variables for storing
    tempRunTime = zeros(nC,1);
    tempBetaMat = zeros(p,nC);
    tempFunMat = zeros(nC,1);
    
    for k=1:nC
        c=cpath(k);
        v = s * c * Xmat(:, j);
        A(n+2*p + 1, (2 * p + 2):end) = -v;
        A(n+2*p + 2, (2 * p + 2):end) = v;
        
        % Initialized to be empty for first c (checked in wrapper code)
        x2=xx;
        y2=yy;
        s2=ss;
        
        tic 
        [xx,yy,ss,info]=scs_indirectWrapper(A,b,c2,x2,y2,s2,K,params);
        tempRunTime(k) = toc;
        
        tempBetaMat(:,k) = xx((p + 2):(2 * p + 1));
        tempFunMat(k) = info.pobj;
    end
    runTimes(i,:) = tempRunTime;
    betaTREXMat(:, i,:) = tempBetaMat;
    funTREXMat(i,:) = tempFunMat ;
end

% Find best solution for each value of c
[~, ihat] = min(funTREXMat,[],1);

betaTREX=zeros(p,nC);
for k=1:nC
    betaTREX(:,k) = betaTREXMat(:, ihat(k),k);
end

if nC==1
    betaTREX = squeeze(betaTREX);
    betaTREXMat = squeeze(betaTREXMat);
    funTREXMat = squeeze(funTREXMat);
end

% So far only options and runTimes are written out in out structure
out.opts = opts;
out.runTimes = runTimes;

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


% ---------------------------------------------------------------
% ---------------------------------------------------------------
function res = isoctave
% any hack to find out whether we are running octave
s = version;
res = 0;
if exist('fflush', 'builtin') && eval(s(1)) < 7
    res = 1;
end



