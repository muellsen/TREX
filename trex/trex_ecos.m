function [betaTREX, betaTREXMat, funTREXMat,out] = trex_ecos(Xmat, Yvec, inopts)
% Tuning-free Regression that adapts to the design matrix X (TREX)
% Exact solution using ecos sdp solver and enumerate all 2p solutions
% (single-thread mode for repetitions)
% Input:  Response Yvec in R^n
%         Data Xmat in R^nxp (n p-dimensional measurements)
%         structure of inopts with free parameters
% Output: TREX solutions: betaTREX (best),betaTREXMat (2*p x p matrix sorted by dimension
%         and +/- signs),funTREXMat (2*p function values)

[n, p] = size(Xmat);

% Options for c-trex
defopts.feastol = 1e-12;
defopts.reltol = 1e-12;
defopts.abstol = 1e-12;
defopts.maxit = 200;
defopts.cpath = 0.5;
defopts.activeSet = 1:p;
defopts.verbose = 0;

% Merge options inopts and defopts
if nargin < 3 || isempty(inopts) % no input options available
    opts = defopts;
else
    opts = getoptions(inopts, defopts);
end

verbose = inopts.verbose;
cpath = inopts.cpath;

[n, p] = size(Xmat);
v = sparse(n, 1); % initialization
cc = [ones(p + 1, 1); zeros(p + n, 1)];
% positive orthant part:
Ip = speye(p); In = speye(n);
GG = [sparse(2 * p,1), [-Ip; -Ip], [Ip; -Ip], sparse(2 * p, n)];
% second-order cone part:
GG = [
    GG;
    [-1, sparse(1, 2 * p), -v'];
    [-1, sparse(1, 2 * p), v'];
    [sparse(n, 2 * p + 1), 2 * In]
    ];
hh = zeros(size(GG, 1), 1);
dims = struct('l', 2 * p, 'q', n + 2);
AA = [sparse(n, p + 1), Xmat, In];
bb = Yvec;
opts = inopts;%ecosoptimset('maxit', maxit, 'verbose', verbose);

% Create storage for 2p beta vectors and function values
betaTREXMat = zeros(p, 2 * p,length(cpath));
funTREXMat = zeros(2 * p, length(cpath));

i = 1;
for j = 1:p
    
    if verbose && mod(j, 100) == 0
        disp(j)
    end
    
    for s = [-1, 1]
        
        for k=1:length(cpath)
            c=cpath(k);
            v = s * c * Xmat(:, j);
            GG(dims.l + 1, (2 * p + 2):end) = -v;
            GG(dims.l + 2, (2 * p + 2):end) = v;
            % https://www.embotech.com/ECOS/Matlab-Interface/Matlab-Native
            [xx, ~, info, ~, ~] = ecos(cc, GG, hh, dims, AA, bb, opts);
            betaTREXMat(:, i,k) = xx((p + 2):(2 * p + 1));
            funTREXMat(i,k) = info.pcost;
        end
        i = i + 1;
        
    end
end

% Find best solution for each value of c
[~, ihat] = min(funTREXMat,[],1);

betaTREX=zeros(p,length(cpath));
for k=1:length(cpath)
    betaTREX(:,k) = betaTREXMat(:, ihat(k),k);
end

% Squeeze if only one c value has been provided
if length(cpath)==1
    betaTREX = squeeze(betaTREX);
    betaTREXMat = squeeze(betaTREXMat);
    funTREXMat = squeeze(funTREXMat);
end

% Additional output
out.opts = opts;

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

