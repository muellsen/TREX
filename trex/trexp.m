function [betaTREX,betaTREXMat,funTREXMat,betaBoot,bestBootMat,betaBootMat,funBootMat,out] = trexp(Xmat,Yvec,inopts)
% Tuning-free Regression that adapts to the design matrix X (TREX)
% Heuristic solution using projected gradient optimization and restarts
% (parallel mode for repetitions)
% Input:  Response Yvec in R^n
%         Data Xmat in R^nxp (n p-dimensional measurements)
%         structure of inopts with free parameters
% Output: TREX solutions: betaTREX (best),betaTREXMat (solution in order of
%         beta0 start vectors),funTREXMat (corresponding function value)
%         Bootstrapped TREX solutions: betaBoot,bestBootMat,betaBootMat,funBootMat

global X
global Y
global normConst % constant c
global qNorm     % q approximation of inf norm

% Data dimension
[n,p]=size(Xmat);

% Options for trex
defopts.q = 40;                          % order of the q-norm to approximate the max-norm
defopts.normConst = 0.5;                 % Normalizing constant
defopts.cpath = [];                      % Regularization path
defopts.rep = 1;                         % number of repetitions
defopts.boots = 0;                       % number of bootstrap samples
defopts.beta0 = zeros(p,defopts.rep);    % Initial beta estimators (if available)
defopts.consensusThresh = 0.5;           % Majority vote threshold for B-TREX
defopts.union = [];                      % Dummy options for graph recovery

% Options for Optimizers
defopts.optTol=1e-6;
defopts.progTol=1e-6;
defopts.MaxIter = max(1000,round(0.2*p));
defopts.verbose = 0;

% Check consistency of the input options
if nargin==3
    
    if ~isfield(inopts,'beta0') && isfield(inopts,'rep')
        % Generate random sparse solutions; the rule of thumb
        % implemented here is that beta0 only contain n/2 non-zero
        % entries
        numBetas = inopts.rep;
        beta0Mat = zeros(p,numBetas);
        % First solution is the zero vector
        for i=2:numBetas
            % Generate random indices between 1 and p
            if p>n
                randInds = randperm(p,ceil(n/10));
                beta0Mat(randInds,i) = randn(ceil(n/10),1);
            else
                randInds = randperm(p,ceil(p/10));
                beta0Mat(randInds,i) = randn(ceil(p/10),1);
            end
        end
        
        defopts.beta0 = beta0Mat; % Initial values of beta values
    end
    
    if isfield(inopts,'beta0') && isfield(inopts,'rep')
        currRep = size(inopts.beta0,2);
        if currRep>=inopts.rep;
            inopts.beta0 = inopts.beta0(:,1:inopts.rep);
        else
            % Generate random sparse solutions; the rule of thumb
            % implemented here is that beta0 only contain n/2 non-zero
            % entries
            numBetas = inopts.rep-currRep;
            beta0Mat = zeros(p,numBetas);
            for i=1:numBetas
                % Generate random indices between 1 and p
                if p>n
                    randInds = randperm(p,ceil(n/10));
                    beta0Mat(randInds,i) = randn(ceil(n/10),1);
                else
                    randInds = randperm(p,ceil(p/10));
                    beta0Mat(randInds,i) = randn(ceil(p/10),1);
                end
            end
            inopts.beta0 = [inopts.beta0,beta0Mat];% Initial values of beta values
        end
    end
    
end

% Merge options inopts and defopts
if nargin < 3 || isempty(inopts) % no input options available
    opts = defopts;
else
    opts = getoptions(inopts, defopts);
end

% Current option display
if opts.verbose
    opts
    pause(0.5)
end

% Set values from options
qNorm = opts.q;
numRep = opts.rep;
numBoots = opts.boots;
consensusThresh = opts.consensusThresh;

% Normalizing constant for trex objective
normConst = opts.normConst;

% Regularization path;
cpath = opts.cpath;

if isempty(cpath)
    cpath = normConst;
end
nC = length(cpath);

% Optimizer settings
l1g2opts.optTol = opts.optTol;
l1g2opts.progTol = opts.progTol;
l1g2opts.MaxIter = opts.MaxIter;
l1g2opts.verbose = opts.verbose;

% Output results
betaTREXMat = zeros(p,numRep,nC);
funTREXMat = zeros(numRep,nC);
funEvalsTREXMat = zeros(numRep,nC);

betaBootMat = zeros(p,numRep,numBoots);
funBootMat = zeros(numRep,numBoots);
funEvalsBootMat = zeros(numRep,numBoots);

X = Xmat;
Y = Yvec;

topts.X = X;
topts.Y = Y;
topts.normConst = normConst;
topts.qNorm = qNorm;

verbose = opts.verbose;
beta0 = opts.beta0;

% Generate the classic TREX solution
parfor i=1:numRep
    
    if verbose
        % Display
        disp(['###### TREX Iteration ',num2str(i),' #######']);
    end
    
    % Vector of the penalty parameters (default 1 for trex)
    lambdaVec = ones(p,1);
    
    % Start points
    xstart = beta0(:,i);
    
    % Mark Schmidt's general L1 solver
    
    tempFunMat = zeros(1,nC);
    tempBetaMat = zeros(p,nC);
    tempFunEvalsMat = zeros(1,nC);
    
    % Loop over regularization path
    for k=1:nC
        
        % Adapt constant
        c_val = cpath(k);
        
        [tempBeta,tempF,temp2] = L1General2_PSSgbP(@objTREXp,xstart,lambdaVec,topts,c_val,l1g2opts);
        tempFunMat(1,k) = tempF;
        tempBetaMat(:,k) = tempBeta;
        tempFunEvalsMat(1,k) = temp2;
        % warm start
        xstart = tempBeta;
    end
    
    betaTREXMat(:,i,:) = tempBetaMat;
    funTREXMat(i,:) = tempFunMat;
    funEvalsTREXMat(i,:) = tempFunEvalsMat;
    
    
    %betaTREXMat(:,i) = L1General2_PSSsp(@objTREX,xstart,lambdaVec,l1g2opts);
    %betaTREXMat(:,i) = L1General2_TMP(@objTREX,xstart,lambdaVec,l1g2opts);
    
    %funTREXMat(i) = objTREXp(temp,topts) + sum(abs(temp));
    
    if verbose && nC==0
        disp(['######################################']);
        disp(['TREX Solution ',num2str(i),': Sparsity: ',num2str(sum(betaTREXMat(:,i,:)~=0)),' Objval: ',num2str(funTREXMat(i,:))]);
        disp(['######################################']);
        disp([''])
    end
    
end

% Sort TREX solutions by function value
[~,sortedInds] = sort(funTREXMat,1,'ascend');

% Return best solution
betaTREX = squeeze(betaTREXMat(:,sortedInds(1),:));

if numBoots>0
    
    % Generate the bootstrap samples
    for b=1:numBoots
        
        % Sequential bootstrap (unique selection of m = (1-1/exp(1))*n
        % indices) to maximize information content of the boostrap
        cnt=0;
        mBoot = round((1-1/exp(1))*n);
        bootInds = randi(n,1,n);
        while length(unique(bootInds))~=mBoot
            bootInds = randi(n,1,n);
            cnt=cnt+1;
        end
        
        X = Xmat(bootInds,:);
        Y = Yvec(bootInds);
        
        if opts.verbose
            % Display
            disp(['###### BOOTSTRAP Percentage ',num2str(b/numBoots*100),' #######']);
        end
        
        parfor i=1:numRep
            
            % Vector of the penalty parameters (default 1 for trex)
            lambdaVec = ones(p,1);
            
            % Start points
            xstart = opts.beta0(:,i);
            
            % Mark Schmidt's general L1 solver
            [betaBootMat(:,i,b),funEvalsBootMat(i,b)] = L1General2_PSSgb(@objTREX,xstart,lambdaVec,l1g2opts);
            %betaBootMat(:,i,b) = L1General2_TMP(@objTREX,xstart,lambdaVec,l1g2opts);
            funBootMat(i,b) = objTREX(squeeze(betaBootMat(:,i,b))) + sum(abs(squeeze(betaBootMat(:,i,b))));
            
            if opts.verbose
                disp(['######################################']);
                disp(['Bootstrap: ',num2str(b),' Repetition: ',num2str(i), ' Sparsity: ',num2str(sum(betaBootMat(:,i,b)~=0)),' Objval: ',num2str(funBootMat(i,b))]);
                disp(['######################################']);
                disp([''])
            end
        end
        temp = funBootMat(:,b);
        
        % Sort the solutions by function value
        [~,sortedInds] = sort(temp,'ascend');
        betaBootMat(:,:,b) = betaBootMat(:,sortedInds,b);
        funBootMat(:,b) = funBootMat(sortedInds,b);
        
    end
    
    % Look only at the best solution and select the support using the bootstrap
    bestBootMat = squeeze(betaBootMat(:,1,:));
    
    % Find all j where beta is non-zero in more than consensusThresh of the bootstraps
    suppInds = find(sum(bestBootMat~=0,2)>(numBoots*consensusThresh));
    betaBoot = zeros(p,1);
    %beta(suppInds) = median(bestBootMat(suppInds,:),[],2);
    %beta(suppInds) = median(bestBootMat(suppInds,:),2);
    
    for s=1:length(suppInds)
        nnzInds = find(bestBootMat(suppInds(s),:)~=0);
        betaBoot(suppInds(s)) = median(bestBootMat(suppInds(s),nnzInds));
    end
    
    out.funEvalsBootMat = funEvalsBootMat;
else
    betaBootMat = [];
    funBootMat = [];
    bestBootMat = [];
    betaBoot = [];
end


out.opts = opts;
out.funEvalsTREXMat = funEvalsTREXMat;

% % Epsilon for calling an entry a true 0;
% epsSp = 1e-9;
%
% % Majority vote for the repetitions
% betaInds = find(sum(abs(betaMat)>epsSp,2)>majorThresh*numRep);
%
% betaMaj = zeros(p,1);
% %betaMaj(betaInds) = sum(betaMat(betaInds,:),2)./sum(betaMat(betaInds,:)~=0,2);
% betaMaj(betaInds) = median(betaMat(betaInds,:),2);
%
%
% % Return the best beta in terms of objective function
% beta = betaMat(:,1);



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
