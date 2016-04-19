function [S, W] = filter(X, y, q, varargin)
% FILTER  Run the knockoff filter on a data set.
%   [S, W] = FILTER(X, y, q, ...)
%
%   This function runs the knockoff procedure from start to finish,
%   creating the knockoffs, computing the test statistics, and selecting
%   variables. It is the main entry point for the knockoff package.
%
%   Inputs:
%       X - n x p scaled predictor matrix or table (n > p)
%       y - n x 1 response vector
%       q - Target false discovery rate (FDR)
%
%   Optional Inputs:
%       'Statistic' - A handle to a function f(X, X_knockoffs, y) that
%                     computes the test statistics. By default, the lasso
%                     signed-max statistic.
%       'Knockoffs' - Method to use for creating knockoffs. See CREATE.
%       'Threshold' - Method to use for thresholding. See SELECTVARS.
%       'Normalize' - Whether to automatically normalize the columns of
%                     the predictor matrix. Default: yes.
%       'Randomize' - Whether to use randomization in the construction of
%                     knockoffs and (when n<2p) augmenting the model with
%                     extra rows. Default: no.
%
%   Ouputs:
%       S - Column vector of selected variables. Contains indices if X is
%           a matrix and variable names if X is a table.
%       W - 1 x 2p vector of test statistics computed for the variables and
%           their knockoffs
%
%   See also CREATE, SELECTVARS.

parser = inputParser;
parser.CaseSensitive = false;
if (~verLessThan('matlab', '8.2')) % R2013b or later
    parser.PartialMatching = false;
end

istable_safe = @(x) ~verLessThan('matlab', '8.2') && istable(x);
parser.addRequired('X', @(x) isnumeric(x) || istable_safe(x));
parser.addRequired('y', @isnumeric);
parser.addRequired('q', @(x) isnumeric(x) && isscalar(x));

parser.addOptional('Statistic', ...
    @knockoff.stats.lassoSignedMax, ...
    @(x) isa(x, 'function_handle'));
parser.addOptional('Knockoffs', []);
parser.addOptional('Threshold', []);
parser.addOptional('Normalize', true, @islogical);
parser.addOptional('Randomize', false, @islogical);
parser.parse(X, y, q, varargin{:});

% Extract variable names.
if (istable_safe(X))
    Xnames = X.Properties.VariableNames;
    X = table2array(X);
else
    Xnames = {};
end

% Verify dimensions.
[n,p] = size(X);
if (n <= p)
    error('knockoff:DimensionError', 'Data matrix must have n > p')
elseif (n < 2*p)
    warning('knockoff:DimensionWarning', ...
        'Data matrix has p < n < 2*p. Augmenting the model with extra rows.');
    [U,~,~] = svd(X);
    U_2 = U(:,(p+1):n);
    sigma = sqrt(mean((U_2'*y).^2)); % = sqrt(RSS/(n-p))
    if (parser.Results.Randomize)
        y_extra = randn(2*p-n,1)*sigma;
    else
        seed = rng(0);
        y_extra = randn(2*p-n,1)*sigma;
        rng(seed);
    end
    y = [y;y_extra];
    X = [X;zeros(2*p-n,p)];
end

% Normalize X, if necessary.
if (parser.Results.Normalize)
    X = knockoff.private.normc(X);
end

% Run the knockoff filter.
X_ko = knockoff.create(X, parser.Results.Knockoffs, parser.Results.Randomize);
W = parser.Results.Statistic(X, X_ko, y);
S = knockoff.selectVars(W, q, parser.Results.Threshold);

if (~isempty(Xnames))
    S = Xnames(S);
end

end