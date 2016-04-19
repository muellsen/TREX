function X_ko = create(X, method, randomize)
% CREATE  Create knockoff variables
%   X_ko = CREATE(X)
%   X_ko = CREATE(X, method)
%   X_ko = CREATE(X, method, randomize)
%
%   By default, creates equi-correleated knockoffs without randomization.
%
%   Inputs:
%       X - n x p scaled data matrix (n >= 2*p)
%       method - either 'equi' (for equi-correlated knockoffs) or 'sdp'
%                (for knockoffs optimized using semi-definite programming).
%       randomize - whether to use randomization in the construction of
%                   the knockoff variables
%   
%   Outputs:
%       X_ko - n x p knockoff variable matrix

if ~exist('method', 'var') || isempty(method), method = 'equi'; end;
if ~exist('randomize', 'var'), randomize = []; end;

method = lower(method);
switch method
    case 'equi'
        X_ko = knockoff.private.createEquicorrelated(X, randomize);
    case 'sdp'
        X_ko = knockoff.private.createSDP(X, randomize);
    otherwise
        error('Invalid knockoff creation method %s', method)
end

end