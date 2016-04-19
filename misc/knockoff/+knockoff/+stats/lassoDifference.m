function [W, Z] = lassoDifference(X, X_ko, y, nlambda)
% LASSODIFFERENCE  The difference lasso statistic W
%   [W, Z] = LASSODIFFERENCE(X, X_ko, y)
%   [W, Z] = LASSODIFFERENCE(X, X_ko, y, nlambda)
%
%   Computes the statistic
%
%     W_j = Z_j - \tilde Z_j,
%
%   where Z_j and \tilde Z_j are the maximum values of the regularization
%   parameter lambda at which the jth variable and its knockoff,
%   respectively, enter the lasso model.
%
%   See also LASSOSIGNEDMAX.

if ~exist('nlambda', 'var'), nlambda = []; end

Z = knockoff.private.lassoMaxLambda([X X_ko], y, nlambda);

p = size(X,2);
orig = 1:p; ko = (p+1):(2*p);
W = Z(orig) - Z(ko);

end