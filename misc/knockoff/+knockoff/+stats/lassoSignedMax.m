function [W, Z] = lassoSignedMax(X, X_ko, y, nlambda)
% LASSOSIGNEDMAX  The signed maximum lasso statistic W
%   [W, Z] = LASSOSIGNEDMAX(X, X_ko, y)
%   [W, Z] = LASSOSIGNEDMAX(X, X_ko, y, nlambda)
%
%   Computes the laso statistic of Equation 1.7:
%
%     W_j = max(Z_j, \tilde Z_j) * sgn(Z_j - \tilde Z_j).
%
%   Here Z_j and \tilde Z_j are the maximum valued of the regularization
%   parameter lambda at which the jth variable and its knockoff,
%   respectively, enter the lasso model.
%
%   Note that the lasso path is not computed exactly, but approximated by
%   a fine grid of lambda values. The optional parameter 'nlambda' controls
%   the number of points in this grid. The default value is 5*p, where p 
%   is the number of columns in X. If the lasso path contains closely
%   spaced knots, it may be useful to increase the value of 'nlambda'.
%   The demo 'FirstExamples' shows how to do this.
%
%   See also LASSODIFFERENCE.

if ~exist('nlambda', 'var'), nlambda = []; end

Z = knockoff.private.lassoMaxLambda([X X_ko], y, nlambda);
%pause

p = size(X,2);
orig = 1:p; ko = (p+1):(2*p);
W = max(Z(orig), Z(ko)) .* sign(Z(orig) - Z(ko));

end