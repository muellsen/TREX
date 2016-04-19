function [W, Z] = forwardSelectionOMP(X, X_ko, y)
% FORWARDSELECTIONOMP  Forward selection statistic with OMP
%   [W, Z] = FORWARDSELECTIONOMP(X, X_ko, y)
%
%   This variant of forward selection uses orthogonal matching pursuit
%   (OMP).
%
%   See also FORWARDSELECTION.

added = knockoff.private.forwardSelectionOMP([X X_ko], y);
[~,order] = sort(added);

p = size(X,2);
Z = 2*p + 1 - order;
orig = 1:p; ko = (p+1):(2*p);
W = max(Z(orig), Z(ko)) .* sign(Z(orig) - Z(ko));

end