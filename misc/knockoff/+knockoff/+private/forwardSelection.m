function varargout = forwardSelection(X, y)
% FORWARDSELECTION  Fast implementation of forward selection
%
%   Assumes that the columns of X are normalized to 1.

[varargout{1:nargout}] = ...
    knockoff.private.sequentialfs(@criterion, @target, X, y);

end

function c = criterion(~, x, residual)
    c = -abs(dot(x, residual));
end

function nextResidual = target(~, x, residual)
    nextResidual = residual - dot(x, residual) .* x;
end