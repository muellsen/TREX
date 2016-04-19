function varargout = forwardSelectionSlow(X, y)
% FORWARDSELECTIONSLOW  Slow reference implementation of forward selection

[varargout{1:nargout}] = ...
    knockoff.private.sequentialfs(@criterion, @target, X, y);

end

function c = criterion(~, x, residual)
    c = -abs(dot(x, residual));
end

function nextResidual = target(~, x, residual)
    [~,~,nextResidual] = regress(residual, x);
end