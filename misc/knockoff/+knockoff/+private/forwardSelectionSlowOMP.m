function varargout = forwardSelectionSlowOMP(X, y)
% FORWARDSELECTIONSLOWOMP  Slow reference implementation of forward 
%   selection with orthogonal matching pursuit (OMP)

function residual = target(X, x, ~)
    warning_state = warning('off', 'stats:regress:RankDefDesignMat');
    [~,~,residual] = regress(y, [X x]);
    warning(warning_state);
end

[varargout{1:nargout}] = ...
    knockoff.private.sequentialfs(@criterion, @target, X, y);

end

function c = criterion(~, x, residual)
    c = -abs(dot(x, residual));
end