function varargout = forwardSelectionOMP(X, y)
% FORWARDSELECTIONOMP  Fast implementation of forward selection with
%   orthogonal mathcing pursuit (OMP)
%
%   Assumes that the columns of X are normalized to 1.

[n,p] = size(X);

Q = zeros(n,p);
i = 1;

function nextResidual = target(~, x, residual)
    % Orthonormalize using modified Gram-Schmidt.
    for j = 1:i-1
        x = x - dot(Q(:,j), x) .* Q(:,j);
    end
    q = x / norm(x);
    
    nextResidual = residual - dot(q,y) .* q;
    Q(:,i) = q;
    i = i+1;
end

[varargout{1:nargout}] = ...
    knockoff.private.sequentialfs(@criterion, @target, X, y);

end

function c = criterion(~, x, residual)
    c = -abs(dot(x, residual));
end