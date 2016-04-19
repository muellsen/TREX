function [funVal,grad] = objTREX_wog(betaVec)
% TREX objective function without gradient
% funVal: function value at current betaVec

% Global data vector X, response Y, normalizing constant normConst, and
% q-norm approximation to inf-norm
global X
global Y
global normConst
global qNorm

[n,p]= size(X);

K = (Y-X*betaVec);
XK = X'*K;
sK2 = sum(K.^2);

u = normConst;
denom = (sum(abs(XK).^qNorm)^(1/qNorm));

% NOTE: projVec==XK
%projVec = (X'*(Y-X*betaVec))
%alphaVec = projVec.^(qNorm-1);

alphaVec = XK.^(qNorm-1);

% Function value
funVal = sK2/(u*sum(abs(XK).^qNorm)^(1/qNorm));
