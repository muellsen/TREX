function [funVal,grad] = objTREX(betaVec)
% TREX objective function
% funVal: function value at current betaVec
% grad:   Gradient at betaVec

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

% Gradient vector
grad1 = - 2*XK/(u*denom);
grad2 = (sK2 * (X'*X)*alphaVec)/(u*denom^(qNorm+1)); 

grad = grad1 + grad2;
