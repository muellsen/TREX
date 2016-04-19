function [funVal,L1beta] = objTREX_SA(betaVec,X,Y,normConst)
% TREX standalone objective function
% funVal: function value at current betaVec
% L1beta: L1 norm of betaVec

% Global data vector X, response Y, normalizing constant normConst, and
% q-norm approximation to inf-norm

[n,p]= size(X);

temp = sum((Y-X*betaVec).^2)/(normConst*max(abs(X'*(Y-X*betaVec))));

L1beta = sum(abs(betaVec));

funVal = temp+L1beta;

