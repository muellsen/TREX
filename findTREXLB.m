function c_min = findTREXLB(X,Y,betaVec,epsVec,sig)
% Function that determines the lower bound on the TREX c-path (if all
% quantitites were known
% cf. Assumption 2... Bien et al, 2018
% max(abs(X'*X*beta)) > (1+2/c)*max(abs(X'*eps))
LH = max(abs(X'*X*betaVec));
RH = max(abs(X'.*sig*epsVec));

c_min = 2/((LH/RH)-1);






