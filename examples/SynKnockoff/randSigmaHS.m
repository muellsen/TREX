function sigVec = randSigmaHS(n,sigma1,sigma2)

sigVec =  sigma1.*ones(n,1);
randInds = randi(2,n,1)-1;
sigVec(randInds==1) = sigma2;