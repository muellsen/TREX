classdef FilterTest < knockoff.tests.KnockoffTestCase
    % Test the main entry point for this package.
    
    methods (Test)
        % Test whether the filter is invariant under permutations of
        % the columns of the design matrix.
        function testPermutationInvariance(self)
            n = 100; p = 50; k = 5; sigma = 0.5; q = 0.20;
            X = randn(n, p);
            beta = zeros(p, 1);
            beta(randsample(p,k)) = 1;
            y = X*beta + sigma .* randn(n,1);
            
            X = array2table(X);
            I = randperm(p);
            S = knockoff.filter(X, y, q);
            S_perm = knockoff.filter(X(:,I), y, q);
            self.verifyEqual(sort(S), sort(S_perm));
        end
    end

end