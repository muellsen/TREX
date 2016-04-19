classdef CreateEquicorrelatedTest < knockoff.tests.KnockoffTestCase
    
    methods (Test)
        function testCovariances(self)
            X = knockoff.private.normc(randn(20,10));
            X_ko = knockoff.private.createEquicorrelated(X, false);
            self.verifyCovariances(X, X_ko);
        end

        function testRandomizedCovariances(self)
            X = knockoff.private.normc(randn(20,10));
            X_ko = knockoff.private.createEquicorrelated(X, true);
            self.verifyCovariances(X, X_ko);
        end

        function testDimensionCheck(self)
            X = knockoff.private.normc(randn(10,10));
            self.verifyError(...
                @() knockoff.private.createEquicorrelated(X), ...
                'knockoff:DimensionError')
        end
        
        function testPermutationInvariance(self)
            X = knockoff.private.normc(randn(20,10));
            I = randperm(10);
            X_ko = knockoff.private.createEquicorrelated(X, false);
            X_perm_ko = knockoff.private.createEquicorrelated(X(:,I), false);
            self.verifyAlmostEqual(X_ko(:,I), X_perm_ko)
        end
    end
    
    methods
        function verifyCovariances(self, X, X_ko)
            G = X'*X;
            s = min(2*min(eig(G)), 1);
            s = repmat(s, [1, size(X,2)]);
            self.verifyAlmostEqual(X_ko'*X_ko, G);
            self.verifyAlmostEqual(X'*X_ko, G - diag(s));
        end
    end

end