classdef CreateSDPTest < knockoff.tests.KnockoffTestCase
    
    methods (Test)
        function testCovariances(self)
            X = knockoff.private.normc(randn(20,10));
            X_ko = knockoff.private.createSDP(X, false);
            self.verifyCovariances(X, X_ko);
        end

        function testRandomizedCovariances(self)
            X = knockoff.private.normc(randn(20,10));
            X_ko = knockoff.private.createSDP(X, true);
            self.verifyCovariances(X, X_ko);
        end

        function testDimensionCheck(self)
            X = knockoff.private.normc(randn(10,10));
            self.verifyError(@() knockoff.private.createSDP(X), ...
                'knockoff:DimensionError')
        end
        
        function testPermutationInvariance(self)
            X = knockoff.private.normc(randn(20,10));
            I = randperm(10);
            X_ko = knockoff.private.createSDP(X, false);
            X_perm_ko = knockoff.private.createSDP(X(:,I), false);
            self.verifyEqual(X_ko(:,I), X_perm_ko, 'AbsTol', 1e-4);
        end
    end
    
    methods
        function verifyCovariances(self, X, X_ko)
            G = X'*X;
            self.verifyAlmostEqual(X_ko'*X_ko, G);
            self.verifyAlmostEqual(offdiag(X'*X_ko), offdiag(G))
            self.verifyLessThan(diag(X'*X_ko), 1 + 1e-5);
        end
    end

end

function B = offdiag(A)
B = A - diag(diag(A));
end