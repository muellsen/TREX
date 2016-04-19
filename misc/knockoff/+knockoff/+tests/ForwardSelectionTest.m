classdef ForwardSelectionTest < knockoff.tests.KnockoffTestCase
    
    methods (Test)
        function testSequentialFS(self)
            X = [ 1 2 3 4 5 ];
            y = 10;
            added = knockoff.private.sequentialfs(...
                @(~,x,t) abs(x-t), @(~,~,t) t-1, X, y);
            self.verifyEqual(added, [ 5 4 3 2 1]);
            
            [added, history] = knockoff.private.sequentialfs(...
                @(~,x,t) (x-t)^2, @(~,~,t) t-1, X, y);
            self.verifyEqual(added, [ 5 4 3 2 1]);
            self.verifyEqual(history.Crit, [ 25 25 25 25 25 ]);
            self.verifyEqual(history.Target, [ 10 9 8 7 6 ]');
            self.verifyEqual(history.In, logical(...
                [ 0 0 0 0 1;
                  0 0 0 1 1;
                  0 0 1 1 1;
                  0 1 1 1 1;
                  1 1 1 1 1; ]));
        end
        
        function testForwardSelectionSlow(self)
            X = [ 1 0 0;
                  1 0 0;
                  1 0 0;
                  0 1 0;
                  0 0 1;
                  0 0 1; ];
            y = [ 1 1 1 1 1 1]';
            added = knockoff.private.forwardSelectionSlow(X, y);
            self.verifyEqual(added, [ 1 3 2 ]);
        end
        
        function testForwardSelectionSlowOMP(self)
            X = [ 1 0 0;
                  1 0 0;
                  1 0 0;
                  0 1 0;
                  0 0 1;
                  0 0 1; ];
            y = [ 1 1 1 1 1 1]';
            added = knockoff.private.forwardSelectionSlowOMP(X, y);
            self.verifyEqual(added, [ 1 3 2 ]);
        end
        
        function testForwardSelectionFast(self)
            % Check the optimized implementation by comparing it to the
            % slow reference implementation.
            n = 200; p = 100;
            X = knockoff.private.normc(randn(n, p));
            y = X * randn(p,1) + 0.1 .* randn(n,1);
            self.verifyResiduals(...
                @knockoff.private.forwardSelection, ...
                @knockoff.private.forwardSelectionSlow, ...
                X, y);
        end
        
        function testForwardSelectionFastOmp(self)
            % Check the optimized implementation by comparing it to the
            % slow reference implementation.
            n = 200; p = 100;
            X = knockoff.private.normc(randn(n, p));
            y = X * randn(p,1) + 0.1 .* randn(n,1);
            self.verifyResiduals(...
                @knockoff.private.forwardSelectionOMP, ...
                @knockoff.private.forwardSelectionSlowOMP, ...
                X, y);
         end
    end
    
    methods
        function verifyResiduals(self, actual_fn, expected_fn, varargin)
            [~,historyActual] = actual_fn(varargin{:});
            [~,historyExpected] = expected_fn(varargin{:});
            self.verifyAlmostEqual(historyActual.Target', ...
                                   historyExpected.Target');
        end
    end

end