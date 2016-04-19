function ForwardSelectionBenchmark

n = 800; p = 400;
X = knockoff.private.normc(randn(n, p));
y = X * randn(p,1) + 0.1 .* randn(n,1);

D = table();
D.call = {'FS_Naive'; 'FS_Opt'; 'FS_OMP_Naive'; 'FS_OMP_Opt'};
D.time = [
    timeit(@() knockoff.private.forwardSelectionSlow(X,y));
    timeit(@() knockoff.private.forwardSelection(X,y));
    timeit(@() knockoff.private.forwardSelectionSlowOMP(X,y));
    timeit(@() knockoff.private.forwardSelectionOMP(X,y));
    ];
D.Properties.VariableUnits = {'' 's'};
display(D);

end

