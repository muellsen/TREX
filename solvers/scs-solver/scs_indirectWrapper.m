function [x, y, s, info] = scs_indirectWrapper(A,b,c,x,y,s, cone, params)
% Wrapper for the parfor implementation of scs

data.A = A;
data.b = b;
data.c = c;

% Warm start solutions used over c-path
if ~isempty(x)
    data.x = x;
    data.y = y;
    data.s = s;
end

[x, y, s, info] = scs_indirect(data, cone, params);

