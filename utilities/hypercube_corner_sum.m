function [S] = hypercube_corner_sum(n, C0)

if n >= 1
    C = C0;
    for iteration = 1: n-1
        C = C0 .* (sum(C) - (n-iteration)*C);
    end
    S = sum(C);
else
    S = 1;
end