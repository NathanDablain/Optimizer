function J = cost_step1(X, N)
    J = abs(X(N) - 2) + abs(X(N*2) - 0);
end
