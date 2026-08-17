function J = cost_step2(X, N)
  J = 0;
  for i = 1:N
    J = J + X(2*N+i)^2;
  endfor

end
