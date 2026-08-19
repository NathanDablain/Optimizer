function [J, grad] = cost_step1_kkt(z, N)
   %(z - 2)^2 -> z^2 - 4*z + 4

   J = (z(N) - 2)^2 + (z(N*2) - 0)^2;
   grad = zeros(length(z),1);
   grad(N) = 2*z(N) - 4;
   grad(2*N) = 2*z(2*N);
end

