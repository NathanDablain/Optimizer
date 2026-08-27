function partial_block = evaluate_second_partial_derivative(states,inputs,lambdas,h)
  partial_block = zeros(3,3);
  % Should return a matrix of states+inputs by states+inputs
  % First row is partial of first state followed by partial of each state and input
  % Second row is partial of second state followed by partial of each state and input
  % ...
  % we are left with lambda_1*(-0.5*h*dx_1) + lambda_2*(-0.5*h*dx_2)
  % This becomes:
  % lambda_1*(-0.5*h*v) +...
  % lambda_2*(-0.5*h*u)
  % For smooth functions the off diagonal terms should equal

  % Because our system is linear, there is no need to evaluate the block
end
