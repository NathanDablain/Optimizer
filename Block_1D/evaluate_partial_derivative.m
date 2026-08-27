 function partial_block = evaluate_partial_derivative(states, inputs, h)

  % Take partial wrt x_k, v_k, a_k, x_k+1, v_k+1, a_k+1
  % c = z_k+1 - z_k - 0.5*h*(f_k + f_k+1)
  % c_1 = x_k+1 - x_k - 0.5*h*(v_k + v_k+1)
  % c_2 = v_k+1 - v_k - 0.5*h*(a_k + a_k+1)
  partial_block = zeros(2,6);

  partial_block(1,1) = -1;
  partial_block(1,2) = -0.5*h;
  partial_block(1,3) = 0;
  partial_block(1,4) = 1;
  partial_block(1,5) = -0.5*h;
  partial_block(1,6) = 0;

  partial_block(2,1) = 0;
  partial_block(2,2) = -1;
  partial_block(2,3) = -0.5*h;
  partial_block(2,4) = 0;
  partial_block(2,5) = 1;
  partial_block(2,6) = -0.5*h;

 end