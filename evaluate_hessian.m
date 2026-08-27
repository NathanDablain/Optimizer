function H = evaluate_hessian(z, lambda, cons)
  % Number of decision variables
  M = length(z);

  % Initialize hessian
  H = zeros(M,M);

  knot_size = (M - cons.N_sv) / cons.N_knots;
  N_inputs = knot_size - cons.N_states;
  [states, inputs] = unpack_trap(z, cons.N_knots, cons.N_states, cons.N_inputs);

  % Load the nonlinear dynamics
  for i = 1 : cons.N_knots-1
    h = cons.t(i+1) - cons.t(i);
    starting_index = (i-1)*knot_size + 1;
    ending_index = starting_index + knot_size - 1;
    starting_con_index = (i-1)*cons.N_states + 1;
    ending_con_index = starting_con_index + cons.N_states - 1;
    % Evaluate for k decision variables
    partial_block1 = evaluate_second_partial_derivative(states(i,:),inputs(i,:),lambda(starting_con_index:ending_con_index),h);
    H(starting_index:ending_index,starting_index:ending_index) = ...
        H(starting_index:ending_index,starting_index:ending_index) +  partial_block1;

    starting_index = i*knot_size + 1;
    ending_index = starting_index + knot_size - 1;
    % Evaluate for k+1 decision variables
    partial_block2 = evaluate_second_partial_derivative(states(i+1,:),inputs(i+1,:),lambda(starting_con_index:ending_con_index),h);
    H(starting_index:ending_index,starting_index:ending_index) = partial_block2;

  end

  % Load the cost weighting
  [~, ~, double_grad] = cost(z, cons);
  H = H + double_grad;
end
