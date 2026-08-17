function x = pack_trap(states, inputs)
  % States matrix should have dimension of (N_knots, N_states)
  % Inputs matrix should have dimension of (N_knots, N_inputs)

  N_knots = height(states);
  N_states = width(states);
  N_inputs = width(inputs);

  x = zeros((N_states + N_inputs)*N_knots,1);

  for i = 0:N_states+N_inputs-1
    if i < N_states
      x(i*N_knots+1:(i+1)*N_knots,1) = states(:,i+1);
    else
      x(i*N_knots+1:(i+1)*N_knots,1) = inputs(:,i-N_states+1);
    end
  endfor

end
