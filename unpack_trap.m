function [states, inputs, slack] = unpack_trap(x, N_knots, N_states, N_inputs)
  % x is packed to have all knot decision variables in order with
  % states sitting on top of inputs

  states = zeros(N_knots, N_states);
  inputs = zeros(N_knots, N_inputs);

  knot_size = N_states + N_inputs;
  for i = 0:N_knots-1
    state_start = i * knot_size + 1;
    input_start = state_start + N_states;
    states(i+1,:) = x(state_start:input_start-1)';
    inputs(i+1,:) = x(input_start:input_start+N_inputs-1)';
  end

  slack = x(input_start+N_inputs:end);
end
