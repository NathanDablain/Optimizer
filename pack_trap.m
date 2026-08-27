function x = pack_trap(states, inputs, slack)
  % States matrix should have dimension of (N_knots, N_states)
  % Inputs matrix should have dimension of (N_knots, N_inputs)
  % slack matrix should have dimension of (N_slack_vars, 1)

  N_knots = height(states);
  N_states = width(states);
  N_inputs = width(inputs);
  N_sv = width(slack);

  x = zeros((N_states + N_inputs)*N_knots + N_sv,1);

  knot_size = N_states + N_inputs;

  for i = 0:N_knots-1
    state_start = i * knot_size + 1;
    input_start = state_start + N_states;
    x(state_start:input_start-1) = states(i+1,:)';
    x(input_start:input_start+N_inputs-1) = inputs(i+1,:)';
  end

  x(input_start+N_inputs:end) = slack;
end

