function [states, inputs] = unpack_trap(x, N_knots, N_states, N_inputs)
  % x is packed to have full vectors of states sit one on top of another
  % followed by vectors of the inputs

  states = zeros(N_knots, N_states);
  inputs = zeros(N_knots, N_inputs);

  for i = 0:N_states+N_inputs-1
    if i < N_states
      states(:,i+1) = x(i*N_knots+1:(i+1)*N_knots,1);
    else
      inputs(:,i-N_states+1) = x(i*N_knots+1:(i+1)*N_knots,1);
    endif
  endfor
end
