% Generic function to the trapezoidal method, utilizes dynamics
function [ineq, eq] = constraints_trap(X, t, N_states, N_inputs, IC)
  ineq = [];
  N_knots = length(t);
  [states, inputs] = unpack_trap(X, N_knots, N_states, N_inputs);
  eq = zeros(N_states*(N_knots-1), 1);

  counter = 0;
  for i = 1:N_knots-1
    h = t(i+1) - t(i);
    f1 = diffeq(states(i,:),inputs(i,:));
    f2 = diffeq(states(i+1,:),inputs(i+1,:));
    for j = 1:N_states
      counter = counter + 1;
      eq(counter) = states(i+1,j) - states(i,j) - 0.5*h*(f1(j) + f2(j));
    endfor
  endfor
  for i = 1:N_states
    eq(end+1) = IC(i) - states(1,i);
  endfor
  tolerance = 1.0e-10;
  eq(abs(eq) < tolerance) = 0.0;
end
