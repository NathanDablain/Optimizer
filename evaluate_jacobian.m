function Jac = evaluate_jacobian(z, cons)
  % Evaluate the jacobian of each constraint wrt all decision variables
  % There are really three categories of constraints in trapezoidal collocation
  % There are the trapezoidal constraints that bind one knot to the next, there
  % are N_knots - 1 of these and they only depend on the decision variables of the
  % knot in question and the next knot.
  % Then there are the initial conditions, there is one for each state and they depend
  % on just the first knot.
  % Then there are the slack variables, they depend on the slack variable itself and the
  % state or input being bound

  % Number of constraints
  N = cons.N_states*(cons.N_knots - 1) + length(cons.IC) + cons.N_sv;

  % Number of decision variables
  M = length(z);

  % Initialize jacobian
  Jac = zeros(N,M);

  knot_size = (M - cons.N_sv) / cons.N_knots;
  [states, inputs, slack] = unpack_trap(z, cons.N_knots, cons.N_states, cons.N_inputs);

  % Load the jacobian with the type 1 constraints
  for i = 1 : cons.N_knots-1
    h = cons.t(i+1) - cons.t(i);
    starting_index_row = cons.N_states*(i - 1) + 1;
    ending_index_row = starting_index_row + cons.N_states - 1;
    starting_index_col = knot_size*(i - 1) + 1;
    ending_index_col = starting_index_col + 2*knot_size - 1;
    Jac(starting_index_row:ending_index_row, starting_index_col:ending_index_col) = ...
        evaluate_partial_derivative(states(i:i+1,:), inputs(i:i+1,:), h);
  end

  % Load the jacobian with the type 2 constraints
  for i = 1 : cons.N_states
    con_index = cons.N_states*(cons.N_knots - 1) + i;
    Jac(con_index,i) = -1;
  end

  % Load the jacobian with the type 3 constraints
  for i = 1 : length(cons.lb_index)
    con_index = cons.N_states*cons.N_knots + i;
    var_index = M - cons.N_sv + i;
    Jac(con_index, cons.lb_index(i)) = -1.0;
    Jac(con_index, var_index) = 1.0;
  end

  for i = 1 : length(cons.ub_index)
    con_index = cons.N_states*cons.N_knots + length(cons.lb_index) + i;
    var_index = M - cons.N_sv + length(cons.lb_index) + i;
    Jac(con_index, cons.ub_index(i)) = 1.0;
    Jac(con_index, var_index) = 1.0;
  end

end

