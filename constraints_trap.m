function [ineq, eq, Jac] = constraints_trap(X, cons)
  % To convert inequality constraints to equality constraints, introduce a new
  % decision variable s which represents the slack. This is the interior point method.
  % The new equality constraint becomes g(x) <= 0 -> g(x) + s = 0. f(x) -> f(x) - mu*sum(ln(s))
  % The natural log forces the slack variable to be positive
  ineq = [];
  [states, inputs, ~] = unpack_trap(X, cons.N_knots, cons.N_states, cons.N_inputs);
  % Add the constraints from collocation, from initial conditions, and from bounds
  N_colloc_constraints = cons.N_states*(cons.N_knots-1);
  N_ic_constraints = cons.N_states;
  N_bound_constraints = cons.N_sv;
  N_constraints = N_colloc_constraints + N_ic_constraints + N_bound_constraints;
  eq = zeros(N_constraints, 1);

  Jac = evaluate_jacobian(X, cons);

  counter = 0;
  for i = 1:cons.N_knots-1
    h = cons.t(i+1) - cons.t(i);
    f1 = diffeq(states(i,:)',inputs(i,:),cons.t(i));
    f2 = diffeq(states(i+1,:)',inputs(i+1,:),cons.t(i+1));
    for j = 1 : cons.N_states
      counter = counter + 1;
      eq(counter) = states(i+1,j) - states(i,j) - 0.5*h*(f1(j) + f2(j));
    end
  end

  for i = 1 : N_ic_constraints
    eq(N_colloc_constraints + i) = cons.IC(i) - states(1,i);
  end

  N_lbs = 0;
  N_ubs = 0;
  % i iterates along the length of the number of non slack decision variables
  for i = 1 : length(cons.lb)
    if cons.lb(i) ~= -inf
      N_lbs = N_lbs + 1;
      con_index = N_colloc_constraints + N_ic_constraints + N_lbs;
      var_index = length(X) - cons.N_sv + N_lbs;
      eq(con_index) = ...
          cons.lb(i) - X(i) + X(var_index);
    end
  end

  for i = 1 : length(cons.ub)
    if cons.ub(i) ~= inf
      N_ubs = N_ubs + 1;
      con_index = N_colloc_constraints + N_ic_constraints + N_lbs + N_ubs;
      var_index = length(X) - cons.N_sv + N_lbs + N_ubs;
      eq(con_index) = ...
          X(i) - cons.ub(i) + X(var_index);
    end
  end
end
