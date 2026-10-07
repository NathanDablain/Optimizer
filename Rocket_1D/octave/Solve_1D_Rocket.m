function O = Solve_1D_Rocket(O)

  option = struct('eval_cost', true, 'eval_grad', true, 'eval_eq', true,...
                  'eval_jac', true, 'eval_hes', true);

  option_step = struct('eval_cost', true, 'eval_grad', false, 'eval_eq', true,...
                  'eval_jac', false, 'eval_hes', false);

  dim = O.N_decision_variables + O.N_constraints;
  A = zeros(dim,dim);

  [cost, grad, eq, jacobian, hessian] = Evaluate_1D_Rocket(O, option);

  for iterations = 1 : O.iterations

    A(1:O.N_decision_variables,1:O.N_decision_variables) = hessian;
    A(1:O.N_decision_variables,(O.N_decision_variables+1):end) = jacobian';
    A((O.N_decision_variables+1):end,1:O.N_decision_variables) = jacobian;

    b = [-grad;-eq];
    x = A \ b;

    delta_z = x(1:O.N_decision_variables,1);
    lambda_new = x(O.N_decision_variables+1:end,1);
    delta_lambda = lambda_new - O.lambda;

    alpha = 1.0;

    % Loop through slack variables, find which will become the most negative, modify
    % delta_z so that it does not become negative
    tau = 0.995;
    z_new = O.z + delta_z;

    for j = 1:O.N_knots
      z_start = O.knot_size*(j-1) + 1;
      z_end = O.knot_size*j;
      z_knot = O.z(z_start:z_end);
      z_knot_new = z_new(z_start:z_end);
      delta_z_knot = delta_z(z_start:z_end);

      s_new = [z_knot_new(3) - O.lb(1);...
               z_knot_new(4) - O.lb(2);...
               O.ub(1) - z_knot_new(4)];
      s_cur = [z_knot(3) - O.lb(1);...
               z_knot(4) - O.lb(2);...
               O.ub(1) - z_knot(4)];
      for k = 1:length(s_new)
        if s_new(k) <= 0
          alpha_check = - s_cur(k) / (s_new(k) - s_cur(k));
          alpha = min(alpha, alpha_check);
        end
      end

    end

    alpha = alpha * tau;

    % alpha is the max step we can take that will not violate our bounds
    % now perform a ternary search from 0 -> alpha to find the lowest cost
    alow = 0;
    ahigh = alpha;
    O1 = O;
    O2 = O;
    for i = 1:12
      amid1 = alow + (1/3)*(ahigh - alow);
      amid2 = alow + (2/3)*(ahigh - alow);
      % Update mid1
      O1.z = O.z + amid1*delta_z;
      [cost1, ~, eq1, ~, ~] = Evaluate_1D_Rocket(O1, option_step);
      p1 = cost1 + max(abs(eq1));

      % Update mid2
      O2.z = O.z + amid2*delta_z;
      [cost2, ~, eq2, ~, ~] = Evaluate_1D_Rocket(O2, option_step);
      p2 = cost2 + max(abs(eq2));

      if p1 <= p2
        ahigh = amid2;
      elseif p1 > p2
        alow = amid1;
      end

    end
    alpha = (amid1 + amid2)/2;

    O.z = O.z + alpha*delta_z;
    O.lambda = O.lambda + alpha*delta_lambda;

    [cost, grad, eq, jacobian, hessian] = Evaluate_1D_Rocket(O, option);

    % Primal feasability checks if we are satisfying our constraints
    O.pf = max(abs(eq));

    % Dual feasability checks if we are at a local minimum while balancing our constraints
    % This is actually the stationary condition
    lagrange_grad = (grad + jacobian'*O.lambda);
    O.df = max(abs(lagrange_grad));

  end

  for i = 1:O.N_knots
    knot_start = O.knot_size*(i-1) + 1;
    knot_end = O.knot_size*i;
    local_knot = O.z(knot_start:knot_end);

    O.states(:,i) = local_knot(1:3);
    O.inputs(i) = local_knot(4);
  end

end

