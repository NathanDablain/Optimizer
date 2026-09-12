clear
clc
close all

O = struct(...
  'N_states', 6,...
  'N_inputs', 2,...
  'N_knots', 10,...
  'N_lb', 3,...
  'N_ub', 2,...
  'ic', [0 0 100 50 0 pi/4 0 0],...
  'xd', [100 100 500],...
  'lb', [0 -0.5 -0.5],...
  'ub', [0.5 0.5],...
  'tf', 10.0,...
  't', [],...
  'knot_size', [],... % derived
  'N_decision_variables', [],... % derived
  'First_knot_constraints', 17,...
  'Middle_knot_constraints', 11,...
  'End_knot_constraints', 5,...
  'N_constraints', [],... % derived
  'z', [],...
  'lambda', []...
  );

O2 = struct(...
  'N_states', 8,...
  'N_inputs', 2,...
  'N_knots', 10,...
  'N_lb', 3,...
  'N_ub', 2,...
  'ic', [0 0 -100 50 0 0 0 0 0 0],...
  'xd', [100 100 500],...
  'lb', [0 -0.5 -0.5],...
  'ub', [0.5 0.5],...
  'tf', 10.0,...
  't', [],...
  'knot_size', [],... % derived
  'N_decision_variables', [],... % derived
  'First_knot_constraints', 17,...
  'Middle_knot_constraints', 11,...
  'End_knot_constraints', 5,...
  'N_constraints', [],... % derived
  'z', [],...
  'lambda', []...
  );
    % Set the derived quantities
O.knot_size = O.N_states +...
              O.N_inputs +...
              O.N_lb +...
              O.N_ub;

O.N_decision_variables = O.knot_size * O.N_knots;

O.N_constraints = O.First_knot_constraints +...
                  O.Middle_knot_constraints*(O.N_knots - 2) +...
                  O.End_knot_constraints;

O.z = zeros(O.N_decision_variables, 1);
O.lambda = zeros(O.N_constraints, 1);
O.t = linspace(0, O.tf, O.N_knots);

O2.knot_size = O2.N_states +...
              O2.N_inputs +...
              O2.N_lb +...
              O2.N_ub;

O2.N_decision_variables = O2.knot_size * O2.N_knots;

O2.z = zeros(O2.N_decision_variables, 1);

O2.t = linspace(0, O2.tf, O2.N_knots);
O2.ic(5:8) = [1*cos(O.ic(6)/2)*cos(O.ic(5)/2) + 0*sin(O.ic(6)/2)*sin(O.ic(5)/2),...
              0*cos(O.ic(6)/2)*cos(O.ic(5)/2) - 1*sin(O.ic(6)/2)*sin(O.ic(5)/2),...
              1*sin(O.ic(6)/2)*cos(O.ic(5)/2) + 0*cos(O.ic(6)/2)*sin(O.ic(5)/2),...
              1*cos(O.ic(6)/2)*sin(O.ic(5)/2) - 0*sin(O.ic(6)/2)*cos(O.ic(5)/2)];
% Initialize the states and slack variables
O.z(1:length(O.ic)) = O.ic';

% Perform a forward simulation
[O, O2] = Simulate_Rocket(O, O2);

% Initialize slack
##for i = 1:O.N_knots
##  knot_start = O.knot_size*(i-1) + 1;
##  knot_end = O.knot_size*i;
##  O.z(knot_start+8) = O.z(knot_start+2) - O.lb(1);
##  O.z(knot_start+9) = O.z(knot_start+6) - O.lb(2);
##  O.z(knot_start+10) = O.z(knot_start+7) - O.lb(3);
##  O.z(knot_start+11) = O.ub(1) - O.z(knot_start+6);
##  O.z(knot_start+12) = O.ub(2) - O.z(knot_start+7);
##end
##
##option = struct('eval_cost', true, 'eval_grad', true, 'eval_eq', true,...
##                'eval_jac', true, 'eval_hes', true);
##
##option_step = struct('eval_cost', true, 'eval_grad', false, 'eval_eq', false,...
##                'eval_jac', false, 'eval_hes', false);
##tic
##[cost, grad, eq, jacobian, hessian] = Evaluate_3D_Rocket(O, option);
##
##  pf_tolerance = 1.0e-5;
##  df_tolerance = 1.0e-4;
##  max_iterations = 20;
##
##  dim = O.N_decision_variables + O.N_constraints;
##  A = zeros(dim,dim);
##
##  for iterations = 1 : max_iterations
##
##    A(1:O.N_decision_variables,1:O.N_decision_variables) = hessian;
##    A(1:O.N_decision_variables,(O.N_decision_variables+1):end) = jacobian';
##    A((O.N_decision_variables+1):end,1:O.N_decision_variables) = jacobian;
##    b = [-grad;-eq];
##    x = A \ b;
##
##    delta_z = x(1:O.N_decision_variables,1);
##    lambda_new = x(O.N_decision_variables+1:end,1);
##    delta_lambda = lambda_new - O.lambda;
##
##    alpha = 1.0;
##
##    % Loop through slack variables, find which will become the most negative, modify
##    % delta_z so that it does not become negative
##    tau = 0.995;
##    z_new = O.z + delta_z;
##
##    for j = 1:O.N_knots
##      z_start = O.knot_size*(j-1) + 1;
##      z_end = O.knot_size*j;
##      z_knot = O.z(z_start:z_end);
##      z_knot_new = z_new(z_start:z_end);
##      delta_z_knot = delta_z(z_start:z_end);
##      if z_knot_new(9) <= 0 && delta_z_knot(9) < 0
##        alpha_check = - z_knot(9) / delta_z_knot(9);
##        alpha = min(alpha, alpha_check);
##      end
##      if z_knot_new(10) <= 0 && delta_z_knot(10) < 0
##        alpha_check = - z_knot(10) / delta_z_knot(10);
##        alpha = min(alpha, alpha_check);
##      end
##      if z_knot_new(11) <= 0 && delta_z_knot(11) < 0
##        alpha_check = - z_knot(11) / delta_z_knot(11);
##        alpha = min(alpha, alpha_check);
##      end
##      if z_knot_new(12) <= 0 && delta_z_knot(12) < 0
##        alpha_check = - z_knot(12) / delta_z_knot(12);
##        alpha = min(alpha, alpha_check);
##      end
##      if z_knot_new(13) <= 0 && delta_z_knot(13) < 0
##        alpha_check = - z_knot(13) / delta_z_knot(13);
##        alpha = min(alpha, alpha_check);
##      end
##    end
##
##    alpha = alpha * tau;
##
##    % alpha is the max step we can take that will not violate our bounds
##    % now perform a ternary search from 0 -> alpha to find the lowest cost
##    alow = 0;
##    ahigh = alpha;
##    O1 = O;
##    O2 = O;
##    for i = 1:12
##      amid1 = alow + (1/3)*(ahigh - alow);
##      amid2 = alow + (2/3)*(ahigh - alow);
##      % Update mid1
##      O1.z = O.z + amid1*delta_z;
##      [cost1, ~, ~, ~, ~] = Evaluate_3D_Rocket(O1, option_step);
##
##      % Update mid2
##      O2.z = O.z + amid2*delta_z;
##      [cost2, ~, ~, ~, ~] = Evaluate_3D_Rocket(O2, option_step);
##
##      if cost1 <= cost2
##        ahigh = amid2;
##      elseif cost1 > cost2
##        alow = amid1;
##      end
##
##    end
##    alpha = (amid1 + amid2)/2;
##
##    O.z = O.z + alpha*delta_z;
##    O.lambda = O.lambda + alpha*delta_lambda;
##
##    [cost, grad, eq, jacobian, hessian] = Evaluate_3D_Rocket(O, option);
##
##    % Primal feasability checks if we are satisfying our constraints
##    [pf, pfindex] = max(abs(eq));
##    if pf < pf_tolerance
##      pf_satisfied = true;
##    else
##      pf_satisfied = false;
##    end
##
##    % Dual feasability checks if we are at a local minimum while balancing our constraints
##    % This is actually the stationary condition
##    lagrange_grad = (grad - jacobian'*O.lambda);
##    [df, dfindex] = max(abs(lagrange_grad));
##    if df < df_tolerance
##      df_satisfied = true;
##    else
##      df_satisfied = false;
##    end
##
##    if pf_satisfied && df_satisfied
##      disp(['Converged in ' num2str(iterations) ' iterations'])
##      break;
##    elseif iterations == max_iterations
##      disp(['Failed to converge with pf of ' num2str(pf) ' and df of ' num2str(df)])
##    end
##  end
##
##toc

states = zeros(6,O.N_knots);
inputs = zeros(2,O.N_knots);
states2 = zeros(8,O2.N_knots);
inputs2 = zeros(2,O2.N_knots);

for i = 1:O.N_knots
  knot_start = O.knot_size*(i-1) + 1;
  knot_end = O.knot_size*i;
  local_knot = O.z(knot_start:knot_end);

  states(:,i) = local_knot(1:6);
  inputs(:,i) = local_knot(7:8);

  knot_start = O2.knot_size*(i-1) + 1;
  knot_end = O2.knot_size*i;
  local_knot = O2.z(knot_start:knot_end);

  states2(:,i) = local_knot(1:8);
  inputs2(:,i) = local_knot(9:10);
end

figure()
subplot(3,1,1)
plot(O.t, states(1,:))
hold on
plot(O2.t, states2(1,:), 'r--')
subplot(3,1,2)
plot(O.t, states(2,:))
hold on
plot(O2.t, states2(2,:), 'r--')
subplot(3,1,3)
plot(O.t, states(3,:))
hold on
plot(O2.t, -states2(3,:), 'r--')

##figure()
##subplot(3,1,1)
##plot(O.t, states(4,:))
##hold on
##plot(O2.t, states2(4,:), 'r--')
##subplot(3,1,2)
##plot(O.t, states(5,:))
##subplot(3,1,3)
##plot(O.t, states(6,:))
##
##figure()
##subplot(3,1,1)
##plot(O.t, inputs(1,:))
##subplot(3,1,2)
##plot(O.t, inputs(2,:))
##subplot(3,1,3)
##plot(O.t, 1./cos(states(6,:)))
