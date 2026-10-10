clear
clc
close all

O = struct(...
  'N_states', 6,...
  'N_inputs', 4,...
  'N_knots', 10,...
  'N_lb', 1,...
  'N_ub', 2,...
  'ic', [0 0 -100 25 0 -25 0 1 0 0.01],...
  'xd', [100 1000 -300],...
  'lb', [0.0],...
  'ub', [0.0 1.5],...
  'tf', 10.0,...
  't', [],...
  'knot_size', [],... % derived
  'N_decision_variables', [],... % derived
  'First_knot_constraints', 14,...
  'Middle_knot_constraints', 8,...
  'End_knot_constraints', 2,...
  'N_constraints', [],... % derived
  'z', [],...
  'lambda', []...
  );


% Set the derived quantities
O.knot_size = O.N_states +...
              O.N_inputs;

O.N_decision_variables = O.knot_size * O.N_knots;

O.N_constraints = O.First_knot_constraints +...
                  O.Middle_knot_constraints*(O.N_knots - 2) +...
                  O.End_knot_constraints;

O.z = zeros(O.N_decision_variables, 1);
O.lambda = zeros(O.N_constraints, 1);
O.t = linspace(0, O.tf, O.N_knots);

% Initialize the states and slack variables
O.z(1:length(O.ic)) = O.ic';

% Perform a forward simulation
O = Simulate_Rocket(O);

option = struct('eval_cost', true, 'eval_grad', true, 'eval_eq', true,...
                'eval_jac', true, 'eval_hes', true, 'full_hes', true);

option_step = struct('eval_cost', true, 'eval_grad', false, 'eval_eq', true,...
                'eval_jac', false, 'eval_hes', false, 'full_hes', false);
tic
[cost, grad, eq, jacobian, hessian] = Evaluate_3D_Rocket(O, option);

  pf_tolerance = 1.0e-2;
  df_tolerance = 1.0e-3;
  max_iterations = 100;

  dim = O.N_decision_variables + O.N_constraints;
  A = zeros(dim,dim);

  for iterations = 1 : max_iterations

    A(1:O.N_decision_variables,1:O.N_decision_variables) = hessian;
    A(1:O.N_decision_variables,(O.N_decision_variables+1):end) = jacobian';
    A((O.N_decision_variables+1):end,1:O.N_decision_variables) = jacobian;
    b = [-grad;-eq];
    x = A \ b;

    delta_z = x(1:O.N_decision_variables,1);
    lambda_new = x(O.N_decision_variables+1:end,1);
    delta_lambda = lambda_new - O.lambda;


    alpha = 1.0;

    % Loop bouunds and calculate slack , find which will become the most negative,
    % modify delta_z so that it does not become negative
    tau = 0.995;
    z_new = O.z + delta_z;

    for j = 1:O.N_knots
      z_start = O.knot_size*(j-1) + 1;
      z_end = O.knot_size*j;
      z_knot = O.z(z_start:z_end);
      z_knot_new = z_new(z_start:z_end);
      delta_z_knot = delta_z(z_start:z_end);

      s_new = [z_knot_new(10) - O.lb(1);...
               O.ub(1) - z_knot_new(3);...
               O.ub(2) - z_knot_new(10);...
               z_knot_new(7:9) + 1.1*ones(3,1);...
               1.1*ones(3,1) - z_knot_new(7:9)];
      s_cur = [z_knot(10) - O.lb(1);...
               O.ub(1) - z_knot(3);...
               O.ub(2) - z_knot(10);...
               z_knot(7:9) + 1.1*ones(3,1);...
               1.1*ones(3,1) - z_knot(7:9)];

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
      [cost1, ~, eq1, ~, ~] = Evaluate_3D_Rocket(O1, option_step);
      p1 = cost1 + max(abs(eq1));

      % Update mid2
      O2.z = O.z + amid2*delta_z;
      [cost2, ~, eq2, ~, ~] = Evaluate_3D_Rocket(O2, option_step);
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

    [cost, grad, eq, jacobian, hessian] = Evaluate_3D_Rocket(O, option);

    % Primal feasability checks if we are satisfying our constraints
    [pf, pfindex] = max(abs(eq));
    if pf < pf_tolerance
      pf_satisfied = true;
    else
      pf_satisfied = false;
    end

    % Dual feasability checks if we are at a local minimum while balancing our constraints
    % This is actually the stationary condition
    lagrange_grad = (grad + jacobian'*O.lambda);
    [df, dfindex] = max(abs(lagrange_grad));
    if df < df_tolerance
      df_satisfied = true;
    else
      df_satisfied = false;
    end

    if pf_satisfied && df_satisfied
      disp(['Converged in ' num2str(iterations) ' iterations'])
      disp(['pf of ' num2str(pf) ' and df of ' num2str(df)])
      break;
    elseif iterations == max_iterations
      disp(['Failed to converge with pf of ' num2str(pf) ' and df of ' num2str(df)])
    end
  end

toc

a_norm = zeros(1,O.N_knots);
states = zeros(6,O.N_knots);
inputs = zeros(4,O.N_knots);

for i = 1:O.N_knots
  knot_start = O.knot_size*(i-1) + 1;
  knot_end = O.knot_size*i;
  local_knot = O.z(knot_start:knot_end);

  states(:,i) = local_knot(1:6);
  inputs(:,i) = local_knot(7:10);
  a_norm(i) = norm(inputs(1:3,i));
end

figure()
subplot(3,1,1)
plot(O.t, states(1,:))
subplot(3,1,2)
plot(O.t, states(2,:))
subplot(3,1,3)
plot(O.t, -states(3,:))

figure()
subplot(3,1,1)
plot(O.t, states(4,:))
subplot(3,1,2)
plot(O.t, states(5,:))
subplot(3,1,3)
plot(O.t, -states(6,:))

figure()
subplot(3,1,1)
hold on
plot(O.t, inputs(1,:))
plot(O.t, inputs(2,:))
plot(O.t, inputs(3,:))
subplot(3,1,2)
plot(O.t, inputs(4,:))
subplot(3,1,3)
plot(O.t, a_norm)

figure()
plot3(states(1,:), states(2,:), -states(3,:))
