clear
clc
close all

O = struct(...
  'N_states', 2,...
  'N_inputs', 1,...
  'N_knots', 30,...
  'N_lb', 1,...
  'N_ub', 1,...
  'ic', [0 0.0],...
  'xd', [2 0.4],...
  'lb', -2.0,...
  'ub', 2.0,...
  'tf', 2.0,...
  't', [],...
  'knot_size', [],... % derived
  'N_decision_variables', [],... % derived
  'First_knot_constraints', 6,...
  'Middle_knot_constraints', 4,...
  'End_knot_constraints', 2,...
  'N_constraints', [],... % derived
  'N_nz', [],...
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

O.N_nz = 2*(14 + 12*(O.N_knots-2) + 4) + (2 + 2*(O.N_knots-2) + 4);

% Initialize the states and slack variables
O.z(1:length(O.ic)) = O.ic';

for i = 1:O.N_knots
  knot_start = O.knot_size*(i-1) + 1;
  knot_end = O.knot_size*i;
  O.z(knot_start+3) = O.z(knot_start+2) - O.lb(1);
  O.z(knot_start+4) = O.ub(1) - O.z(knot_start+2);
end

option = struct('eval_cost', true, 'eval_grad', true, 'eval_eq', true,...
                'eval_jac', true, 'eval_hes', true);

option_step = struct('eval_cost', true, 'eval_grad', false, 'eval_eq', false,...
                'eval_jac', false, 'eval_hes', false);
tic
[cost, grad, eq, jacobian, hessian] = Evaluate_1D_Block(O, option);

pf_tolerance = 1.0e-5;
df_tolerance = 1.0e-2;
max_iterations = 10;

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
    if z_knot_new(4) <= 0 && delta_z_knot(4) < 0
      alpha_check = - z_knot(4) / delta_z_knot(4);
      alpha = min(alpha, alpha_check);
    end
    if z_knot_new(5) <= 0 && delta_z_knot(5) < 0
      alpha_check = - z_knot(5) / delta_z_knot(5);
      alpha = min(alpha, alpha_check);
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
    [cost1, ~, ~, ~, ~] = Evaluate_1D_Block(O1, option_step);

    % Update mid2
    O2.z = O.z + amid2*delta_z;
    [cost2, ~, ~, ~, ~] = Evaluate_1D_Block(O2, option_step);

    if cost1 <= cost2
      ahigh = amid2;
    elseif cost1 > cost2
      alow = amid1;
    end

  end
  alpha = (amid1 + amid2)/2;

  O.z = O.z + alpha*delta_z;
  O.lambda = O.lambda + alpha*delta_lambda;

  [cost, grad, eq, jacobian, hessian] = Evaluate_1D_Block(O, option);

  % Primal feasability checks if we are satisfying our constraints
  [pf, pfindex] = max(abs(eq));
  if pf < pf_tolerance
    pf_satisfied = true;
  else
    pf_satisfied = false;
  end

  % Dual feasability checks if we are at a local minimum while balancing our constraints
  % This is actually the stationary condition
  lagrange_grad = (grad - jacobian'*O.lambda);
  [df, dfindex] = max(abs(lagrange_grad));
  if df < df_tolerance
    df_satisfied = true;
  else
    df_satisfied = false;
  end

  if pf_satisfied && df_satisfied
    disp(['Converged in ' num2str(iterations) ' iterations'])
    disp(['with pf of ' num2str(pf) ' and df of ' num2str(df)])
    break;
  elseif iterations == max_iterations
    disp(['Failed to converge with pf of ' num2str(pf) ' and df of ' num2str(df)])
  end
end
toc

states = zeros(2,O.N_knots);
inputs = zeros(1,O.N_knots);

for i = 1:O.N_knots
  knot_start = O.knot_size*(i-1) + 1;
  knot_end = O.knot_size*i;
  local_knot = O.z(knot_start:knot_end);

  states(:,i) = local_knot(1:2);
  inputs(i) = local_knot(3);
end

figure()
subplot(3,1,1)
plot(O.t, states(1,:))
xlabel('Time (s)')
ylabel('Position (m)')
subplot(3,1,2)
plot(O.t, states(2,:))
xlabel('Time (s)')
ylabel('Velocity (m/s)')
subplot(3,1,3)
plot(O.t, inputs)
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
