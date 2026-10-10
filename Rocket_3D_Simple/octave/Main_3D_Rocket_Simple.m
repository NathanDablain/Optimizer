clear
##clc
close all

O = struct(...
  'N_states', 10,...
  'N_inputs', 2,...
  'N_knots', 30,...
  'N_lb', 2,...
  'N_ub', 3,...
  'ic', [0 0 -100 50 0 0 0 0 0 0 0.1 pi],...
  'xd', [-3000 400 -5000],...
  'lb', [0.0 (0.0 - pi/12)],...
  'ub', [0 1.0 (2*pi + pi/12)],...
  'tf', 19.0,...
  't', [],...
  'knot_size', [],... % derived
  'N_decision_variables', [],... % derived
  'First_knot_constraints', 18,...
  'Middle_knot_constraints', 10,...
  'End_knot_constraints', 0,...
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

gam_ic = pi/4;
chi_ic = 0.0;
O.ic(5:8) = [cos(gam_ic/2)*cos(chi_ic/2),...
            -sin(gam_ic/2)*sin(chi_ic/2),...
             sin(gam_ic/2)*cos(chi_ic/2),...
             cos(gam_ic/2)*sin(chi_ic/2)];
% Initialize the states and slack variables
O.z(1:length(O.ic)) = O.ic';

% Perform a forward simulation
O = Simulate_Rocket(O);

option = struct('eval_cost', true, 'eval_grad', true, 'eval_eq', true,...
                'eval_jac', true, 'eval_hes', true, 'full_hes', true);

option_step = struct('eval_cost', true, 'eval_grad', false, 'eval_eq', false,...
                'eval_jac', false, 'eval_hes', false, 'full_hes', false);
tic
[cost, grad, eq, jacobian, hessian] = Evaluate_3D_Rocket_Simple(O, option);

  pf_tolerance = 1.0e-2;
  df_tolerance = 1.0e-3;
  max_iterations = 50;

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

       s_new = [z_knot_new(11) - O.lb(1);...
               z_knot_new(12) - O.lb(2);...
               O.ub(1) - z_knot_new(3);...
               O.ub(2) - z_knot_new(11);...
               O.ub(3) - z_knot_new(12)];
      s_cur = [z_knot(11) - O.lb(1);...
               z_knot(12) - O.lb(2);...
               O.ub(1) - z_knot(3);...
               O.ub(2) - z_knot(11);...
               O.ub(3) - z_knot(12)];

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
      [cost1, ~, ~, ~, ~] = Evaluate_3D_Rocket_Simple(O1, option_step);
      p1 = cost1;

      % Update mid2
      O2.z = O.z + amid2*delta_z;
      [cost2, ~, ~, ~, ~] = Evaluate_3D_Rocket_Simple(O2, option_step);
      p2 = cost2;

      if p1 <= p2
        ahigh = amid2;
      elseif p1 > p2
        alow = amid1;
      end

    end
    alpha = (amid1 + amid2)/2;

    O.z = O.z + alpha*delta_z;
    O.lambda = O.lambda + alpha*delta_lambda;

    [cost, grad, eq, jacobian, hessian] = Evaluate_3D_Rocket_Simple(O, option);

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

q_norm = zeros(1,O.N_knots);
states = zeros(10,O.N_knots);
inputs = zeros(2,O.N_knots);

for i = 1:O.N_knots
  knot_start = O.knot_size*(i-1) + 1;
  knot_end = O.knot_size*i;
  local_knot = O.z(knot_start:knot_end);

  states(:,i) = local_knot(1:10);
  inputs(:,i) = local_knot(11:12);
  q_norm(i) = norm(local_knot(5:8));
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
plot(O.t, inputs(1,:))
subplot(3,1,3)
plot(O.t, inputs(2,:))

figure()
subplot(2,1,1)
plot(O.t, q_norm)
subplot(2,1,2)
plot3(states(1,:), states(2,:), -states(3,:))
xlabel('North (m)')
ylabel('East (m)')
zlabel('Up (m)')
##daspect([1 0.1 1])
