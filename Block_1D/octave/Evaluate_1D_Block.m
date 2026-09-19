function [cost, grad, eq, jacobian, hessian, sparse_data] = Evaluate_1D_Block(O, option)
% Initialize outputs
cost = [];
grad = [];
eq = [];
jacobian = [];
hessian = [];
sparse_data = struct('rows', [], 'cols', [], 'data', []);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                       Evaluate cost                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  mu = 1.0e-7;
  function contribution = get_first_knot_cost(z_knot)
    contribution = -mu*(log(z_knot(4)) + log(z_knot(5)));
  end

  function contribution = get_middle_knot_cost(z_knot)
    contribution = -mu*(log(z_knot(4)) + log(z_knot(5)));
  end

  function contribution = get_end_knot_cost(z_knot, xd)
    contribution = (z_knot(1) - xd(1))^2 + (z_knot(2) - xd(2))^2 - mu*(log(z_knot(4)) + log(z_knot(5)));
  end

if option.eval_cost
  cost = 0;
  for i = 1:O.N_knots
    z_start = O.knot_size*(i-1) + 1;
    z_end = O.knot_size*i;
    z_knot = O.z(z_start:z_end);
    if i == 1
      knot_cost_contribution = get_first_knot_cost(z_knot);
    elseif i == O.N_knots
      knot_cost_contribution = get_end_knot_cost(z_knot, O.xd);
    else
      knot_cost_contribution = get_middle_knot_cost(z_knot);
    end
    cost = cost + knot_cost_contribution;
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Evaluate cost gradient wrt decision variables     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function grad_vec = get_first_knot_grad(z_knot)
  grad_vec = zeros(length(z_knot), 1);
  grad_vec(4) = -mu / z_knot(4);
  grad_vec(5) = -mu / z_knot(5);
end

function grad_vec = get_middle_knot_grad(z_knot)
  grad_vec = zeros(length(z_knot), 1);
  grad_vec(4) = -mu / z_knot(4);
  grad_vec(5) = -mu / z_knot(5);
end

function grad_vec = get_end_knot_grad(z_knot, xd)
  grad_vec = zeros(length(z_knot), 1);
  grad_vec(1) = 2.0 * z_knot(1) - 2.0 * xd(1);
  grad_vec(2) = 2.0 * z_knot(2) - 2.0 * xd(2);
  grad_vec(4) = -mu / z_knot(4);
  grad_vec(5) = -mu / z_knot(5);
end

if option.eval_grad
  grad = zeros(O.N_decision_variables, 1);
  for i = 1:O.N_knots
    z_start = O.knot_size*(i-1) + 1;
    z_end = O.knot_size*i;
    z_knot = O.z(z_start:z_end);
    if i == 1
      grad_vec = get_first_knot_grad(z_knot);
    elseif i == O.N_knots
      grad_vec = get_end_knot_grad(z_knot, O.xd);
    else
      grad_vec = get_middle_knot_grad(z_knot);
    end
    grad(z_start:z_end) = grad_vec;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%               Evaluate equality constraints            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dx = get_dx(z_knot)
  dx(1) = z_knot(2);
  dx(2) = z_knot(3);
end


function eq_vec = get_first_knot_eq(z_knot, t_knot, z_next, t_next, lb, ub, ic)

  h = t_next - t_knot;
  dx1 = get_dx(z_knot(1:3));
  dx2 = get_dx(z_next(1:3));

  eq_vec = [...
            % the trapezoidal constraints for the knot
            z_next(1) - z_knot(1) - 0.5*h*(dx1(1) + dx2(1));...
            z_next(2) - z_knot(2) - 0.5*h*(dx1(2) + dx2(2));...
            % the slack variable constraints for the knot
            % first lower bound
            lb(1) - z_knot(3) + z_knot(4);...
            % then upper bound
            z_knot(3) - ub(1) + z_knot(5);...
            % the initial condition constraints for the knot
            ic(1) - z_knot(1);...
            ic(2) - z_knot(2)];
end

function eq_vec = get_middle_knot_eq(z_knot, t_knot, z_next, t_next, lb, ub)

  h = t_next - t_knot;
  dx1 = get_dx(z_knot(1:3));
  dx2 = get_dx(z_next(1:3));

  eq_vec = [...
            % the trapezoidal constraints for the knot
            z_next(1) - z_knot(1) - 0.5*h*(dx1(1) + dx2(1));...
            z_next(2) - z_knot(2) - 0.5*h*(dx1(2) + dx2(2));...
            % the slack variable constraints for the knot
            % first lower bound
            lb(1) - z_knot(3) + z_knot(4);...
            % then upper bound
            z_knot(3) - ub(1) + z_knot(5)];
end

function eq_vec = get_end_knot_eq(z_knot, t_knot, lb, ub)

  eq_vec = [...
            % the slack variable constraints for the knot
            % first lower bound
            lb(1) - z_knot(3) + z_knot(4);...
            % then upper bound
            z_knot(3) - ub(1) + z_knot(5)];
end
if option.eval_eq
  eq = zeros(O.N_constraints, 1);
  con_start = 1;
  for i = 1:O.N_knots
    z_start = O.knot_size*(i-1) + 1;
    z_end = O.knot_size*i;
    z_knot = O.z(z_start:z_end);

    if i == 1
      z_start_next = z_end + 1;
      z_end_next = z_start_next + O.knot_size - 1;
      z_next = O.z(z_start_next:z_end_next);
      eq_vec = get_first_knot_eq(z_knot, O.t(i), z_next, O.t(i+1), O.lb, O.ub, O.ic);
    elseif i == O.N_knots
      eq_vec = get_end_knot_eq(z_knot, O.t(i), O.lb, O.ub);
    else
      z_start_next = z_end + 1;
      z_end_next = z_start_next + O.knot_size - 1;
      z_next = O.z(z_start_next:z_end_next);
      eq_vec = get_middle_knot_eq(z_knot, O.t(i), z_next, O.t(i+1), O.lb, O.ub);
    end
    con_end = con_start + length(eq_vec) - 1;
    eq(con_start:con_end) = eq_vec;
    con_start = con_end + 1;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Evaluate jacobian of constraints wrt decision variables%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function jac_block = get_first_knot_jac(z_knot, t_knot, z_next, t_next)
  jac_block = zeros(6, 8);
  h = t_next - t_knot;

  jac_block(1,1) = -1;
  jac_block(1,2) = -h/2;
  jac_block(1,6) = 1;
  jac_block(1,7) = -h/2;

  jac_block(2,2) = -1;
  jac_block(2,3) = -h/2;
  jac_block(2,7) = 1;
  jac_block(2,8) = -h/2;

  jac_block(3,3) = -1;
  jac_block(3,4) = 1;

  jac_block(4,3) = 1;
  jac_block(4,5) = 1;

  jac_block(5,1) = -1;

  jac_block(6,2) = -1;

end

function jac_block = get_middle_knot_jac(z_knot, t_knot, z_next, t_next)
  jac_block = zeros(4, 8);
  h = t_next - t_knot;

  jac_block(1,1) = -1;
  jac_block(1,2) = -h/2;
  jac_block(1,6) = 1;
  jac_block(1,7) = -h/2;

  jac_block(2,2) = -1;
  jac_block(2,3) = -h/2;
  jac_block(2,7) = 1;
  jac_block(2,8) = -h/2;

  jac_block(3,3) = -1;
  jac_block(3,4) = 1;

  jac_block(4,3) = 1;
  jac_block(4,5) = 1;
end

function jac_block = get_end_knot_jac(z_knot, t_knot)
  jac_block = zeros(2, 5);

  jac_block(1,3) = -1;
  jac_block(1,4) = 1;
  
  jac_block(2,3) = 1;
  jac_block(2,5) = 1;
end

if option.eval_jac
  jacobian = zeros(O.N_constraints, O.N_decision_variables);
  con_start = 1;
  for i = 1:O.N_knots
    z_start = O.knot_size*(i-1) + 1;
    z_end = O.knot_size*i;
    z_knot = O.z(z_start:z_end);

    if i == 1
      z_start_next = z_end + 1;
      z_end_next = z_start_next + O.knot_size - 1;
      z_next = O.z(z_start_next:z_end_next);
      jac_block = get_first_knot_jac(z_knot, O.t(i), z_next, O.t(i+1));
    elseif i == O.N_knots
      jac_block = get_end_knot_jac(z_knot, O.t(i));
    else
      z_start_next = z_end + 1;
      z_end_next = z_start_next + O.knot_size - 1;
      z_next = O.z(z_start_next:z_end_next);
      jac_block = get_middle_knot_jac(z_knot, O.t(i), z_next, O.t(i+1));
    end
    con_end = con_start + height(jac_block) - 1;
    z_end = z_start + width(jac_block) - 1;
    jacobian(con_start:con_end,z_start:z_end) = jac_block;
##    if option.sparse
##      [J_rows, J_cols] =
##    end
    con_start = con_end + 1;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Evaluate hessian of lagrangian wrt decision variables  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hes_block = get_first_knot_hes(z_knot, t_knot, t_next)
  hes_block = zeros(5, 5);
  h = t_next - t_knot;

  hes_block(4,4) = mu/(z_knot(4)^2);
  hes_block(5,5) = mu/(z_knot(5)^2);
end

function hes_block = get_middle_knot_hes(z_knot, t_knot, t_next)
  hes_block = zeros(5, 5);
  h = t_next - t_knot;

  hes_block(4,4) = mu/(z_knot(4)^2);
  hes_block(5,5) = mu/(z_knot(5)^2);
end

function hes_block = get_end_knot_hes(z_knot, t_knot)
  hes_block = zeros(5, 5);

  hes_block(1,1) = 2;
  hes_block(2,2) = 2;
  hes_block(4,4) = mu/(z_knot(4)^2);
  hes_block(5,5) = mu/(z_knot(5)^2);
end

if option.eval_hes
  hessian = zeros(O.N_decision_variables, O.N_decision_variables);
  for i = 1:O.N_knots
    z_start = O.knot_size*(i-1) + 1;
    z_end = O.knot_size*i;
    z_knot = O.z(z_start:z_end);

    if i == 1
      hes_block = get_first_knot_hes(z_knot, O.t(i), O.t(i+1));
    elseif i == O.N_knots
      hes_block = get_end_knot_hes(z_knot, O.t(i));
    else
      hes_block = get_middle_knot_hes(z_knot, O.t(i), O.t(i+1));
    end

    hessian(z_start:z_end,z_start:z_end) = hes_block;
  end
end

end


