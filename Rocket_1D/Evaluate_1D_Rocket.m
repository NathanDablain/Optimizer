function [cost, grad, eq, jacobian, hessian] = Evaluate_1D_Rocket(O, option)
% Initialize outputs
cost = [];
grad = [];
eq = [];
jacobian = [];
hessian = [];
% knot vars are :
% 1 - p,
% 2 - v,
% 3 - m,
% 4 - T,
% 5 - s_m_lb,
% 6 - s_T_lb,
% 7 - s_T_ub
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                       Evaluate cost                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  mu = 3.0e-2;
  function contribution = get_first_knot_cost(z_knot)
    contribution = -mu*(log(z_knot(5)) + log(z_knot(6)) + log(z_knot(7)));
  end

  function contribution = get_middle_knot_cost(z_knot)
    contribution = -mu*(log(z_knot(5)) + log(z_knot(6)) + log(z_knot(7)));
  end

  function contribution = get_end_knot_cost(z_knot, xd)
    contribution = (z_knot(1) - xd(1))^2 -...
                    mu*(log(z_knot(5)) + log(z_knot(6)) + log(z_knot(7)));
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
  grad_vec(5) = -mu / z_knot(5);
  grad_vec(6) = -mu / z_knot(6);
  grad_vec(7) = -mu / z_knot(7);
end

function grad_vec = get_middle_knot_grad(z_knot)
  grad_vec = zeros(length(z_knot), 1);
  grad_vec(5) = -mu / z_knot(5);
  grad_vec(6) = -mu / z_knot(6);
  grad_vec(7) = -mu / z_knot(7);
end

function grad_vec = get_end_knot_grad(z_knot, xd)
  grad_vec = zeros(length(z_knot), 1);
  grad_vec(1) = 2.0 * z_knot(1) - 2.0 * xd(1);
  grad_vec(5) = -mu / z_knot(5);
  grad_vec(6) = -mu / z_knot(6);
  grad_vec(7) = -mu / z_knot(7);
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
C = 3000.0;
A = 10.52;
function knot_params = get_knot_params(z_knot)
  gravity0 = 9.8065;
  Re = 6371e3;

  knot_params.g = gravity0*((Re/(Re + z_knot(1)))^2);
  knot_params.g_part_h = -2*Re*gravity0/((Re + z_knot(1))^3);
  knot_params.g_part_h2 = (6*Re*gravity0)/((Re + z_knot(1))^4);

  if z_knot(1) < 11000
    c1 = 15.04;
    c2 = 0.00649;
    c3 = 101.29;
    c4 = 273.1;
    c5 = 288.08;
    c6 = 5.256;
    c7 = 0.2869;
    Temperature = c1 - c2*z_knot(1);
    Pressure =c3*(((Temperature+c4)/c5)^c6);

    knot_params.rho_part_h = (c2*c3*((c1 - c2*z_knot(1) + c4)/c5)^c6*...
                             (1 - c6))/(c7*(c1-c2*z_knot(1)+c4)^2);

    knot_params.rho_part_h2 = (c2^2*c3*((c1 - c2*z_knot(1) + c4)/c5)^c6*...
                             (c6 - 1)*(c6 - 2))/(c7*(c1-c2*z_knot(1)+c4)^3);
  elseif z_knot(1) >= 11000 && z_knot(1) < 25000
    c1 = -56.46;
    c2 = 22.65;
    c3 = 1.73;
    c4 = 0.000157;
    c5 = 0.2869;
    c6 = 273.1;
    Temperature = c1;
    Pressure = c2*exp(c3 - c4*z_knot(1));

    knot_params.rho_part_h = (-c2*c4*exp(c3 - c4*z_knot(1)))/(c5*(c1+c6));
    knot_params.rho_part_h2 = (c2*c4^2*exp(c3 - c4*z_knot(1)))/(c5*(c1+c6));
  else
    c1 = -131.21;
    c2 = 0.00299;
    c3 = 2.488;
    c4 = 273.1;
    c5 = 216.6;
    c6 = -11.388;
    c7 = 0.2869;
    Temperature = c1 + c2*z_knot(1);
    Pressure = c3*(((Temperature+c4)/c5)^c6);

    knot_params.rho_part_h = (c2*c3*((c1 + c2*z_knot(1) + c4)/c5)^c6*...
                             (c6 - 1))/(c7*(c1+c2*z_knot(1)+c4)^2);
    knot_params.rho_part_h2 = (c2^2*c3*((c1 + c2*z_knot(1) + c4)/c5)^c6*...
                             (c6 - 1)*(c6 - 2))/(c7*(c1+c2*z_knot(1)+c4)^3);
  end
  knot_params.rho = Pressure/(0.2869*(Temperature+273.1));

  speed_of_sound = sqrt(1.4*287*(Temperature+273.1));
  Mach = z_knot(2) / speed_of_sound;
  if Mach < 0.8
    knot_params.c_D = 0.22;

    knot_params.c_D_part_s = 0;
    knot_params.c_D_part_s2 = 0;
  elseif Mach >= 0.8 && Mach < 1.2
    c1 = 0.22;
    c2 = 0.48;
    c3 = 0.8;
    knot_params.c_D = c1 + c2*sin(pi*(Mach - c3)/c3)^2;

    knot_params.c_D_part_s = (pi*c2*sin(2*pi*z_knot(2)/(speed_of_sound*c3)))/(speed_of_sound*c3);
    knot_params.c_D_part_s2 = (2*pi^2*c2*cos((2*pi*z_knot(2))/(speed_of_sound*c3)))/(speed_of_sound^2*c3^2);
  elseif Mach >= 1.2
    c1 = 0.25;
    c2 = 0.54;
    c3 = 1.2;
    knot_params.c_D = c1 + c2/(Mach^c3);

    knot_params.c_D_part_s = (-c2*c3*((z_knot(2)/speed_of_sound)^-c3))/z_knot(2);
    knot_params.c_D_part_s2 = (c2*c3*((z_knot(2)/speed_of_sound)^-c3)*(c3 + 1))/(z_knot(2)^2);
  end

  Drag = 0.5*knot_params.c_D*knot_params.rho*A*z_knot(2)*z_knot(2);

  knot_params.dx(1) = z_knot(2);
  knot_params.dx(2) = (z_knot(4) - Drag)/z_knot(3) - knot_params.g;
  knot_params.dx(3) = -z_knot(4) / C;
end

function eq_vec = get_first_knot_eq(z_knot, t_knot, knot_params,...
                                    z_next, t_next, next_params, lb, ub, ic)

  h = t_next - t_knot;

  eq_vec = [...
            % the trapezoidal constraints for the knot
            z_next(1) - z_knot(1) - 0.5*h*(knot_params.dx(1) + next_params.dx(1));...
            z_next(2) - z_knot(2) - 0.5*h*(knot_params.dx(2) + next_params.dx(2));...
            z_next(3) - z_knot(3) - 0.5*h*(knot_params.dx(3) + next_params.dx(3));...
            % the slack variable constraints for the knot
            % first lower bound
            lb(1) - z_knot(3) + z_knot(5);...
            lb(2) - z_knot(4) + z_knot(6);...
            % then upper bound
            z_knot(4) - ub(1) + z_knot(7);...
            % the initial condition constraints for the knot
            ic(1) - z_knot(1);...
            ic(2) - z_knot(2);...
            ic(3) - z_knot(3)];
end

function eq_vec = get_middle_knot_eq(z_knot, t_knot, knot_params,...
                                     z_next, t_next, next_params, lb, ub)

  h = t_next - t_knot;

  eq_vec = [...
            % the trapezoidal constraints for the knot
            z_next(1) - z_knot(1) - 0.5*h*(knot_params.dx(1) + next_params.dx(1));...
            z_next(2) - z_knot(2) - 0.5*h*(knot_params.dx(2) + next_params.dx(2));...
            z_next(3) - z_knot(3) - 0.5*h*(knot_params.dx(3) + next_params.dx(3));...
            % the slack variable constraints for the knot
            % first lower bound
            lb(1) - z_knot(3) + z_knot(5);...
            lb(2) - z_knot(4) + z_knot(6);...
            % then upper bound
            z_knot(4) - ub(1) + z_knot(7)];
end

function eq_vec = get_end_knot_eq(z_knot, t_knot, lb, ub)

  eq_vec = [...
            % the slack variable constraints for the knot
            % first lower bound
            lb(1) - z_knot(3) + z_knot(5);...
            lb(2) - z_knot(4) + z_knot(6);...
            % then upper bound
            z_knot(4) - ub(1) + z_knot(7)];
end

if option.eval_eq
  params(1:O.N_knots) = struct('dx', zeros(1,length(O.N_states)), 'g', 0, 'rho', 0, 'c_D', 0,...
                               'g_part_h', 0, 'rho_part_h', 0, 'c_D_part_s', 0,...
                               'g_part_h2', 0, 'rho_part_h2', 0, 'c_D_part_s2', 0);
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
      params(i) = get_knot_params(z_knot);
      params(i+1) = get_knot_params(z_next);
      eq_vec = get_first_knot_eq(z_knot, O.t(i), params(i), z_next, O.t(i+1), params(i+1), O.lb, O.ub, O.ic);
    elseif i == O.N_knots
      eq_vec = get_end_knot_eq(z_knot, O.t(i), O.lb, O.ub);
    else
      z_start_next = z_end + 1;
      z_end_next = z_start_next + O.knot_size - 1;
      z_next = O.z(z_start_next:z_end_next);
      params(i+1) = get_knot_params(z_next);
      eq_vec = get_middle_knot_eq(z_knot, O.t(i), params(i), z_next, O.t(i+1), params(i+1), O.lb, O.ub);
    end
    con_end = con_start + length(eq_vec) - 1;
    eq(con_start:con_end) = eq_vec;
    con_start = con_end + 1;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Evaluate jacobian of constraints wrt decision variables%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function jac_block = get_first_knot_jac(z_knot, t_knot, knot_params, z_next, t_next, next_params)
  jac_block = zeros(9, 11);
  h = t_next - t_knot;
  p1 = z_knot(1);
  v1 = z_knot(2);
  m1 = z_knot(3);
  T1 = z_knot(4);
  c_D1 = knot_params.c_D;
  rho1 = knot_params.rho;
  g_part_h1 = knot_params.g_part_h;
  rho_part_h1 = knot_params.rho_part_h;
  c_D_part_s1 = knot_params.c_D_part_s;

  p2 = z_next(1);
  v2 = z_next(2);
  m2 = z_next(3);
  T2 = z_next(4);
  c_D2 = next_params.c_D;
  rho2 = next_params.rho;
  g_part_h2 = next_params.g_part_h;
  rho_part_h2 = next_params.rho_part_h;
  c_D_part_s2 = next_params.c_D_part_s;

  jac_block(1,1) = -1;
  jac_block(1,2) = -h/2;
  jac_block(1,8) = 1;
  jac_block(1,9) = -h/2;
  jac_block(2,1) = h*(A*v1^2*c_D1*rho_part_h1 + 2*m1*g_part_h1)/...
  (4*m1);
  jac_block(2,2) = (A*h*v1*(v1*c_D_part_s1 + 2*c_D1)*rho1/4 - m1)/m1;
  jac_block(2,3) = h*(-A*v1^2*c_D1*rho1 + 2*T1)/(4*m1^2);
  jac_block(2,4) = -h/(2*m1);
  jac_block(2,8) = h*(A*v2^2*c_D2*rho_part_h2 + 2*m2*g_part_h2)/...
  (4*m2);
  jac_block(2,9) = (A*h*v2*(v2*c_D_part_s2 + 2*c_D2)*rho2/4 + m2)/m2;
  jac_block(2,10) = h*(-A*v2^2*c_D2*rho2 + 2*T2)/(4*m2^2);
  jac_block(2,11) = -h/(2*m2);
  jac_block(3,3) = -1;
  jac_block(3,4) = h/(2*C);
  jac_block(3,10) = 1;
  jac_block(3,11) = h/(2*C);
  jac_block(4,3) = -1;
  jac_block(4,5) = 1;
  jac_block(5,4) = -1;
  jac_block(5,6) = 1;
  jac_block(6,4) = 1;
  jac_block(6,7) = 1;
  jac_block(7,1) = -1;
  jac_block(8,2) = -1;
  jac_block(9,3) = -1;

end

function jac_block = get_middle_knot_jac(z_knot, t_knot, knot_params, z_next, t_next, next_params)
  jac_block = zeros(6, 11);
  h = t_next - t_knot;
  p1 = z_knot(1);
  v1 = z_knot(2);
  m1 = z_knot(3);
  T1 = z_knot(4);
  c_D1 = knot_params.c_D;
  rho1 = knot_params.rho;
  g_part_h1 = knot_params.g_part_h;
  rho_part_h1 = knot_params.rho_part_h;
  c_D_part_s1 = knot_params.c_D_part_s;

  p2 = z_next(1);
  v2 = z_next(2);
  m2 = z_next(3);
  T2 = z_next(4);
  c_D2 = next_params.c_D;
  rho2 = next_params.rho;
  g_part_h2 = next_params.g_part_h;
  rho_part_h2 = next_params.rho_part_h;
  c_D_part_s2 = next_params.c_D_part_s;

  jac_block(1,1) = -1;
  jac_block(1,2) = -h/2;
  jac_block(1,8) = 1;
  jac_block(1,9) = -h/2;
  jac_block(2,1) = h*(A*v1^2*c_D1*rho_part_h1 + 2*m1*g_part_h1)/...
  (4*m1);
  jac_block(2,2) = (A*h*v1*(v1*c_D_part_s1 + 2*c_D1)*rho1/4 - m1)/m1;
  jac_block(2,3) = h*(-A*v1^2*c_D1*rho1 + 2*T1)/(4*m1^2);
  jac_block(2,4) = -h/(2*m1);
  jac_block(2,8) = h*(A*v2^2*c_D2*rho_part_h2 + 2*m2*g_part_h2)/...
  (4*m2);
  jac_block(2,9) = (A*h*v2*(v2*c_D_part_s2 + 2*c_D2)*rho2/4 + m2)/m2;
  jac_block(2,10) = h*(-A*v2^2*c_D2*rho2 + 2*T2)/(4*m2^2);
  jac_block(2,11) = -h/(2*m2);
  jac_block(3,3) = -1;
  jac_block(3,4) = h/(2*C);
  jac_block(3,10) = 1;
  jac_block(3,11) = h/(2*C);
  jac_block(4,3) = -1;
  jac_block(4,5) = 1;
  jac_block(5,4) = -1;
  jac_block(5,6) = 1;
  jac_block(6,4) = 1;
  jac_block(6,7) = 1;

end

function jac_block = get_end_knot_jac
  jac_block = zeros(3, 7);

  jac_block(1,3) = -1;
  jac_block(1,5) = 1;
  jac_block(2,4) = -1;
  jac_block(2,6) = 1;
  jac_block(3,4) = 1;
  jac_block(3,7) = 1;
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
      jac_block = get_first_knot_jac(z_knot, O.t(i), params(i), z_next, O.t(i+1), params(i+1));
    elseif i == O.N_knots
      jac_block = get_end_knot_jac;
    else
      z_start_next = z_end + 1;
      z_end_next = z_start_next + O.knot_size - 1;
      z_next = O.z(z_start_next:z_end_next);
      jac_block = get_middle_knot_jac(z_knot, O.t(i), params(i), z_next, O.t(i+1), params(i+1));
    end
    con_end = con_start + height(jac_block) - 1;
    z_end = z_start + width(jac_block) - 1;
    jacobian(con_start:con_end,z_start:z_end) = jac_block;
    con_start = con_end + 1;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Evaluate hessian of lagrangian wrt decision variables  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hes_block = get_first_knot_hes(z_knot, t_knot, knot_params, lambda_knot, t_next)
  hes_block = zeros(7, 7);
  h = t_next - t_knot;
  p1 = z_knot(1);
  v1 = z_knot(2);
  m1 = z_knot(3);
  T1 = z_knot(4);
  s_m_lb = z_knot(5);
  s_T_lb = z_knot(6);
  s_T_ub = z_knot(7);
  c_D1 = knot_params.c_D;
  rho1 = knot_params.rho;
  g_part_h = knot_params.g_part_h;
  rho_part_h = knot_params.rho_part_h;
  c_D_part_s = knot_params.c_D_part_s;
  g_part_h2 = knot_params.g_part_h2;
  rho_part_h2 = knot_params.rho_part_h2;
  c_D_part_s2 = knot_params.c_D_part_s2;
  lambda2 = lambda_knot(2);

  hes_block(1,1) = (h*lambda2*(A*v1^2*c_D1*rho_part_h2 +...
  2*m1*g_part_h2)/(4*m1));
  hes_block(1,2) = A*h*lambda2*v1*(v1*c_D_part_s + 2*c_D1)*...
  rho_part_h/(4*m1);
  hes_block(1,3) = -A*h*lambda2*v1^2*c_D1*rho_part_h/(4*m1^2);
  hes_block(2,1) = hes_block(1,2);
  hes_block(2,2) = A*h*lambda2*(v1^2*c_D_part_s2 + 4*v1*...
  c_D_part_s + 2*c_D1)*rho1/(4*m1);
  hes_block(2,3) = -A*h*lambda2*v1*(v1*c_D_part_s + 2*c_D1)*rho1/(4*m1^2);
  hes_block(3,1) = hes_block(1,3);
  hes_block(3,2) = hes_block(2,3);
  hes_block(3,3) = h*lambda2*(A*v1^2*c_D1*rho1 - 2*T1)/(2*m1^3);
  hes_block(3,4) = h*lambda2/(2*m1^2);
  hes_block(4,3) = hes_block(3,4);
  hes_block(5,5) = mu/(s_m_lb^2);
  hes_block(6,6) = mu/(s_T_lb^2);
  hes_block(7,7) = mu/(s_T_ub^2);
end

function hes_block = get_middle_knot_hes(z_knot, t_knot, knot_params, lambda_knot, t_next, lambda_last)
  hes_block = zeros(7, 7);
  h = t_next - t_knot;
  p1 = z_knot(1);
  v1 = z_knot(2);
  m1 = z_knot(3);
  T1 = z_knot(4);
  s_m_lb = z_knot(5);
  s_T_lb = z_knot(6);
  s_T_ub = z_knot(7);
  c_D1 = knot_params.c_D;
  rho1 = knot_params.rho;
  g_part_h = knot_params.g_part_h;
  rho_part_h = knot_params.rho_part_h;
  c_D_part_s = knot_params.c_D_part_s;
  g_part_h2 = knot_params.g_part_h2;
  rho_part_h2 = knot_params.rho_part_h2;
  c_D_part_s2 = knot_params.c_D_part_s2;
  % This only works if h is constant
  lambda2 = lambda_knot(2) + lambda_last(2);

  hes_block(1,1) = (h*lambda2*(A*v1^2*c_D1*rho_part_h2 +...
  2*m1*g_part_h2)/(4*m1));
  hes_block(1,2) = A*h*lambda2*v1*(v1*c_D_part_s + 2*c_D1)*...
  rho_part_h/(4*m1);
  hes_block(1,3) = -A*h*lambda2*v1^2*c_D1*rho_part_h/(4*m1^2);
  hes_block(2,1) = hes_block(1,2);
  hes_block(2,2) = A*h*lambda2*(v1^2*c_D_part_s2 + 4*v1*...
  c_D_part_s + 2*c_D1)*rho1/(4*m1);
  hes_block(2,3) = -A*h*lambda2*v1*(v1*c_D_part_s + 2*c_D1)*rho1/(4*m1^2);
  hes_block(3,1) = hes_block(1,3);
  hes_block(3,2) = hes_block(2,3);
  hes_block(3,3) = h*lambda2*(A*v1^2*c_D1*rho1 - 2*T1)/(2*m1^3);
  hes_block(4,3) = hes_block(3,4);
  hes_block(4,3) = h*lambda2/(2*m1^2);
  hes_block(5,5) = mu/(s_m_lb^2);
  hes_block(6,6) = mu/(s_T_lb^2);
  hes_block(7,7) = mu/(s_T_ub^2);
end

function hes_block = get_end_knot_hes(z_knot, t_knot, knot_params, lambda_knot, t_next, lambda_last)
  hes_block = zeros(7, 7);
  h = t_next - t_knot;
  p1 = z_knot(1);
  v1 = z_knot(2);
  m1 = z_knot(3);
  T1 = z_knot(4);
  s_m_lb = z_knot(5);
  s_T_lb = z_knot(6);
  s_T_ub = z_knot(7);
  c_D1 = knot_params.c_D;
  rho1 = knot_params.rho;
  g_part_h = knot_params.g_part_h;
  rho_part_h = knot_params.rho_part_h;
  c_D_part_s = knot_params.c_D_part_s;
  g_part_h2 = knot_params.g_part_h2;
  rho_part_h2 = knot_params.rho_part_h2;
  c_D_part_s2 = knot_params.c_D_part_s2;
  lambda2 = lambda_last(2);

  hes_block(1,1) = 2 + (h*lambda2*(A*v1^2*c_D1*rho_part_h2 +...
  2*m1*g_part_h2)/(4*m1));
  hes_block(1,2) = A*h*lambda2*v1*(v1*c_D_part_s + 2*c_D1)*...
  rho_part_h/(4*m1);
  hes_block(1,3) = -A*h*lambda2*v1^2*c_D1*rho_part_h/(4*m1^2);
  hes_block(2,1) = hes_block(1,2);
  hes_block(2,2) = A*h*lambda2*(v1^2*c_D_part_s2 + 4*v1*...
  c_D_part_s + 2*c_D1)*rho1/(4*m1);
  hes_block(2,3) = -A*h*lambda2*v1*(v1*c_D_part_s + 2*c_D1)*rho1/(4*m1^2);
  hes_block(3,1) = hes_block(1,3);
  hes_block(3,2) = hes_block(2,3);
  hes_block(3,3) = h*lambda2*(A*v1^2*c_D1*rho1 - 2*T1)/(2*m1^3);
  hes_block(3,4) = h*lambda2/(2*m1^2);
  hes_block(4,3) = hes_block(3,4);
  hes_block(5,5) = mu/(s_m_lb^2);
  hes_block(6,6) = mu/(s_T_lb^2);
  hes_block(7,7) = mu/(s_T_ub^2);
end

if option.eval_hes
  hessian = zeros(O.N_decision_variables, O.N_decision_variables);
  con_start = 1;
  for i = 1:O.N_knots
    z_start = O.knot_size*(i-1) + 1;
    z_end = O.knot_size*i;
    z_knot = O.z(z_start:z_end);

    if i == 1
      con_end = con_start + O.First_knot_constraints - 1;
      lambda_knot = O.lambda(con_start:con_end);
      hes_block = get_first_knot_hes(z_knot, O.t(i), params(i), lambda_knot, O.t(i+1));
      con_start = con_end + 1;
    elseif i == O.N_knots
      con_end = con_start + O.End_knot_constraints - 1;
      lambda_knot = O.lambda(con_start:con_end);
      hes_block = get_end_knot_hes(z_knot, O.t(i-1), params(i), lambda_knot, O.t(i), lambda_last);
    else
      con_end = con_start + O.Middle_knot_constraints - 1;
      lambda_knot = O.lambda(con_start:con_end);
      hes_block = get_middle_knot_hes(z_knot, O.t(i), params(i), lambda_knot, O.t(i+1), lambda_last);
      con_start = con_end + 1;
    end

    hessian(z_start:z_end,z_start:z_end) = hes_block;
    lambda_last = lambda_knot;
  end
end

end
