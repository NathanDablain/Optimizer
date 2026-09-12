function [cost, grad, eq, jacobian, hessian] = Evaluate_3D_Rocket(O, option)
% Initialize outputs
cost = [];
grad = [];
eq = [];
jacobian = [];
hessian = [];
% knot vars are :
% 1 - p_n,
% 2 - p_e,
% 3 - p_u,
% 4 - s,
% 5 - psi
% 6 - gam
% 7 - Clp
% 8 - Clg
% 9 - s_pu_lb,
% 10 - s_Cl_ub,
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                       Evaluate cost                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  mu = 1.0e-10;
  function contribution = get_first_knot_cost(z_knot)
    contribution = -mu*(log(z_knot(9)) + log(z_knot(10)) + log(z_knot(11)) +...
                        log(z_knot(12)) + log(z_knot(13)));
  end

  function contribution = get_middle_knot_cost(z_knot)
    contribution = -mu*(log(z_knot(9)) + log(z_knot(10)) + log(z_knot(11)) +...
                        log(z_knot(12)) + log(z_knot(13)));
  end

  function contribution = get_end_knot_cost(z_knot, xd)
    contribution = (z_knot(1) - xd(1))^2 +...
                   (z_knot(2) - xd(2))^2 +...
                   (z_knot(3) - xd(3))^2 -...
                    mu*(log(z_knot(9)) + log(z_knot(10)) + log(z_knot(11)) +...
                        log(z_knot(12)) + log(z_knot(13)));
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
  grad_vec(9) = -mu / z_knot(9);
  grad_vec(10) = -mu / z_knot(10);
  grad_vec(11) = -mu / z_knot(11);
  grad_vec(12) = -mu / z_knot(12);
  grad_vec(13) = -mu / z_knot(13);
end

function grad_vec = get_middle_knot_grad(z_knot)
  grad_vec = zeros(length(z_knot), 1);
  grad_vec(9) = -mu / z_knot(9);
  grad_vec(10) = -mu / z_knot(10);
  grad_vec(11) = -mu / z_knot(11);
  grad_vec(12) = -mu / z_knot(12);
  grad_vec(13) = -mu / z_knot(13);
end

function grad_vec = get_end_knot_grad(z_knot, xd)
  grad_vec = zeros(length(z_knot), 1);
  grad_vec(1) = 2.0 * z_knot(1) - 2.0 * xd(1);
  grad_vec(2) = 2.0 * z_knot(2) - 2.0 * xd(2);
  grad_vec(3) = 2.0 * z_knot(3) - 2.0 * xd(3);
  grad_vec(9) = -mu / z_knot(9);
  grad_vec(10) = -mu / z_knot(10);
  grad_vec(11) = -mu / z_knot(11);
  grad_vec(12) = -mu / z_knot(12);
  grad_vec(13) = -mu / z_knot(13);
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
A = 0.25;
function knot_params = get_knot_params(z_knot, t_knot)
  gravity0 = 9.8065;
  Re = 6371e3;

  knot_params.g = gravity0*((Re/(Re + z_knot(3)))^2);
  knot_params.g_part_h = -2*Re*gravity0/((Re + z_knot(3))^3);
  knot_params.g_part_h2 = (6*Re*gravity0)/((Re + z_knot(3))^4);

  if z_knot(3) < 11000
    c1 = 15.04;
    c2 = 0.00649;
    c3 = 101.29;
    c4 = 273.1;
    c5 = 288.08;
    c6 = 5.256;
    c7 = 0.2869;
    Temperature = c1 - c2*z_knot(3);
    Pressure =c3*(((Temperature+c4)/c5)^c6);

    knot_params.rho_part_h = (c2*c3*((c1 - c2*z_knot(3) + c4)/c5)^c6*...
                             (1 - c6))/(c7*(c1-c2*z_knot(3)+c4)^2);

    knot_params.rho_part_h2 = (c2^2*c3*((c1 - c2*z_knot(3) + c4)/c5)^c6*...
                             (c6 - 1)*(c6 - 2))/(c7*(c1-c2*z_knot(3)+c4)^3);
  elseif z_knot(3) >= 11000 && z_knot(3) < 25000
    c1 = -56.46;
    c2 = 22.65;
    c3 = 1.73;
    c4 = 0.000157;
    c5 = 0.2869;
    c6 = 273.1;
    Temperature = c1;
    Pressure = c2*exp(c3 - c4*z_knot(3));

    knot_params.rho_part_h = (-c2*c4*exp(c3 - c4*z_knot(3)))/(c5*(c1+c6));
    knot_params.rho_part_h2 = (c2*c4^2*exp(c3 - c4*z_knot(3)))/(c5*(c1+c6));
  else
    c1 = -131.21;
    c2 = 0.00299;
    c3 = 2.488;
    c4 = 273.1;
    c5 = 216.6;
    c6 = -11.388;
    c7 = 0.2869;
    Temperature = c1 + c2*z_knot(3);
    Pressure = c3*(((Temperature+c4)/c5)^c6);

    knot_params.rho_part_h = (c2*c3*((c1 + c2*z_knot(3) + c4)/c5)^c6*...
                             (c6 - 1))/(c7*(c1+c2*z_knot(3)+c4)^2);
    knot_params.rho_part_h2 = (c2^2*c3*((c1 + c2*z_knot(3) + c4)/c5)^c6*...
                             (c6 - 1)*(c6 - 2))/(c7*(c1+c2*z_knot(3)+c4)^3);
  end
  knot_params.rho = Pressure/(0.2869*(Temperature+273.1));

  speed_of_sound = sqrt(1.4*287*(Temperature+273.1));
  Mach = z_knot(4) / speed_of_sound;
  if Mach < 0.8
    knot_params.c_D = 0.22;

    knot_params.c_D_part_s = 0;
    knot_params.c_D_part_s2 = 0;
  elseif Mach >= 0.8 && Mach < 1.2
    c1 = 0.22;
    c2 = 0.48;
    c3 = 0.8;
    knot_params.c_D = c1 + c2*sin(pi*(Mach - c3)/c3)^2;

    knot_params.c_D_part_s = (pi*c2*sin(2*pi*z_knot(4)/(speed_of_sound*c3)))/(speed_of_sound*c3);
    knot_params.c_D_part_s2 = (2*pi^2*c2*cos((2*pi*z_knot(4))/(speed_of_sound*c3)))/(speed_of_sound^2*c3^2);
  elseif Mach >= 1.2
    c1 = 0.25;
    c2 = 0.54;
    c3 = 1.2;
    knot_params.c_D = c1 + c2/(Mach^c3);

    knot_params.c_D_part_s = (-c2*c3*((z_knot(4)/speed_of_sound)^-c3))/z_knot(4);
    knot_params.c_D_part_s2 = (c2*c3*((z_knot(4)/speed_of_sound)^-c3)*(c3 + 1))/(z_knot(4)^2);
  end

  % Made up thrust mass tables to take us supersonic
  t_table = [0 0.2  0.5  2.5 3   3.25 4   6   8  10  11  12   13 13.5];
  T_table = [0 300 1000 1000 800 600 550 525 500 450 350 250 100 0];
  m_table = [15 14.92 14.52 11.8533 11.32 11.12 10.57 9.17 7.8367 6.6367 6.17 5.8367 5.7033 5.7033];
  if t_knot >= t_table(end)
    knot_params.T = T_table(end);
    knot_params.m = m_table(end);
  else
    knot_params.T = interp1(t_table, T_table, t_knot);
    knot_params.m = interp1(t_table, m_table, t_knot);
  end

end


function dx = get_dx(z_knot, knot_params)

  Drag = 0.5*knot_params.c_D*knot_params.rho*A*z_knot(4)*z_knot(4);

  dx(1) = z_knot(4) * cos(z_knot(6)) * cos(z_knot(5));
  dx(2) = z_knot(4) * cos(z_knot(6)) * sin(z_knot(5));
  dx(3) = z_knot(4) * sin(z_knot(6));
  dx(4) = ((knot_params.T - Drag)/knot_params.m) - (knot_params.g*sin(z_knot(6)));
  dx(5) = (knot_params.rho * A * z_knot(4) * z_knot(7)) / (knot_params.m * 2 * cos(z_knot(6)));
  dx(6) = ((0.5 * knot_params.rho * A * z_knot(4) * z_knot(8))/knot_params.m) - ((knot_params.g * cos(z_knot(6))) / z_knot(4));
end


function eq_vec = get_first_knot_eq(z_knot, t_knot, knot_params,...
                                    z_next, t_next, next_params, lb, ub, ic)

  h = t_next - t_knot;
  dx1 = get_dx(z_knot(1:8), knot_params);
  dx2 = get_dx(z_next(1:8), next_params);

  eq_vec = [...
            % the trapezoidal constraints for the knot
            z_next(1) - z_knot(1) - 0.5*h*(dx1(1) + dx2(1));...
            z_next(2) - z_knot(2) - 0.5*h*(dx1(2) + dx2(2));...
            z_next(3) - z_knot(3) - 0.5*h*(dx1(3) + dx2(3));...
            z_next(4) - z_knot(4) - 0.5*h*(dx1(4) + dx2(4));...
            z_next(5) - z_knot(5) - 0.5*h*(dx1(5) + dx2(5));...
            z_next(6) - z_knot(6) - 0.5*h*(dx1(6) + dx2(6));...
            % the slack variable constraints for the knot
            % first lower bound
            lb(1) - z_knot(3) + z_knot(9);...
            lb(2) - z_knot(7) + z_knot(10);...
            lb(3) - z_knot(8) + z_knot(11);...
            % then upper bound
            z_knot(7) - ub(1) + z_knot(12);...
            z_knot(8) - ub(2) + z_knot(13);...
            % the initial condition constraints for the knot
            ic(1) - z_knot(1);...
            ic(2) - z_knot(2);...
            ic(3) - z_knot(3);...
            ic(4) - z_knot(4);...
            ic(5) - z_knot(5);...
            ic(6) - z_knot(6)];
end

function eq_vec = get_middle_knot_eq(z_knot, t_knot, knot_params,...
                                     z_next, t_next, next_params, lb, ub)

  h = t_next - t_knot;
  dx1 = get_dx(z_knot(1:8), knot_params);
  dx2 = get_dx(z_next(1:8), next_params);

  eq_vec = [...
            % the trapezoidal constraints for the knot
            z_next(1) - z_knot(1) - 0.5*h*(dx1(1) + dx2(1));...
            z_next(2) - z_knot(2) - 0.5*h*(dx1(2) + dx2(2));...
            z_next(3) - z_knot(3) - 0.5*h*(dx1(3) + dx2(3));...
            z_next(4) - z_knot(4) - 0.5*h*(dx1(4) + dx2(4));...
            z_next(5) - z_knot(5) - 0.5*h*(dx1(5) + dx2(5));...
            z_next(6) - z_knot(6) - 0.5*h*(dx1(6) + dx2(6));...
            % the slack variable constraints for the knot
            % first lower bound
            lb(1) - z_knot(3) + z_knot(9);...
            lb(2) - z_knot(7) + z_knot(10);...
            lb(3) - z_knot(8) + z_knot(11);...
            % then upper bound
            z_knot(7) - ub(1) + z_knot(12);...
            z_knot(8) - ub(2) + z_knot(13)];
end

function eq_vec = get_end_knot_eq(z_knot, t_knot, lb, ub)

  eq_vec = [...
            % the slack variable constraints for the knot
            % first lower bound
            lb(1) - z_knot(3) + z_knot(9);...
            lb(2) - z_knot(7) + z_knot(10);...
            lb(3) - z_knot(8) + z_knot(11);...
            % then upper bound
            z_knot(7) - ub(1) + z_knot(12);...
            z_knot(8) - ub(2) + z_knot(13)];
end

if option.eval_eq
  params(1:O.N_knots) = struct('T', 0, 'm', 0, 'g', 0, 'rho', 0, 'c_D', 0,...
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
      params(i) = get_knot_params(z_knot, O.t(i));
      params(i+1) = get_knot_params(z_next, O.t(i+1));
      eq_vec = get_first_knot_eq(z_knot, O.t(i), params(i), z_next, O.t(i+1), params(i+1), O.lb, O.ub, O.ic);
    elseif i == O.N_knots
      eq_vec = get_end_knot_eq(z_knot, O.t(i), O.lb, O.ub);
    else
      z_start_next = z_end + 1;
      z_end_next = z_start_next + O.knot_size - 1;
      z_next = O.z(z_start_next:z_end_next);
      params(i+1) = get_knot_params(z_next,O.t(i+1));
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
  jac_block = zeros(17, 21);
  h = t_next - t_knot;
  pn1 = z_knot(1);
  pe1 = z_knot(2);
  pu1 = z_knot(3);
  s1 = z_knot(4);
  psi1 = z_knot(5);
  gam1 = z_knot(6);
  Clp1 = z_knot(7);
  Clg1 = z_knot(8);
  g1 = knot_params.g;
  m1 = knot_params.m;
  c_D1 = knot_params.c_D;
  rho1 = knot_params.rho;
  g_part_h1 = knot_params.g_part_h;
  rho_part_h1 = knot_params.rho_part_h;
  c_D_part_s1 = knot_params.c_D_part_s;

  pn2 = z_next(1);
  pe2 = z_next(2);
  pu2 = z_next(3);
  s2 = z_next(4);
  psi2 = z_next(5);
  gam2 = z_next(6);
  Clp2 = z_next(7);
  Clg2 = z_next(8);
  g2 = next_params.g;
  m2 = next_params.m;
  c_D2 = next_params.c_D;
  rho2 = next_params.rho;
  g_part_h2 = next_params.g_part_h;
  rho_part_h2 = next_params.rho_part_h;
  c_D_part_s2 = next_params.c_D_part_s;

  jac_block(1,1) = -1;
  jac_block(1,4) = -h*cos(gam1)*cos(psi1)/2;
  jac_block(1,5) = h*s1*sin(psi1)*cos(gam1)/2;
  jac_block(1,6) = h*s1*sin(gam1)*cos(psi1)/2;
  jac_block(1,14) = 1;
  jac_block(1,17) = -h*cos(gam2)*cos(psi2)/2;
  jac_block(1,18) = h*s2*sin(psi2)*cos(gam2)/2;
  jac_block(1,19) = h*s2*sin(gam2)*cos(psi2)/2;
  jac_block(2,2) = -1;
  jac_block(2,4) = -h*sin(psi1)*cos(gam1)/2;
  jac_block(2,5) = -h*s1*cos(gam1)*cos(psi1)/2;
  jac_block(2,6) = h*s1*sin(gam1)*sin(psi1)/2;
  jac_block(2,15) = 1;
  jac_block(2,17) = -h*sin(psi2)*cos(gam2)/2;
  jac_block(2,18) = -h*s2*cos(gam2)*cos(psi2)/2;
  jac_block(2,19) = h*s2*sin(gam2)*sin(psi2)/2;
  jac_block(3,3) = -1;
  jac_block(3,4) = -h*sin(gam1)/2;
  jac_block(3,6) = -h*s1*cos(gam1)/2;
  jac_block(3,16) = 1;
  jac_block(3,17) = -h*sin(gam2)/2;
  jac_block(3,19) = -h*s2*cos(gam2)/2;
  jac_block(4,3) = h*(A*s1^2*c_D1*rho_part_h1 + 2*m1*sin(gam1)*g_part_h1)/(4*m1);
  jac_block(4,4) = (A*h*s1*(s1*c_D_part_s1 + 2*c_D1)*rho1/4 - m1)/m1;
  jac_block(4,6) = h*g1*cos(gam1)/2;
  jac_block(4,16) = h*(A*s2^2*c_D2*rho_part_h2 + 2*m2*sin(gam2)*g_part_h2)/(4*m2);
  jac_block(4,17) = (A*h*s2*(s2*c_D_part_s2 + 2*c_D2)*rho2/4 + m2)/m2;
  jac_block(4,19) = h*g2*cos(gam2)/2;
  jac_block(5,3) = -A*Clp1*h*s1*rho_part_h1/(4*m1*cos(gam1));
  jac_block(5,4) = -A*Clp1*h*rho1/(4*m1*cos(gam1));
  jac_block(5,5) = -1;
  jac_block(5,6) = -A*Clp1*h*s1*rho1*sin(gam1)/(4*m1*cos(gam1)^2);
  jac_block(5,7) = -A*h*s1*rho1/(4*m1*cos(gam1));
  jac_block(5,16) = -A*Clp2*h*s2*rho_part_h2/(4*m2*cos(gam2));
  jac_block(5,17) = -A*Clp2*h*rho2/(4*m2*cos(gam2));
  jac_block(5,18) = 1;
  jac_block(5,19) = -A*Clp2*h*s2*rho2*sin(gam2)/(4*m2*cos(gam2)^2);
  jac_block(5,20) = -A*h*s2*rho2/(4*m2*cos(gam2));
  jac_block(6,3) = -A*Clg1*h*s1*rho_part_h1/(4*m1) + h*cos(gam1)*g_part_h1/(2*s1);
  jac_block(6,4) = -A*Clg1*h*rho1/(4*m1) - h*g1*cos(gam1)/(2*s1^2);
  jac_block(6,6) = (-h*g1*sin(gam1)/2 - s1)/s1;
  jac_block(6,8) = -A*h*s1*rho1/(4*m1);
  jac_block(6,16) = -A*Clg2*h*s2*rho_part_h2/(4*m2) + h*cos(gam2)*g_part_h2/(2*s2);
  jac_block(6,17) = -A*Clg2*h*rho2/(4*m2) - h*g2*cos(gam2)/(2*s2^2);
  jac_block(6,19) = (-h*g2*sin(gam2)/2 + s2)/s2;
  jac_block(6,21) = -A*h*s2*rho2/(4*m2);
  jac_block(7,3) = -1;
  jac_block(7,9) = 1;
  jac_block(8,7) = -1;
  jac_block(8,10) = 1;
  jac_block(9,8) = -1;
  jac_block(9,11) = 1;
  jac_block(10,7) = 1;
  jac_block(10,12) = 1;
  jac_block(11,8) = 1;
  jac_block(11,13) = 1;
  jac_block(12,1) = -1;
  jac_block(13,2) = -1;
  jac_block(14,3) = -1;
  jac_block(15,4) = -1;
  jac_block(16,5) = -1;
  jac_block(17,6) = -1;
end

function jac_block = get_middle_knot_jac(z_knot, t_knot, knot_params, z_next, t_next, next_params)
  jac_block = zeros(11, 21);
  h = t_next - t_knot;
  pn1 = z_knot(1);
  pe1 = z_knot(2);
  pu1 = z_knot(3);
  s1 = z_knot(4);
  psi1 = z_knot(5);
  gam1 = z_knot(6);
  Clp1 = z_knot(7);
  Clg1 = z_knot(8);
  g1 = knot_params.g;
  m1 = knot_params.m;
  c_D1 = knot_params.c_D;
  rho1 = knot_params.rho;
  g_part_h1 = knot_params.g_part_h;
  rho_part_h1 = knot_params.rho_part_h;
  c_D_part_s1 = knot_params.c_D_part_s;

  pn2 = z_next(1);
  pe2 = z_next(2);
  pu2 = z_next(3);
  s2 = z_next(4);
  psi2 = z_next(5);
  gam2 = z_next(6);
  Clp2 = z_next(7);
  Clg2 = z_next(8);
  g2 = next_params.g;
  m2 = next_params.m;
  c_D2 = next_params.c_D;
  rho2 = next_params.rho;
  g_part_h2 = next_params.g_part_h;
  rho_part_h2 = next_params.rho_part_h;
  c_D_part_s2 = next_params.c_D_part_s;

  jac_block(1,1) = -1;
  jac_block(1,4) = -h*cos(gam1)*cos(psi1)/2;
  jac_block(1,5) = h*s1*sin(psi1)*cos(gam1)/2;
  jac_block(1,6) = h*s1*sin(gam1)*cos(psi1)/2;
  jac_block(1,14) = 1;
  jac_block(1,17) = -h*cos(gam2)*cos(psi2)/2;
  jac_block(1,18) = h*s2*sin(psi2)*cos(gam2)/2;
  jac_block(1,19) = h*s2*sin(gam2)*cos(psi2)/2;
  jac_block(2,2) = -1;
  jac_block(2,4) = -h*sin(psi1)*cos(gam1)/2;
  jac_block(2,5) = -h*s1*cos(gam1)*cos(psi1)/2;
  jac_block(2,6) = h*s1*sin(gam1)*sin(psi1)/2;
  jac_block(2,15) = 1;
  jac_block(2,17) = -h*sin(psi2)*cos(gam2)/2;
  jac_block(2,18) = -h*s2*cos(gam2)*cos(psi2)/2;
  jac_block(2,19) = h*s2*sin(gam2)*sin(psi2)/2;
  jac_block(3,3) = -1;
  jac_block(3,4) = -h*sin(gam1)/2;
  jac_block(3,6) = -h*s1*cos(gam1)/2;
  jac_block(3,16) = 1;
  jac_block(3,17) = -h*sin(gam2)/2;
  jac_block(3,19) = -h*s2*cos(gam2)/2;
  jac_block(4,3) = h*(A*s1^2*c_D1*rho_part_h1 + 2*m1*sin(gam1)*g_part_h1)/(4*m1);
  jac_block(4,4) = (A*h*s1*(s1*c_D_part_s1 + 2*c_D1)*rho1/4 - m1)/m1;
  jac_block(4,6) = h*g1*cos(gam1)/2;
  jac_block(4,16) = h*(A*s2^2*c_D2*rho_part_h2 + 2*m2*sin(gam2)*g_part_h2)/(4*m2);
  jac_block(4,17) = (A*h*s2*(s2*c_D_part_s2 + 2*c_D2)*rho2/4 + m2)/m2;
  jac_block(4,19) = h*g2*cos(gam2)/2;
  jac_block(5,3) = -A*Clp1*h*s1*rho_part_h1/(4*m1*cos(gam1));
  jac_block(5,4) = -A*Clp1*h*rho1/(4*m1*cos(gam1));
  jac_block(5,5) = -1;
  jac_block(5,6) = -A*Clp1*h*s1*rho1*sin(gam1)/(4*m1*cos(gam1)^2);
  jac_block(5,7) = -A*h*s1*rho1/(4*m1*cos(gam1));
  jac_block(5,16) = -A*Clp2*h*s2*rho_part_h2/(4*m2*cos(gam2));
  jac_block(5,17) = -A*Clp2*h*rho2/(4*m2*cos(gam2));
  jac_block(5,18) = 1;
  jac_block(5,19) = -A*Clp2*h*s2*rho2*sin(gam2)/(4*m2*cos(gam2)^2);
  jac_block(5,20) = -A*h*s2*rho2/(4*m2*cos(gam2));
  jac_block(6,3) = -A*Clg1*h*s1*rho_part_h1/(4*m1) + h*cos(gam1)*g_part_h1/(2*s1);
  jac_block(6,4) = -A*Clg1*h*rho1/(4*m1) - h*g1*cos(gam1)/(2*s1^2);
  jac_block(6,6) = (-h*g1*sin(gam1)/2 - s1)/s1;
  jac_block(6,8) = -A*h*s1*rho1/(4*m1);
  jac_block(6,16) = -A*Clg2*h*s2*rho_part_h2/(4*m2) + h*cos(gam2)*g_part_h2/(2*s2);
  jac_block(6,17) = -A*Clg2*h*rho2/(4*m2) - h*g2*cos(gam2)/(2*s2^2);
  jac_block(6,19) = (-h*g2*sin(gam2)/2 + s2)/s2;
  jac_block(6,21) = -A*h*s2*rho2/(4*m2);
  jac_block(7,3) = -1;
  jac_block(7,9) = 1;
  jac_block(8,7) = -1;
  jac_block(8,10) = 1;
  jac_block(9,8) = -1;
  jac_block(9,11) = 1;
  jac_block(10,7) = 1;
  jac_block(10,12) = 1;
  jac_block(11,8) = 1;
  jac_block(11,13) = 1;

end

function jac_block = get_end_knot_jac(z_knot)
  jac_block = zeros(5, 13);

  jac_block(1,3) = -1;
  jac_block(1,9) = 1;
  jac_block(2,7) = -1;
  jac_block(2,10) = 1;
  jac_block(3,8) = -1;
  jac_block(3,11) = 1;
  jac_block(4,7) = 1;
  jac_block(4,12) = 1;
  jac_block(5,8) = 1;
  jac_block(5,13) = 1;
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
      jac_block = get_end_knot_jac(z_knot);
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
  hes_block = zeros(10, 10);
  h = t_next - t_knot;
  pn1 = z_knot(1);
  pe1 = z_knot(2);
  pu1 = z_knot(3);
  s1 = z_knot(4);
  psi1 = z_knot(5);
  gam1 = z_knot(6);
  Clp1 = z_knot(7);
  Clg1 = z_knot(8);
  s_pu_lb = z_knot(9);
  s_Clp_lb = z_knot(10);
  s_Clg_lb = z_knot(11);
  s_Clp_ub = z_knot(12);
  s_Clg_ub = z_knot(13);
  g1 = knot_params.g;
  m1 = knot_params.m;
  c_D1 = knot_params.c_D;
  rho1 = knot_params.rho;
  g_part_h = knot_params.g_part_h;
  rho_part_h = knot_params.rho_part_h;
  c_D_part_s = knot_params.c_D_part_s;
  g_part_h2 = knot_params.g_part_h2;
  rho_part_h2 = knot_params.rho_part_h2;
  c_D_part_s2 = knot_params.c_D_part_s2;
  lambda1 = lambda_knot(1);
  lambda2 = lambda_knot(2);
  lambda3 = lambda_knot(3);
  lambda4 = lambda_knot(4);
  lambda5 = lambda_knot(5);
  lambda6 = lambda_knot(6);
  lambda7 = lambda_knot(7);
  lambda8 = lambda_knot(8);

  hes_block(3,3) = h*(-A*Clg1*lambda6*s1^2*rho_part_h2 - A*Clp1*lambda5*s1^2*rho_part_h2/cos(gam1) +...
                   A*lambda4*s1^3*c_D1*rho_part_h2 + 2*lambda4*m1*s1*sin(gam1)*g_part_h2 +...
                   2*lambda6*m1*cos(gam1)*g_part_h2)/(4*m1*s1);
  hes_block(3,4) = h*(-A*Clg1*lambda6*s1^2*rho_part_h - A*Clp1*lambda5*s1^2*rho_part_h/cos(gam1) +...
                   A*lambda4*s1^4*c_D_part_s*rho_part_h + 2*A*lambda4*s1^3*c_D1*rho_part_h -...
                   2*lambda6*m1*cos(gam1)*g_part_h)/(4*m1*s1^2);
  hes_block(3,6) = -A*Clp1*h*lambda5*s1*sin(gam1)*rho_part_h/(4*m1*cos(gam1)^2) + h*lambda4*cos(gam1)*...
                    g_part_h/2 - h*lambda6*sin(gam1)*g_part_h/(2*s1);
  hes_block(3,7) = -A*h*lambda5*s1*rho_part_h/(4*m1*cos(gam1));
  hes_block(3,8) = -A*h*lambda6*s1*rho_part_h/(4*m1);
  hes_block(4,3) = hes_block(3,4);
  hes_block(4,4) = h*(A*lambda4*s1^3*(s1^2*c_D_part_s2 + 4*s1*c_D_part_s + 2*c_D1)*rho1 + 4*lambda6*m1*...
                   g1*cos(gam1))/(4*m1*s1^3);
  hes_block(4,5) = h*(lambda1*sin(psi1) - lambda2*cos(psi1))*cos(gam1)/2;
  hes_block(4,6) = -A*Clp1*h*lambda5*rho1*sin(gam1)/(4*m1*cos(gam1)^2) + h*lambda1*sin(gam1)*cos(psi1)/2 +...
                    h*lambda2*sin(gam1)*sin(psi1)/2 - h*lambda3*cos(gam1)/2 + h*lambda6*g1*sin(gam1)/(2*s1^2);
  hes_block(4,7) = -A*h*lambda5*rho1/(4*m1*cos(gam1));
  hes_block(4,8) = -A*h*lambda6*rho1/(4*m1);
  hes_block(5,4) = hes_block(4,5);
  hes_block(5,5) = h*s1*(lambda1*cos(psi1) + lambda2*sin(psi1))*cos(gam1)/2;
  hes_block(5,6) = h*s1*(-lambda1*sin(psi1) + lambda2*cos(psi1))*sin(gam1)/2;
  hes_block(6,3) = hes_block(3,6);
  hes_block(6,4) = hes_block(4,6);
  hes_block(6,5) = hes_block(5,6);
  hes_block(6,6) = h*(-2*A*Clp1*lambda5*s1^2*rho1*sin(gam1)^2 - A*Clp1*lambda5*s1^2*rho1*cos(gam1)^2 -...
                   2*lambda6*m1*g1*cos(gam1)^4 + 2*m1*s1*(lambda1*s1*cos(gam1)*cos(psi1) + lambda2*s1*...
                   sin(psi1)*cos(gam1) + lambda3*s1*sin(gam1) - lambda4*g1*sin(gam1))*cos(gam1)^3)/(4*m1*s1*cos(gam1)^3);
  hes_block(6,7) = -A*h*lambda5*s1*rho1*sin(gam1)/(4*m1*cos(gam1)^2);
  hes_block(7,3) = hes_block(3,7);
  hes_block(7,4) = hes_block(4,7);
  hes_block(7,6) = hes_block(6,7);
  hes_block(8,3) = hes_block(3,8);
  hes_block(8,4) = hes_block(4,8);
  hes_block(9,9) = mu/s_pu_lb^2;
  hes_block(10,10) = mu/s_Clp_lb^2;
  hes_block(11,11) = mu/s_Clg_lb^2;
  hes_block(12,12) = mu/s_Clp_ub^2;
  hes_block(13,13) = mu/s_Clg_ub^2;
end

function hes_block = get_middle_knot_hes(z_knot, t_knot, knot_params, lambda_knot, t_next, lambda_last)
  hes_block = zeros(7, 7);
  h = t_next - t_knot;
  pn1 = z_knot(1);
  pe1 = z_knot(2);
  pu1 = z_knot(3);
  s1 = z_knot(4);
  psi1 = z_knot(5);
  gam1 = z_knot(6);
  Clp1 = z_knot(7);
  Clg1 = z_knot(8);
  s_pu_lb = z_knot(9);
  s_Clp_lb = z_knot(10);
  s_Clg_lb = z_knot(11);
  s_Clp_ub = z_knot(12);
  s_Clg_ub = z_knot(13);
  g1 = knot_params.g;
  m1 = knot_params.m;
  c_D1 = knot_params.c_D;
  rho1 = knot_params.rho;
  g_part_h = knot_params.g_part_h;
  rho_part_h = knot_params.rho_part_h;
  c_D_part_s = knot_params.c_D_part_s;
  g_part_h2 = knot_params.g_part_h2;
  rho_part_h2 = knot_params.rho_part_h2;
  c_D_part_s2 = knot_params.c_D_part_s2;
  lambda1 = lambda_knot(1) + lambda_last(1);
  lambda2 = lambda_knot(2) + lambda_last(2);
  lambda3 = lambda_knot(3) + lambda_last(3);
  lambda4 = lambda_knot(4) + lambda_last(4);
  lambda5 = lambda_knot(5) + lambda_last(5);
  lambda6 = lambda_knot(6) + lambda_last(6);
  lambda7 = lambda_knot(7);
  lambda8 = lambda_knot(8);

  hes_block(3,3) = h*(-A*Clg1*lambda6*s1^2*rho_part_h2 - A*Clp1*lambda5*s1^2*rho_part_h2/cos(gam1) +...
                   A*lambda4*s1^3*c_D1*rho_part_h2 + 2*lambda4*m1*s1*sin(gam1)*g_part_h2 +...
                   2*lambda6*m1*cos(gam1)*g_part_h2)/(4*m1*s1);
  hes_block(3,4) = h*(-A*Clg1*lambda6*s1^2*rho_part_h - A*Clp1*lambda5*s1^2*rho_part_h/cos(gam1) +...
                   A*lambda4*s1^4*c_D_part_s*rho_part_h + 2*A*lambda4*s1^3*c_D1*rho_part_h -...
                   2*lambda6*m1*cos(gam1)*g_part_h)/(4*m1*s1^2);
  hes_block(3,6) = -A*Clp1*h*lambda5*s1*sin(gam1)*rho_part_h/(4*m1*cos(gam1)^2) + h*lambda4*cos(gam1)*...
                    g_part_h/2 - h*lambda6*sin(gam1)*g_part_h/(2*s1);
  hes_block(3,7) = -A*h*lambda5*s1*rho_part_h/(4*m1*cos(gam1));
  hes_block(3,8) = -A*h*lambda6*s1*rho_part_h/(4*m1);
  hes_block(4,3) = hes_block(3,4);
  hes_block(4,4) = h*(A*lambda4*s1^3*(s1^2*c_D_part_s2 + 4*s1*c_D_part_s + 2*c_D1)*rho1 + 4*lambda6*m1*...
                   g1*cos(gam1))/(4*m1*s1^3);
  hes_block(4,5) = h*(lambda1*sin(psi1) - lambda2*cos(psi1))*cos(gam1)/2;
  hes_block(4,6) = -A*Clp1*h*lambda5*rho1*sin(gam1)/(4*m1*cos(gam1)^2) + h*lambda1*sin(gam1)*cos(psi1)/2 +...
                    h*lambda2*sin(gam1)*sin(psi1)/2 - h*lambda3*cos(gam1)/2 + h*lambda6*g1*sin(gam1)/(2*s1^2);
  hes_block(4,7) = -A*h*lambda5*rho1/(4*m1*cos(gam1));
  hes_block(4,8) = -A*h*lambda6*rho1/(4*m1);
  hes_block(5,4) = hes_block(4,5);
  hes_block(5,5) = h*s1*(lambda1*cos(psi1) + lambda2*sin(psi1))*cos(gam1)/2;
  hes_block(5,6) = h*s1*(-lambda1*sin(psi1) + lambda2*cos(psi1))*sin(gam1)/2;
  hes_block(6,3) = hes_block(3,6);
  hes_block(6,4) = hes_block(4,6);
  hes_block(6,5) = hes_block(5,6);
  hes_block(6,6) = h*(-2*A*Clp1*lambda5*s1^2*rho1*sin(gam1)^2 - A*Clp1*lambda5*s1^2*rho1*cos(gam1)^2 -...
                   2*lambda6*m1*g1*cos(gam1)^4 + 2*m1*s1*(lambda1*s1*cos(gam1)*cos(psi1) + lambda2*s1*...
                   sin(psi1)*cos(gam1) + lambda3*s1*sin(gam1) - lambda4*g1*sin(gam1))*cos(gam1)^3)/(4*m1*s1*cos(gam1)^3);
  hes_block(6,7) = -A*h*lambda5*s1*rho1*sin(gam1)/(4*m1*cos(gam1)^2);
  hes_block(7,3) = hes_block(3,7);
  hes_block(7,4) = hes_block(4,7);
  hes_block(7,6) = hes_block(6,7);
  hes_block(8,3) = hes_block(3,8);
  hes_block(8,4) = hes_block(4,8);
  hes_block(9,9) = mu/s_pu_lb^2;
  hes_block(10,10) = mu/s_Clp_lb^2;
  hes_block(11,11) = mu/s_Clg_lb^2;
  hes_block(12,12) = mu/s_Clp_ub^2;
  hes_block(13,13) = mu/s_Clg_ub^2;
end

function hes_block = get_end_knot_hes(z_knot, t_knot, knot_params, lambda_knot, t_next, lambda_last)
  hes_block = zeros(7, 7);
  h = t_next - t_knot;
  pn1 = z_knot(1);
  pe1 = z_knot(2);
  pu1 = z_knot(3);
  s1 = z_knot(4);
  psi1 = z_knot(5);
  gam1 = z_knot(6);
  Clp1 = z_knot(7);
  Clg1 = z_knot(8);
  s_pu_lb = z_knot(9);
  s_Clp_lb = z_knot(10);
  s_Clg_lb = z_knot(11);
  s_Clp_ub = z_knot(12);
  s_Clg_ub = z_knot(13);
  g1 = knot_params.g;
  m1 = knot_params.m;
  c_D1 = knot_params.c_D;
  rho1 = knot_params.rho;
  g_part_h = knot_params.g_part_h;
  rho_part_h = knot_params.rho_part_h;
  c_D_part_s = knot_params.c_D_part_s;
  g_part_h2 = knot_params.g_part_h2;
  rho_part_h2 = knot_params.rho_part_h2;
  c_D_part_s2 = knot_params.c_D_part_s2;
  lambda1 = lambda_last(1);
  lambda2 = lambda_last(2);
  lambda3 = lambda_last(3);
  lambda4 = lambda_last(4);
  lambda5 = lambda_last(5);
  lambda6 = lambda_last(6);
  lambda7 = lambda_knot(1);
  lambda8 = lambda_knot(2);

  hes_block(1,1) = 2;
  hes_block(2,2) = 2;
  hes_block(3,3) = 2 + h*(-A*Clg1*lambda6*s1^2*rho_part_h2 - A*Clp1*lambda5*s1^2*rho_part_h2/cos(gam1) +...
                   A*lambda4*s1^3*c_D1*rho_part_h2 + 2*lambda4*m1*s1*sin(gam1)*g_part_h2 +...
                   2*lambda6*m1*cos(gam1)*g_part_h2)/(4*m1*s1);
  hes_block(3,4) = h*(-A*Clg1*lambda6*s1^2*rho_part_h - A*Clp1*lambda5*s1^2*rho_part_h/cos(gam1) +...
                   A*lambda4*s1^4*c_D_part_s*rho_part_h + 2*A*lambda4*s1^3*c_D1*rho_part_h -...
                   2*lambda6*m1*cos(gam1)*g_part_h)/(4*m1*s1^2);
  hes_block(3,6) = -A*Clp1*h*lambda5*s1*sin(gam1)*rho_part_h/(4*m1*cos(gam1)^2) + h*lambda4*cos(gam1)*...
                    g_part_h/2 - h*lambda6*sin(gam1)*g_part_h/(2*s1);
  hes_block(3,7) = -A*h*lambda5*s1*rho_part_h/(4*m1*cos(gam1));
  hes_block(3,8) = -A*h*lambda6*s1*rho_part_h/(4*m1);
  hes_block(4,3) = hes_block(3,4);
  hes_block(4,4) = h*(A*lambda4*s1^3*(s1^2*c_D_part_s2 + 4*s1*c_D_part_s + 2*c_D1)*rho1 + 4*lambda6*m1*...
                   g1*cos(gam1))/(4*m1*s1^3);
  hes_block(4,5) = h*(lambda1*sin(psi1) - lambda2*cos(psi1))*cos(gam1)/2;
  hes_block(4,6) = -A*Clp1*h*lambda5*rho1*sin(gam1)/(4*m1*cos(gam1)^2) + h*lambda1*sin(gam1)*cos(psi1)/2 +...
                    h*lambda2*sin(gam1)*sin(psi1)/2 - h*lambda3*cos(gam1)/2 + h*lambda6*g1*sin(gam1)/(2*s1^2);
  hes_block(4,7) = -A*h*lambda5*rho1/(4*m1*cos(gam1));
  hes_block(4,8) = -A*h*lambda6*rho1/(4*m1);
  hes_block(5,4) = hes_block(4,5);
  hes_block(5,5) = h*s1*(lambda1*cos(psi1) + lambda2*sin(psi1))*cos(gam1)/2;
  hes_block(5,6) = h*s1*(-lambda1*sin(psi1) + lambda2*cos(psi1))*sin(gam1)/2;
  hes_block(6,3) = hes_block(3,6);
  hes_block(6,4) = hes_block(4,6);
  hes_block(6,5) = hes_block(5,6);
  hes_block(6,6) = h*(-2*A*Clp1*lambda5*s1^2*rho1*sin(gam1)^2 - A*Clp1*lambda5*s1^2*rho1*cos(gam1)^2 -...
                   2*lambda6*m1*g1*cos(gam1)^4 + 2*m1*s1*(lambda1*s1*cos(gam1)*cos(psi1) + lambda2*s1*...
                   sin(psi1)*cos(gam1) + lambda3*s1*sin(gam1) - lambda4*g1*sin(gam1))*cos(gam1)^3)/(4*m1*s1*cos(gam1)^3);
  hes_block(6,7) = -A*h*lambda5*s1*rho1*sin(gam1)/(4*m1*cos(gam1)^2);
  hes_block(7,3) = hes_block(3,7);
  hes_block(7,4) = hes_block(4,7);
  hes_block(7,6) = hes_block(6,7);
  hes_block(8,3) = hes_block(3,8);
  hes_block(8,4) = hes_block(4,8);
  hes_block(9,9) = mu/s_pu_lb^2;
  hes_block(10,10) = mu/s_Clp_lb^2;
  hes_block(11,11) = mu/s_Clg_lb^2;
  hes_block(12,12) = mu/s_Clp_ub^2;
  hes_block(13,13) = mu/s_Clg_ub^2;
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
