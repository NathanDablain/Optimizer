function [cost, grad, eq, jacobian, hessian] = Evaluate_3D_Rocket_Simple(O, option)
% Initialize outputs
cost = [];
grad = [];
eq = [];
jacobian = [];
hessian = [];
% knot vars are :
% 1 - p_n,
% 2 - p_e,
% 3 - p_d,
% 4 - s,
% 5 - q0
% 6 - q1
% 7 - q2
% 8 - q3
% 9 - wy
% 10 - wz
% 11 - c_L
% 12 - sig
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                       Evaluate cost                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  mu = 5.0e-2;
  mu_q = 1.0;
  q_speed = 0.001;
  function contribution = get_first_knot_cost(z_knot, lb, ub)

    contribution = -mu*(log(z_knot(11) - lb(1)) +...
                        log(z_knot(12) - lb(2)) +...
                        log(ub(2) - z_knot(11)) +...
                        log(ub(3) - z_knot(12)));
  end

  function contribution = get_middle_knot_cost(z_knot, lb, ub)

    contribution = -mu*(log(z_knot(11) - lb(1)) +...
                        log(z_knot(12) - lb(2)) +...
                        log(ub(1) - z_knot(3)) +...
                        log(ub(2) - z_knot(11)) +...
                        log(ub(3) - z_knot(12))) +...
                    (mu_q*(norm(z_knot(5:8)) - 1.0))^2;
  end

  function contribution = get_end_knot_cost(z_knot, xd, lb, ub)
    contribution = (z_knot(1) - xd(1))^2 +...
                   (z_knot(2) - xd(2))^2 +...
                   (z_knot(3) - xd(3))^2 -...
                   q_speed*z_knot(4)^2 -...
                   mu*(log(z_knot(11) - lb(1)) +...
                       log(z_knot(12) - lb(2)) +...
                       log(ub(1) - z_knot(3)) +...
                       log(ub(2) - z_knot(11)) +...
                       log(ub(3) - z_knot(12))) +...
                    (mu_q*(norm(z_knot(5:8)) - 1.0))^2;

  end

if option.eval_cost
  cost = 0;
  for i = 1:O.N_knots
    z_start = O.knot_size*(i-1) + 1;
    z_end = O.knot_size*i;
    z_knot = O.z(z_start:z_end);
    if i == 1
      knot_cost_contribution = get_first_knot_cost(z_knot, O.lb, O.ub);
    elseif i == O.N_knots
      knot_cost_contribution = get_end_knot_cost(z_knot, O.xd, O.lb, O.ub);
    else
      knot_cost_contribution = get_middle_knot_cost(z_knot, O.lb, O.ub);
    end
    cost = cost + knot_cost_contribution;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Evaluate cost gradient wrt decision variables     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function grad_vec = get_first_knot_grad(z_knot, lb, ub)
  grad_vec = zeros(length(z_knot), 1);
  grad_vec(11) = (-mu/(z_knot(11) - ub(2))) + (mu/(lb(1) - z_knot(11)));
  grad_vec(12) = (-mu/(z_knot(12) - ub(3))) + (mu/(lb(2) - z_knot(12)));

end

function grad_vec = get_middle_knot_grad(z_knot, lb, ub)
  grad_vec = zeros(length(z_knot), 1);
  grad_vec(3) = -mu/(z_knot(3) - ub(1));
  grad_vec(5) = mu_q*4.0*z_knot(5)*(norm(z_knot(5:8)) - 1);
  grad_vec(6) = mu_q*4.0*z_knot(6)*(norm(z_knot(5:8)) - 1);
  grad_vec(7) = mu_q*4.0*z_knot(7)*(norm(z_knot(5:8)) - 1);
  grad_vec(8) = mu_q*4.0*z_knot(8)*(norm(z_knot(5:8)) - 1);
  grad_vec(11) = (-mu/(z_knot(11) - ub(2))) + (mu/(lb(1) - z_knot(11)));
  grad_vec(12) = (-mu/(z_knot(12) - ub(3))) + (mu/(lb(2) - z_knot(12)));

end

function grad_vec = get_end_knot_grad(z_knot, xd, lb, ub)
  grad_vec = zeros(length(z_knot), 1);
  grad_vec(1) = 2.0 * z_knot(1) - 2.0 * xd(1);
  grad_vec(2) = 2.0 * z_knot(2) - 2.0 * xd(2);
  grad_vec(3) = 2.0 * z_knot(3) - 2.0 * xd(3) - mu/(z_knot(3) - ub(1));
  grad_vec(4) = -2.0*q_speed*z_knot(4);
  grad_vec(5) = mu_q*4.0*z_knot(5)*(norm(z_knot(5:8)) - 1);
  grad_vec(6) = mu_q*4.0*z_knot(6)*(norm(z_knot(5:8)) - 1);
  grad_vec(7) = mu_q*4.0*z_knot(7)*(norm(z_knot(5:8)) - 1);
  grad_vec(8) = mu_q*4.0*z_knot(8)*(norm(z_knot(5:8)) - 1);
  grad_vec(11) = (-mu/(z_knot(11) - ub(2))) + (mu/(lb(1) - z_knot(11)));
  grad_vec(12) = (-mu/(z_knot(12) - ub(3))) + (mu/(lb(2) - z_knot(12)));

end

if option.eval_grad
  grad = zeros(O.N_decision_variables, 1);
  for i = 1:O.N_knots
    z_start = O.knot_size*(i-1) + 1;
    z_end = O.knot_size*i;
    z_knot = O.z(z_start:z_end);
    if i == 1
      grad_vec = get_first_knot_grad(z_knot, O.lb, O.ub);
    elseif i == O.N_knots
      grad_vec = get_end_knot_grad(z_knot, O.xd, O.lb, O.ub);
    else
      grad_vec = get_middle_knot_grad(z_knot, O.lb, O.ub);
    end
    grad(z_start:z_end) = grad_vec;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%               Evaluate equality constraints            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = 0.25;
pi_AR_e = pi * 2.0 * 0.5;
function knot_params = get_knot_params(z_knot, t_knot)
  gravity0 = 9.8065;
  Re = 6371e3;
  h = -z_knot(3);
  s = z_knot(4);

  knot_params.g = gravity0*((Re/(Re + h))^2);

  if h < 11000
    c1 = 15.04;
    c2 = 0.00649;
    c3 = 101.29;
    c4 = 273.1;
    c5 = 288.08;
    c6 = 5.256;
    c7 = 0.2869;
    Temperature = c1 - c2*h;
    Pressure =c3*(((Temperature+c4)/c5)^c6);

    knot_params.rho_part_h = (c2*c3*((c1 - c2*h + c4)/c5)^c6*...
                             (1 - c6))/(c7*(c1-c2*h+c4)^2);

    knot_params.rho_part_h2 = (c2^2*c3*((c1 - c2*h + c4)/c5)^c6*...
                             (c6 - 1)*(c6 - 2))/(c7*(c1-c2*h+c4)^3);
  elseif h >= 11000 && h < 25000
    c1 = -56.46;
    c2 = 22.65;
    c3 = 1.73;
    c4 = 0.000157;
    c5 = 0.2869;
    c6 = 273.1;
    Temperature = c1;
    Pressure = c2*exp(c3 - c4*h);

    knot_params.rho_part_h = (-c2*c4*exp(c3 - c4*h))/(c5*(c1+c6));
    knot_params.rho_part_h2 = (c2*c4^2*exp(c3 - c4*h))/(c5*(c1+c6));
  else
    c1 = -131.21;
    c2 = 0.00299;
    c3 = 2.488;
    c4 = 273.1;
    c5 = 216.6;
    c6 = -11.388;
    c7 = 0.2869;
    Temperature = c1 + c2*h;
    Pressure = c3*(((Temperature+c4)/c5)^c6);

    knot_params.rho_part_h = (c2*c3*((c1 + c2*h + c4)/c5)^c6*...
                             (c6 - 1))/(c7*(c1+c2*h+c4)^2);
    knot_params.rho_part_h2 = (c2^2*c3*((c1 + c2*h + c4)/c5)^c6*...
                             (c6 - 1)*(c6 - 2))/(c7*(c1+c2*h+c4)^3);
  end
  knot_params.rho = Pressure/(0.2869*(Temperature+273.1));

  speed_of_sound = sqrt(1.4*287*(Temperature+273.1));
  Mach = s / speed_of_sound;
  if Mach < 0.8
    knot_params.c_D = 0.22;

    knot_params.c_D_part_s = 0;
    knot_params.c_D_part_s2 = 0;
  elseif Mach >= 0.8 && Mach < 1.2
    c1 = 0.22;
    c2 = 0.48;
    c3 = 0.8;
    knot_params.c_D = c1 + c2*sin(pi*(Mach - c3)/c3)^2;

    knot_params.c_D_part_s = (pi*c2*sin(2*pi*s/(speed_of_sound*c3)))/(speed_of_sound*c3);
    knot_params.c_D_part_s2 = (2*pi^2*c2*cos((2*pi*s)/(speed_of_sound*c3)))/(speed_of_sound^2*c3^2);
  elseif Mach >= 1.2
    c1 = 0.25;
    c2 = 0.54;
    c3 = 1.2;
    knot_params.c_D = c1 + c2/(Mach^c3);

    knot_params.c_D_part_s = (-c2*c3*((s/speed_of_sound)^-c3))/s;
    knot_params.c_D_part_s2 = (c2*c3*((s/speed_of_sound)^-c3)*(c3 + 1))/(s^2);
  end

  knot_params.rho_part_h = - knot_params.rho_part_h;
  knot_params.rho_part_h2 = - knot_params.rho_part_h2;
  % Made up thrust mass tables to take us supersonic
  t_table = [0 0.2  0.5  2.5 3   3.25 4   6   8  10  11  12   13 13.5];
  T_table = 30.*[0 300 1000 1000 800 600 550 525 500 450 350 250 100 0];
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

  s = z_knot(4);
  q0 = z_knot(5);
  q1 = z_knot(6);
  q2 = z_knot(7);
  q3 = z_knot(8);
  wy = z_knot(9);
  wz = z_knot(10);
  c_L = z_knot(11);
  sig = z_knot(12);
  q = [q0 q1 q2 q3];
  D = 0.5*knot_params.rho*A*(knot_params.c_D + (c_L^2/pi_AR_e))*s*s;
  az = 0.5*knot_params.rho*A*cos(sig)*c_L*s*s/knot_params.m;
  ay = 0.5*knot_params.rho*A*sin(sig)*c_L*s*s/knot_params.m;

  v_NED = [s*(-2*q2^2 - 2*q3^2 + 1); 2*s*(q0*q3 + q1*q2); 2*s*(-q0*q2 + q1*q3)];
  omega = [0 0 -wy -wz; 0 0 wz -wy; wy -wz 0 0; wz wy 0 0];

  dx(1:3) = v_NED;
  dx(4) = ((knot_params.T - D)/knot_params.m) + 2*knot_params.g*(-q0*q2 + q1*q3);
  dx(5:8) = 0.5*omega*q';
  dx(9:10) = [(-az - knot_params.g*(-2*q1^2 - 2*q2^2 + 1))/s;...
              (ay + 2*knot_params.g*(q0*q1 + q2*q3))/s];
end

function eq_vec = get_first_knot_eq(z_knot, t_knot, knot_params,...
                                    z_next, t_next, next_params, lb, ub, ic)

  h = t_next - t_knot;
  dx1 = get_dx(z_knot(1:12), knot_params);
  dx2 = get_dx(z_next(1:12), next_params);
  wy = z_knot(9);
  wz = z_knot(10);
  w = [wy wz];

  theta = norm(w)*h;
  if norm(w) < 1.0e-5
    delta_q = zeros(4,1);
  else
    delta_q = [cos(theta/2);0;(w(1)/norm(w))*sin(theta/2);(w(2)/norm(w))*sin(theta/2)];
  end
  constraint = z_next(5:8) - quaternion_multiply(z_knot(5:8), delta_q, 'left');

  eq_vec = [...
            % the trapezoidal constraints for the knot
            z_next(1) - z_knot(1) - 0.5*h*(dx1(1) + dx2(1));...
            z_next(2) - z_knot(2) - 0.5*h*(dx1(2) + dx2(2));...
            z_next(3) - z_knot(3) - 0.5*h*(dx1(3) + dx2(3));...
            z_next(4) - z_knot(4) - 0.5*h*(dx1(4) + dx2(4));...
            constraint(1);...
            constraint(2);...
            constraint(3);...
            constraint(4);...
            % The angular rate constraints
            z_knot(9) - 0.5*(dx1(9) + dx2(9));...
            z_knot(10) - 0.5*(dx1(10) + dx2(10));...
            % the initial condition constraints for the knot
            ic(1) - z_knot(1);...
            ic(2) - z_knot(2);...
            ic(3) - z_knot(3);...
            ic(4) - z_knot(4);...
            ic(5) - z_knot(5);...
            ic(6) - z_knot(6);...
            ic(7) - z_knot(7);...
            ic(8) - z_knot(8)];
end

function eq_vec = get_middle_knot_eq(z_knot, t_knot, knot_params,...
                                     z_next, t_next, next_params, lb, ub)

  h = t_next - t_knot;
  dx1 = get_dx(z_knot(1:12), knot_params);
  dx2 = get_dx(z_next(1:12), next_params);

  wy = z_knot(9);
  wz = z_knot(10);
  w = [wy wz];

  theta = norm(w)*h;
  if norm(w) < 1.0e-5
    delta_q = zeros(4,1);
  else
    delta_q = [cos(theta/2);0;(w(1)/norm(w))*sin(theta/2);(w(2)/norm(w))*sin(theta/2)];
  end

  constraint = z_next(5:8) - quaternion_multiply(z_knot(5:8), delta_q, 'left');

  eq_vec = [...
            % the trapezoidal constraints for the knot
            z_next(1) - z_knot(1) - 0.5*h*(dx1(1) + dx2(1));...
            z_next(2) - z_knot(2) - 0.5*h*(dx1(2) + dx2(2));...
            z_next(3) - z_knot(3) - 0.5*h*(dx1(3) + dx2(3));...
            z_next(4) - z_knot(4) - 0.5*h*(dx1(4) + dx2(4));...
            constraint(1);...
            constraint(2);...
            constraint(3);...
            constraint(4);...
            % The angular rate constraints
            z_knot(9) - 0.5*(dx1(9) + dx2(9));...
            z_knot(10) - 0.5*(dx1(10) + dx2(10))];

end

function eq_vec = get_end_knot_eq(z_knot, t_knot, knot_params, lb, ub)

  eq_vec = [];
end

if option.eval_eq
  params(1:O.N_knots) = struct('T', 0, 'm', 0, 'g', 0, 'rho', 0, 'c_D', 0,...
                               'rho_part_h', 0, 'c_D_part_s', 0,...
                               'rho_part_h2', 0, 'c_D_part_s2', 0);
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
      eq_vec = get_end_knot_eq(z_knot, O.t(i), params(i), O.lb, O.ub);
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
  jac_block = zeros(18, 24);
  h = t_next - t_knot;
  pn1 = z_knot(1);
  pe1 = z_knot(2);
  pd1 = z_knot(3);
  s1 = z_knot(4);
  q01 = z_knot(5);
  q11 = z_knot(6);
  q21 = z_knot(7);
  q31 = z_knot(8);
  wy1 = z_knot(9);
  wz1 = z_knot(10);
  c_L1 = z_knot(11);
  sig1 = z_knot(12);
  g = knot_params.g;
  m1 = knot_params.m;
  c_D1 = knot_params.c_D;
  rho1 = knot_params.rho;
  rho_part_h1 = knot_params.rho_part_h;
  c_D_part_s1 = knot_params.c_D_part_s;

  pn2 = z_next(1);
  pe2 = z_next(2);
  pd2 = z_next(3);
  s2 = z_next(4);
  q02 = z_next(5);
  q12 = z_next(6);
  q22 = z_next(7);
  q32 = z_next(8);
  wy2 = z_next(9);
  wz2 = z_next(10);
  c_L2 = z_next(11);
  sig2 = z_next(12);
##  g2 = next_params.g;
  m2 = next_params.m;
  c_D2 = next_params.c_D;
  rho2 = next_params.rho;
  rho_part_h2 = next_params.rho_part_h;
  c_D_part_s2 = next_params.c_D_part_s;

  jac_block(1,1) = -1;
  jac_block(1,4) = h*(q21^2 + q31^2 - 1/2);
  jac_block(1,7) = 2*h*q21*s1;
  jac_block(1,8) = 2*h*q31*s1;
  jac_block(1,13) = 1;
  jac_block(1,16) = h*(q22^2 + q32^2 - 1/2);
  jac_block(1,19) = 2*h*q22*s2;
  jac_block(1,20) = 2*h*q32*s2;
  jac_block(2,2) = -1;
  jac_block(2,4) = -h*(q01*q31 + q11*q21);
  jac_block(2,5) = -h*q31*s1;
  jac_block(2,6) = -h*q21*s1;
  jac_block(2,7) = -h*q11*s1;
  jac_block(2,8) = -h*q01*s1;
  jac_block(2,14) = 1;
  jac_block(2,16) = -h*(q02*q32 + q12*q22);
  jac_block(2,17) = -h*q32*s2;
  jac_block(2,18) = -h*q22*s2;
  jac_block(2,19) = -h*q12*s2;
  jac_block(2,20) = -h*q02*s2;
  jac_block(3,3) = -1;
  jac_block(3,4) = h*(q01*q21 - q11*q31);
  jac_block(3,5) = h*q21*s1;
  jac_block(3,6) = -h*q31*s1;
  jac_block(3,7) = h*q01*s1;
  jac_block(3,8) = -h*q11*s1;
  jac_block(3,15) = 1;
  jac_block(3,16) = h*(q02*q22 - q12*q32);
  jac_block(3,17) = h*q22*s2;
  jac_block(3,18) = -h*q32*s2;
  jac_block(3,19) = h*q02*s2;
  jac_block(3,20) = -h*q12*s2;

  jac_block(4,3) = A*h*s1^2*(c_L1^2 + pi_AR_e*c_D1)*rho_part_h1/(4*m1*pi_AR_e);
  jac_block(4,4) = (A*h*s1*(2*c_L1^2 + pi_AR_e*s1*c_D_part_s1 + 2*pi_AR_e*c_D1)*rho1/4 - m1*pi_AR_e)/(m1*pi_AR_e);
  jac_block(4,5) = g*h*q21;
  jac_block(4,6) = -g*h*q31;
  jac_block(4,7) = g*h*q01;
  jac_block(4,8) = -g*h*q11;
  jac_block(4,11) = A*c_L1*h*s1^2*rho1/(2*m1*pi_AR_e);
  jac_block(4,15) = A*h*s2^2*(c_L2^2 + pi_AR_e*c_D2)*rho_part_h2/(4*m2*pi_AR_e);
  jac_block(4,16) = (A*h*s2*(2*c_L2^2 + pi_AR_e*s2*c_D_part_s2 + 2*pi_AR_e*c_D2)*rho2/4 + m2*pi_AR_e)/(m2*pi_AR_e);
  jac_block(4,17) = g*h*q22;
  jac_block(4,18) = -g*h*q32;
  jac_block(4,19) = g*h*q02;
  jac_block(4,20) = -g*h*q12;
  jac_block(4,23) = A*c_L2*h*s2^2*rho2/(2*m2*pi_AR_e);

  c1 = sqrt(wy1^2 + wz1^2);
  c2 = cos(h*c1/2);
  c3 = sin(h*c1/2);

  if c1 > 1.0e-5
    jac_block(5,5) = -c2;
    jac_block(5,7) = wy1*c3/c1;
    jac_block(5,8) = wz1*c3/c1;
    jac_block(5,9) = (h*q01*wy1^3*c3 + h*q01*wy1*wz1^2*c3 + h*q21*wy1^2*c1*c2 + h*q31*wy1*wz1*c1*c2 + 2*q21*wz1^2*c3 - 2*q31*wy1*wz1*c3)/(2*c1^3);
    jac_block(5,10) = (h*q01*wy1^2*wz1*c3 + h*q01*wz1^3*c3 + h*q21*wy1*wz1*c1*c2 + h*q31*wz1^2*c1*c2 - 2*q21*wy1*wz1*c3 + 2*q31*wy1^2*c3)/(2*c1^3);
    jac_block(5,17) = 1;
    jac_block(6,6) = -c2;
    jac_block(6,7) = -wz1*c3/c1;
    jac_block(6,8) = wy1*c3/c1;
    jac_block(6,9) = (h*q11*wy1^3*c3 + h*q11*wy1*wz1^2*c3 - h*q21*wy1*wz1*c1*c2 + h*q31*wy1^2*c1*c2 + 2*q21*wy1*wz1*c3 + 2*q31*wz1^2*c3)/(2*c1^3);
    jac_block(6,10) = (h*q11*wy1^2*wz1*c3 + h*q11*wz1^3*c3 - h*q21*wz1^2*c1*c2 + h*q31*wy1*wz1*c1*c2 - 2*q21*wy1^2*c3 - 2*q31*wy1*wz1*c3)/(2*c1^3);
    jac_block(6,18) = 1;
    jac_block(7,5) = -wy1*c3/c1;
    jac_block(7,6) = wz1*c3/c1;
    jac_block(7,7) = -c2;
    jac_block(7,9) = (-h*q01*wy1^2*c1*c2 + h*q11*wy1*wz1*c1*c2 + h*q21*wy1^3*c3 + h*q21*wy1*wz1^2*c3 - 2*q01*wz1^2*c3 - 2*q11*wy1*wz1*c3)/(2*c1^3);
    jac_block(7,10) = (-h*q01*wy1*wz1*c1*c2 + h*q11*wz1^2*c1*c2 + h*q21*wy1^2*wz1*c3 + h*q21*wz1^3*c3 + 2*q01*wy1*wz1*c3 + 2*q11*wy1^2*c3)/(2*c1^3);
    jac_block(7,19) = 1;
    jac_block(8,5) = -wz1*c3/c1;
    jac_block(8,6) = -wy1*c3/c1;
    jac_block(8,8) = -c2;
    jac_block(8,9) = (-h*q01*wy1*wz1*c1*c2 - h*q11*wy1^2*c1*c2 + h*q31*wy1^3*c3 + h*q31*wy1*wz1^2*c3 + 2*q01*wy1*wz1*c3 - 2*q11*wz1^2*c3)/(2*c1^3);
    jac_block(8,10) = (-h*q01*wz1^2*c1*c2 - h*q11*wy1*wz1*c1*c2 + h*q31*wy1^2*wz1*c3 + h*q31*wz1^3*c3 - 2*q01*wy1^2*c3 + 2*q11*wy1*wz1*c3)/(2*c1^3);
    jac_block(8,20) = 1;
  else
    jac_block(5,5) = -c2;
    jac_block(5,17) = 1;
    jac_block(6,6) = -c2;
    jac_block(6,18) = 1;
    jac_block(7,7) = -c2;
    jac_block(7,19) = 1;
    jac_block(8,8) = -c2;
    jac_block(8,20) = 1;
  end

  jac_block(9,3) = A*c_L1*s1*cos(sig1)*rho_part_h1/(4*m1);
  jac_block(9,4) = (A*c_L1*s1^2*rho1*cos(sig1) + 2*g*m1*(2*q11^2 + 2*q21^2 - 1))/(4*m1*s1^2);
  jac_block(9,6) = -2*g*q11/s1;
  jac_block(9,7) = -2*g*q21/s1;
  jac_block(9,9) = 1;
  jac_block(9,11) = A*s1*rho1*cos(sig1)/(4*m1);
  jac_block(9,12) = -A*c_L1*s1*rho1*sin(sig1)/(4*m1);
  jac_block(9,15) = A*c_L2*s2*cos(sig2)*rho_part_h2/(4*m2);
  jac_block(9,16) = (A*c_L2*s2^2*rho2*cos(sig2) + 2*g*m2*(2*q12^2 + 2*q22^2 - 1))/(4*m2*s2^2);
  jac_block(9,18) = -2*g*q12/s2;
  jac_block(9,19) = -2*g*q22/s2;
  jac_block(9,23) = A*s2*rho2*cos(sig2)/(4*m2);
  jac_block(9,24) = -A*c_L2*s2*rho2*sin(sig2)/(4*m2);

  jac_block(10,3) = -A*c_L1*s1*sin(sig1)*rho_part_h1/(4*m1);
  jac_block(10,4) = (-A*c_L1*s1^2*rho1*sin(sig1)/4 + g*m1*(q01*q11 + q21*q31))/(m1*s1^2);
  jac_block(10,5) = -g*q11/s1;
  jac_block(10,6) = -g*q01/s1;
  jac_block(10,7) = -g*q31/s1;
  jac_block(10,8) = -g*q21/s1;
  jac_block(10,10) = 1;
  jac_block(10,11) = -A*s1*rho1*sin(sig1)/(4*m1);
  jac_block(10,12) = -A*c_L1*s1*rho1*cos(sig1)/(4*m1);
  jac_block(10,15) = -A*c_L2*s2*sin(sig2)*rho_part_h2/(4*m2);
  jac_block(10,16) = (-A*c_L2*s2^2*rho2*sin(sig2)/4 + g*m2*(q02*q12 + q22*q32))/(m2*s2^2);
  jac_block(10,17) = -g*q12/s2;
  jac_block(10,18) = -g*q02/s2;
  jac_block(10,19) = -g*q32/s2;
  jac_block(10,20) = -g*q22/s2;
  jac_block(10,23) = -A*s2*rho2*sin(sig2)/(4*m2);
  jac_block(10,24) = -A*c_L2*s2*rho2*cos(sig2)/(4*m2);

  jac_block(11,1) = -1;
  jac_block(12,2) = -1;
  jac_block(13,3) = -1;
  jac_block(14,4) = -1;
  jac_block(15,5) = -1;
  jac_block(16,6) = -1;
  jac_block(17,7) = -1;
  jac_block(18,8) = -1;

end

function jac_block = get_middle_knot_jac(z_knot, t_knot, knot_params, z_next, t_next, next_params)
  jac_block = zeros(10, 24);
  h = t_next - t_knot;
  pn1 = z_knot(1);
  pe1 = z_knot(2);
  pd1 = z_knot(3);
  s1 = z_knot(4);
  q01 = z_knot(5);
  q11 = z_knot(6);
  q21 = z_knot(7);
  q31 = z_knot(8);
  wy1 = z_knot(9);
  wz1 = z_knot(10);
  c_L1 = z_knot(11);
  sig1 = z_knot(12);
  g = knot_params.g;
  m1 = knot_params.m;
  c_D1 = knot_params.c_D;
  rho1 = knot_params.rho;
  rho_part_h1 = knot_params.rho_part_h;
  c_D_part_s1 = knot_params.c_D_part_s;

  pn2 = z_next(1);
  pe2 = z_next(2);
  pd2 = z_next(3);
  s2 = z_next(4);
  q02 = z_next(5);
  q12 = z_next(6);
  q22 = z_next(7);
  q32 = z_next(8);
  wy2 = z_next(9);
  wz2 = z_next(10);
  c_L2 = z_next(11);
  sig2 = z_next(12);
##  g2 = next_params.g;
  m2 = next_params.m;
  c_D2 = next_params.c_D;
  rho2 = next_params.rho;
  rho_part_h2 = next_params.rho_part_h;
  c_D_part_s2 = next_params.c_D_part_s;

  jac_block(1,1) = -1;
  jac_block(1,4) = h*(q21^2 + q31^2 - 1/2);
  jac_block(1,7) = 2*h*q21*s1;
  jac_block(1,8) = 2*h*q31*s1;
  jac_block(1,13) = 1;
  jac_block(1,16) = h*(q22^2 + q32^2 - 1/2);
  jac_block(1,19) = 2*h*q22*s2;
  jac_block(1,20) = 2*h*q32*s2;
  jac_block(2,2) = -1;
  jac_block(2,4) = -h*(q01*q31 + q11*q21);
  jac_block(2,5) = -h*q31*s1;
  jac_block(2,6) = -h*q21*s1;
  jac_block(2,7) = -h*q11*s1;
  jac_block(2,8) = -h*q01*s1;
  jac_block(2,14) = 1;
  jac_block(2,16) = -h*(q02*q32 + q12*q22);
  jac_block(2,17) = -h*q32*s2;
  jac_block(2,18) = -h*q22*s2;
  jac_block(2,19) = -h*q12*s2;
  jac_block(2,20) = -h*q02*s2;
  jac_block(3,3) = -1;
  jac_block(3,4) = h*(q01*q21 - q11*q31);
  jac_block(3,5) = h*q21*s1;
  jac_block(3,6) = -h*q31*s1;
  jac_block(3,7) = h*q01*s1;
  jac_block(3,8) = -h*q11*s1;
  jac_block(3,15) = 1;
  jac_block(3,16) = h*(q02*q22 - q12*q32);
  jac_block(3,17) = h*q22*s2;
  jac_block(3,18) = -h*q32*s2;
  jac_block(3,19) = h*q02*s2;
  jac_block(3,20) = -h*q12*s2;

  jac_block(4,3) = A*h*s1^2*(c_L1^2 + pi_AR_e*c_D1)*rho_part_h1/(4*m1*pi_AR_e);
  jac_block(4,4) = (A*h*s1*(2*c_L1^2 + pi_AR_e*s1*c_D_part_s1 + 2*pi_AR_e*c_D1)*rho1/4 - m1*pi_AR_e)/(m1*pi_AR_e);
  jac_block(4,5) = g*h*q21;
  jac_block(4,6) = -g*h*q31;
  jac_block(4,7) = g*h*q01;
  jac_block(4,8) = -g*h*q11;
  jac_block(4,11) = A*c_L1*h*s1^2*rho1/(2*m1*pi_AR_e);
  jac_block(4,15) = A*h*s2^2*(c_L2^2 + pi_AR_e*c_D2)*rho_part_h2/(4*m2*pi_AR_e);
  jac_block(4,16) = (A*h*s2*(2*c_L2^2 + pi_AR_e*s2*c_D_part_s2 + 2*pi_AR_e*c_D2)*rho2/4 + m2*pi_AR_e)/(m2*pi_AR_e);
  jac_block(4,17) = g*h*q22;
  jac_block(4,18) = -g*h*q32;
  jac_block(4,19) = g*h*q02;
  jac_block(4,20) = -g*h*q12;
  jac_block(4,23) = A*c_L2*h*s2^2*rho2/(2*m2*pi_AR_e);

  c1 = sqrt(wy1^2 + wz1^2);
  c2 = cos(h*c1/2);
  c3 = sin(h*c1/2);

  if c1 > 1.0e-5
    jac_block(5,5) = -c2;
    jac_block(5,7) = wy1*c3/c1;
    jac_block(5,8) = wz1*c3/c1;
    jac_block(5,9) = (h*q01*wy1^3*c3 + h*q01*wy1*wz1^2*c3 + h*q21*wy1^2*c1*c2 + h*q31*wy1*wz1*c1*c2 + 2*q21*wz1^2*c3 - 2*q31*wy1*wz1*c3)/(2*c1^3);
    jac_block(5,10) = (h*q01*wy1^2*wz1*c3 + h*q01*wz1^3*c3 + h*q21*wy1*wz1*c1*c2 + h*q31*wz1^2*c1*c2 - 2*q21*wy1*wz1*c3 + 2*q31*wy1^2*c3)/(2*c1^3);
    jac_block(5,17) = 1;
    jac_block(6,6) = -c2;
    jac_block(6,7) = -wz1*c3/c1;
    jac_block(6,8) = wy1*c3/c1;
    jac_block(6,9) = (h*q11*wy1^3*c3 + h*q11*wy1*wz1^2*c3 - h*q21*wy1*wz1*c1*c2 + h*q31*wy1^2*c1*c2 + 2*q21*wy1*wz1*c3 + 2*q31*wz1^2*c3)/(2*c1^3);
    jac_block(6,10) = (h*q11*wy1^2*wz1*c3 + h*q11*wz1^3*c3 - h*q21*wz1^2*c1*c2 + h*q31*wy1*wz1*c1*c2 - 2*q21*wy1^2*c3 - 2*q31*wy1*wz1*c3)/(2*c1^3);
    jac_block(6,18) = 1;
    jac_block(7,5) = -wy1*c3/c1;
    jac_block(7,6) = wz1*c3/c1;
    jac_block(7,7) = -c2;
    jac_block(7,9) = (-h*q01*wy1^2*c1*c2 + h*q11*wy1*wz1*c1*c2 + h*q21*wy1^3*c3 + h*q21*wy1*wz1^2*c3 - 2*q01*wz1^2*c3 - 2*q11*wy1*wz1*c3)/(2*c1^3);
    jac_block(7,10) = (-h*q01*wy1*wz1*c1*c2 + h*q11*wz1^2*c1*c2 + h*q21*wy1^2*wz1*c3 + h*q21*wz1^3*c3 + 2*q01*wy1*wz1*c3 + 2*q11*wy1^2*c3)/(2*c1^3);
    jac_block(7,19) = 1;
    jac_block(8,5) = -wz1*c3/c1;
    jac_block(8,6) = -wy1*c3/c1;
    jac_block(8,8) = -c2;
    jac_block(8,9) = (-h*q01*wy1*wz1*c1*c2 - h*q11*wy1^2*c1*c2 + h*q31*wy1^3*c3 + h*q31*wy1*wz1^2*c3 + 2*q01*wy1*wz1*c3 - 2*q11*wz1^2*c3)/(2*c1^3);
    jac_block(8,10) = (-h*q01*wz1^2*c1*c2 - h*q11*wy1*wz1*c1*c2 + h*q31*wy1^2*wz1*c3 + h*q31*wz1^3*c3 - 2*q01*wy1^2*c3 + 2*q11*wy1*wz1*c3)/(2*c1^3);
    jac_block(8,20) = 1;
  else
    jac_block(5,5) = -c2;
    jac_block(5,17) = 1;
    jac_block(6,6) = -c2;
    jac_block(6,18) = 1;
    jac_block(7,7) = -c2;
    jac_block(7,19) = 1;
    jac_block(8,8) = -c2;
    jac_block(8,20) = 1;
  end


  jac_block(9,3) = A*c_L1*s1*cos(sig1)*rho_part_h1/(4*m1);
  jac_block(9,4) = (A*c_L1*s1^2*rho1*cos(sig1) + 2*g*m1*(2*q11^2 + 2*q21^2 - 1))/(4*m1*s1^2);
  jac_block(9,6) = -2*g*q11/s1;
  jac_block(9,7) = -2*g*q21/s1;
  jac_block(9,9) = 1;
  jac_block(9,11) = A*s1*rho1*cos(sig1)/(4*m1);
  jac_block(9,12) = -A*c_L1*s1*rho1*sin(sig1)/(4*m1);
  jac_block(9,15) = A*c_L2*s2*cos(sig2)*rho_part_h2/(4*m2);
  jac_block(9,16) = (A*c_L2*s2^2*rho2*cos(sig2) + 2*g*m2*(2*q12^2 + 2*q22^2 - 1))/(4*m2*s2^2);
  jac_block(9,18) = -2*g*q12/s2;
  jac_block(9,19) = -2*g*q22/s2;
  jac_block(9,23) = A*s2*rho2*cos(sig2)/(4*m2);
  jac_block(9,24) = -A*c_L2*s2*rho2*sin(sig2)/(4*m2);

  jac_block(10,3) = -A*c_L1*s1*sin(sig1)*rho_part_h1/(4*m1);
  jac_block(10,4) = (-A*c_L1*s1^2*rho1*sin(sig1)/4 + g*m1*(q01*q11 + q21*q31))/(m1*s1^2);
  jac_block(10,5) = -g*q11/s1;
  jac_block(10,6) = -g*q01/s1;
  jac_block(10,7) = -g*q31/s1;
  jac_block(10,8) = -g*q21/s1;
  jac_block(10,10) = 1;
  jac_block(10,11) = -A*s1*rho1*sin(sig1)/(4*m1);
  jac_block(10,12) = -A*c_L1*s1*rho1*cos(sig1)/(4*m1);
  jac_block(10,15) = -A*c_L2*s2*sin(sig2)*rho_part_h2/(4*m2);
  jac_block(10,16) = (-A*c_L2*s2^2*rho2*sin(sig2)/4 + g*m2*(q02*q12 + q22*q32))/(m2*s2^2);
  jac_block(10,17) = -g*q12/s2;
  jac_block(10,18) = -g*q02/s2;
  jac_block(10,19) = -g*q32/s2;
  jac_block(10,20) = -g*q22/s2;
  jac_block(10,23) = -A*s2*rho2*sin(sig2)/(4*m2);
  jac_block(10,24) = -A*c_L2*s2*rho2*cos(sig2)/(4*m2);

end

function jac_block = get_end_knot_jac(z_knot, knot_params)
  jac_block = [];

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
      jac_block = get_end_knot_jac(z_knot, params(i));
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
function hes_block = get_first_knot_hes(z_knot, t_knot, knot_params, lambda_knot, z_next, t_next, lb, ub)
  hes_block = zeros(12, 12);
  h = t_next - t_knot;
  pn1 = z_knot(1);
  pe1 = z_knot(2);
  pd1 = z_knot(3);
  s1 = z_knot(4);
  q01 = z_knot(5);
  q11 = z_knot(6);
  q21 = z_knot(7);
  q31 = z_knot(8);
  wy1 = z_knot(9);
  wz1 = z_knot(10);
  c_L1 = z_knot(11);
  sig1 = z_knot(12);

  wy2 = z_next(9);
  wz2 = z_next(10);

  lb_c_L = lb(1);
  lb_sig = lb(2);
  ub_c_L = ub(2);
  ub_sig = ub(3);

  g = knot_params.g;
  m1 = knot_params.m;
  c_D1 = knot_params.c_D;
  rho1 = knot_params.rho;
  rho_part_h = knot_params.rho_part_h;
  c_D_part_s = knot_params.c_D_part_s;
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
  lambda9 = lambda_knot(9);
  lambda10 = lambda_knot(10);
  lambda11 = lambda_knot(11);
  lambda12 = lambda_knot(12);
  lambda13 = lambda_knot(13);
  lambda14 = lambda_knot(14);
  lambda15 = lambda_knot(15);
  lambda16 = lambda_knot(16);
  lambda17 = lambda_knot(17);
  lambda18 = lambda_knot(18);

  hes_block(3,3) = A*s1*(-c_L1*pi_AR_e*(lambda10*sin(sig1) - lambda9*cos(sig1)) + h*lambda4*s1*(c_L1^2 + pi_AR_e*c_D1))*rho_part_h2/(4*m1*pi_AR_e);
  hes_block(3,4) = A*(-c_L1*pi_AR_e*(lambda10*sin(sig1) - lambda9*cos(sig1)) + h*lambda4*pi_AR_e*s1^2*c_D_part_s + 2*h*lambda4*s1*(c_L1^2+ pi_AR_e*c_D1))*rho_part_h/(4*m1*pi_AR_e);
  hes_block(3,11) = A*s1*(2*c_L1*h*lambda4*s1 - pi_AR_e*(lambda10*sin(sig1) - lambda9*cos(sig1)))*rho_part_h/(4*m1*pi_AR_e);
  hes_block(3,12) = -A*c_L1*s1*(lambda10*cos(sig1) + lambda9*sin(sig1))*rho_part_h/(4*m1);
  hes_block(4,3) = hes_block(3,4);
  hes_block(4,4) = A*c_L1^2*h*lambda4*rho1/(2*m1*pi_AR_e) + A*h*lambda4*s1^2*rho1*c_D_part_s2/(4*m1) + A*h*lambda4*s1*rho1*c_D_part_s/m1 + A*h*lambda4*c_D1*rho1/(2*m1) - 2*g*lambda10*q01*q11/s1^3 - 2*g*lambda10*q21*q31/s1^3 - 2*g*lambda9*q11^2/s1^3- 2*g*lambda9*q21^2/s1^3 + g*lambda9/s1^3;
  hes_block(4,5) = g*lambda10*q11/s1^2 - h*lambda2*q31 + h*lambda3*q21;
  hes_block(4,6) = (g*(lambda10*q01 + 2*lambda9*q11) - h*s1^2*(lambda2*q21 + lambda3*q31))/s1^2;
  hes_block(4,7) = (g*(lambda10*q31 + 2*lambda9*q21) + h*s1^2*(2*lambda1*q21 - lambda2*q11 + lambda3*q01))/s1^2;
  hes_block(4,8) = (g*lambda10*q21 - h*s1^2*(-2*lambda1*q31 + lambda2*q01 + lambda3*q11))/s1^2;
  hes_block(4,11) = A*(4*c_L1*h*lambda4*s1 - pi_AR_e*(lambda10*sin(sig1) - lambda9*cos(sig1)))*rho1/(4*m1*pi_AR_e);
  hes_block(4,12) = -A*c_L1*(lambda10*cos(sig1) + lambda9*sin(sig1))*rho1/(4*m1);

  c1 = sqrt(wy1^2 + wz1^2);
  c2 = cos(h*c1/2);
  c3 = sin(h*c1/2);

  hes_block(5,4) = hes_block(4,5);
  hes_block(5,6) = -g*lambda10/s1;
  hes_block(5,7) = h*(g*lambda4 + lambda3*s1);
  hes_block(5,8) = -h*lambda2*s1;
  hes_block(5,9) = (h*lambda5*wy1^3*c3 + h*lambda5*wy1*wz1^2*c3 - h*lambda7*wy1^2*c1*c2 - h*lambda8*wy1*wz1*c1*c2 - 2*lambda7*wz1^2*c3 + 2*lambda8*wy1*wz1*c3)/(2*c1^3);
  hes_block(5,10) = (h*lambda5*wy1^2*wz1*c3 + h*lambda5*wz1^3*c3 - h*lambda7*wy1*wz1*c1*c2 - h*lambda8*wz1^2*c1*c2 + 2*lambda7*wy1*wz1*c3 - 2*lambda8*wy1^2*c3)/(2*c1^3);
  hes_block(6,4) = hes_block(4,6);
  hes_block(6,5) = hes_block(5,6);
  hes_block(6,6) = -2*g*lambda9/s1;
  hes_block(6,7) = -h*lambda2*s1;
  hes_block(6,8) = h*(-g*lambda4 - lambda3*s1);
  hes_block(6,9) = (h*lambda6*wy1^3*c3 + h*lambda6*wy1*wz1^2*c3 + h*lambda7*wy1*wz1*c1*c2 - h*lambda8*wy1^2*c1*c2 - 2*lambda7*wy1*wz1*c3 - 2*lambda8*wz1^2*c3)/(2*c1^3);
  hes_block(6,10) = (h*lambda6*wy1^2*wz1*c3 + h*lambda6*wz1^3*c3 + h*lambda7*wz1^2*c1*c2 - h*lambda8*wy1*wz1*c1*c2 + 2*lambda7*wy1^2*c3 + 2*lambda8*wy1*wz1*c3)/(2*c1^3);
  hes_block(7,4) = hes_block(4,7);
  hes_block(7,5) = hes_block(5,7);
  hes_block(7,6) = hes_block(6,7);
  hes_block(7,7) = -2*g*lambda9/s1 + 2*h*lambda1*s1;
  hes_block(7,8) = -g*lambda10/s1;
  hes_block(7,9) = (h*lambda5*wy1^2*c1*c2 - h*lambda6*wy1*wz1*c1*c2 + h*lambda7*wy1^3*c3 + h*lambda7*wy1*wz1^2*c3 + 2*lambda5*wz1^2*c3 + 2*lambda6*wy1*wz1*c3)/(2*c1^3);
  hes_block(7,10) = (h*lambda5*wy1*wz1*c1*c2 - h*lambda6*wz1^2*c1*c2 + h*lambda7*wy1^2*wz1*c3 + h*lambda7*wz1^3*c3 - 2*lambda5*wy1*wz1*c3 - 2*lambda6*wy1^2*c3)/(2*c1^3);
  hes_block(8,4) = hes_block(4,8);
  hes_block(8,5) = hes_block(5,8);
  hes_block(8,6) = hes_block(6,8);
  hes_block(8,7) = hes_block(7,8);
  hes_block(8,8) = 2*h*lambda1*s1;
  hes_block(8,9) = (h*lambda5*wy1*wz1*c1*c2 + h*lambda6*wy1^2*c1*c2 + h*lambda8*wy1^3*c3 + h*lambda8*wy1*wz1^2*c3 - 2*lambda5*wy1*wz1*c3 + 2*lambda6*wz1^2*c3)/(2*c1^3);
  hes_block(8,10) = (h*lambda5*wz1^2*c1*c2 + h*lambda6*wy1*wz1*c1*c2 + h*lambda8*wy1^2*wz1*c3 + h*lambda8*wz1^3*c3 + 2*lambda5*wy1^2*c3 - 2*lambda6*wy1*wz1*c3)/(2*c1^3);

  hes_block(9,5) = hes_block(5,9);
  hes_block(9,6) = hes_block(6,9);
  hes_block(9,7) = hes_block(7,9);
  hes_block(9,8) = hes_block(8,9);
  hes_block(9,9) = (-12*wy1*(lambda5*(h*wy1*(wy1^2 + wz1^2)^2*(q21*wy1 + q31*wz1)*c2 - 2*wy1*c1^3*(q21*wy1 + q31*wz1)*c3 + c1^5*(h*q01*wy1 + 2*q21)*c3) + lambda6*(-h*wy1*(wy1^2 + wz1^2)^2*(q21*wz1 - q31*wy1)*c2 + 2*wy1*c1^3*(q21*wz1 - q31*wy1)*c3 +c1^5*(h*q11*wy1 + 2*q31)*c3) + lambda7*(-h*wy1*(wy1^2 + wz1^2)^2*(q01*wy1 - q11*wz1)*c2 + 2*wy1*c1^3*(q01*wy1 - q11*wz1)*c3 + c1^5*(h*q21*wy1 - 2*q01)*c3) + lambda8*(-h*wy1*(wy1^2 + wz1^2)^2*(q01*wz1 + q11*wy1)*c2 + 2*wy1*c1^3*(q01*wz1 + q11*wy1)*c3 + c1^5*(h*q31*wy1 - 2*q11)*c3) + 2*lambda9*(wy1^2 + wz1^2)^3) + (wy1^2 + wz1^2)*(lambda5*(-h^2*wy1^2*c1^3*(q21*wy1 + q31*wz1)*c3 + 2*h*q01*c1^5*c3 + 2*h*q21*wy1*(wy1^2 + wz1^2)^2*c2 + 6*h*wy1^2*(wy1^2 + wz1^2)*(q21*wy1 + q31*wz1)*c2 + h*wy1*(wy1^2 + wz1^2)^2*(h*q01*wy1 + 2*q21)*c2 + 2*h*(wy1^2 + wz1^2)^2*(q21*wy1 + q31*wz1)*c2 - 4*q21*wy1*c1^3*c3 - 12*wy1^2*c1*(q21*wy1 + q31*wz1)*c3 + 10*wy1*c1^3*(h*q01*wy1 + 2*q21)*c3 - 4*c1^3*(q21*wy1 + q31*wz1)*c3) + lambda6*(h^2*wy1^2*c1^3*(q21*wz1 - q31*wy1)*c3 + 2*h*q11*c1^5*c3 + 2*h*q31*wy1*(wy1^2 + wz1^2)^2*c2 - 6*h*wy1^2*(wy1^2 + wz1^2)*(q21*wz1 - q31*wy1)*c2 + h*wy1*(wy1^2 + wz1^2)^2*(h*q11*wy1 + 2*q31)*c2 - 2*h*(wy1^2 + wz1^2)^2*(q21*wz1 - q31*wy1)*c2 - 4*q31*wy1*c1^3*c3 + 12*wy1^2*c1*(q21*wz1 - q31*wy1)*c3 + 10*wy1*c1^3*(h*q11*wy1 + 2*q31)*c3 + 4*c1^3*(q21*wz1 - q31*wy1)*c3) + lambda7*(h^2*wy1^2*c1^3*(q01*wy1 - q11*wz1)*c3 - 2*h*q01*wy1*(wy1^2 + wz1^2)^2*c2 + 2*h*q21*c1^5*c3 - 6*h*wy1^2*(wy1^2 + wz1^2)*(q01*wy1 - q11*wz1)*c2 + h*wy1*(wy1^2 + wz1^2)^2*(h*q21*wy1 - 2*q01)*c2 - 2*h*(wy1^2 + wz1^2)^2*(q01*wy1 - q11*wz1)*c2 + 4*q01*wy1*c1^3*c3 + 12*wy1^2*c1*(q01*wy1 - q11*wz1)*c3 + 10*wy1*c1^3*(h*q21*wy1 - 2*q01)*c3 + 4*c1^3*(q01*wy1 - q11*wz1)*c3) + lambda8*(h^2*wy1^2*c1^3*(q01*wz1 + q11*wy1)*c3 - 2*h*q11*wy1*(wy1^2 + wz1^2)^2*c2 + 2*h*q31*c1^5*c3 - 6*h*wy1^2*(wy1^2 + wz1^2)*(q01*wz1 + q11*wy1)*c2 +h*wy1*(wy1^2 + wz1^2)^2*(h*q31*wy1 - 2*q11)*c2 - 2*h*(wy1^2 + wz1^2)^2*(q01*wz1 + q11*wy1)*c2 + 4*q11*wy1*c1^3*c3 + 12*wy1^2*c1*(q01*wz1 + q11*wy1)*c3 + 10*wy1*c1^3*(h*q31*wy1 - 2*q11)*c3 + 4*c1^3*(q01*wz1 + q11*wy1)*c3) + 24*lambda9*wy1*(wy1^2 + wz1^2)^2))/(4*(wy1^2 + wz1^2)^4);
  hes_block(9,10) = (-12*wz1*(lambda5*(h*wy1*(wy1^2 + wz1^2)^2*(q21*wy1 + q31*wz1)*c2 - 2*wy1*c1^3*(q21*wy1 + q31*wz1)*c3 + c1^5*(h*q01*wy1 + 2*q21)*c3) + lambda6*(-h*wy1*(wy1^2 + wz1^2)^2*(q21*wz1 - q31*wy1)*c2 + 2*wy1*c1^3*(q21*wz1 - q31*wy1)*c3 + c1^5*(h*q11*wy1 + 2*q31)*c3) + lambda7*(-h*wy1*(wy1^2 + wz1^2)^2*(q01*wy1 - q11*wz1)*c2 + 2*wy1*c1^3*(q01*wy1 - q11*wz1)*c3 + c1^5*(h*q21*wy1 - 2*q01)*c3) + lambda8*(-h*wy1*(wy1^2 + wz1^2)^2*(q01*wz1 + q11*wy1)*c2 + 2*wy1*c1^3*(q01*wz1 + q11*wy1)*c3 + c1^5*(h*q31*wy1 - 2*q11)*c3) + 2*lambda9*(wy1^2+ wz1^2)^3) + (wy1^2 + wz1^2)*(lambda5*(-h^2*wy1*wz1*c1^3*(q21*wy1 + q31*wz1)*c3 + 2*h*q31*wy1*(wy1^2 + wz1^2)^2*c2 + 6*h*wy1*wz1*(wy1^2 + wz1^2)*(q21*wy1 + q31*wz1)*c2 + h*wz1*(wy1^2 + wz1^2)^2*(h*q01*wy1 + 2*q21)*c2 - 4*q31*wy1*c1^3*c3 - 12*wy1*wz1*c1*(q21*wy1 + q31*wz1)*c3 + 10*wz1*c1^3*(h*q01*wy1 + 2*q21)*c3) + lambda6*(h^2*wy1*wz1*c1^3*(q21*wz1 - q31*wy1)*c3 - 2*h*q21*wy1*(wy1^2 + wz1^2)^2*c2 - 6*h*wy1*wz1*(wy1^2 + wz1^2)*(q21*wz1 - q31*wy1)*c2 + h*wz1*(wy1^2 + wz1^2)^2*(h*q11*wy1 + 2*q31)*c2 + 4*q21*wy1*c1^3*c3 + 12*wy1*wz1*c1*(q21*wz1 - q31*wy1)*c3 + 10*wz1*c1^3*(h*q11*wy1 + 2*q31)*c3) + lambda7*(h^2*wy1*wz1*c1^3*(q01*wy1 - q11*wz1)*c3 + 2*h*q11*wy1*(wy1^2 + wz1^2)^2*c2 - 6*h*wy1*wz1*(wy1^2 + wz1^2)*(q01*wy1 - q11*wz1)*c2 + h*wz1*(wy1^2 + wz1^2)^2*(h*q21*wy1 - 2*q01)*c2 - 4*q11*wy1*c1^3*c3 + 12*wy1*wz1*c1*(q01*wy1 - q11*wz1)*c3 + 10*wz1*c1^3*(h*q21*wy1 - 2*q01)*c3) + lambda8*(h^2*wy1*wz1*c1^3*(q01*wz1 + q11*wy1)*c3 - 2*h*q01*wy1*(wy1^2 + wz1^2)^2*c2 - 6*h*wy1*wz1*(wy1^2 + wz1^2)*(q01*wz1 + q11*wy1)*c2 + h*wz1*(wy1^2 + wz1^2)^2*(h*q31*wy1 - 2*q11)*c2 + 4*q01*wy1*c1^3*c3 + 12*wy1*wz1*c1*(q01*wz1 + q11*wy1)*c3 + 10*wz1*c1^3*(h*q31*wy1 - 2*q11)*c3) + 24*lambda9*wz1*(wy1^2 + wz1^2)^2))/(4*(wy1^2 + wz1^2)^4);

  hes_block(10,5) = hes_block(5,10);
  hes_block(10,6) = hes_block(6,10);
  hes_block(10,7) = hes_block(7,10);
  hes_block(10,8) = hes_block(8,10);
  hes_block(10,9) = hes_block(9,10);
  hes_block(10,10) = (-12*wz1*(2*lambda10*(wy1^2 + wz1^2)^3 + lambda5*(h*wz1*(wy1^2 + wz1^2)^2*(q21*wy1 + q31*wz1)*c2 - 2*wz1*c1^3*(q21*wy1 + q31*wz1)*c3 + c1^5*(h*q01*wz1 + 2*q31)*c3) + lambda6*(-h*wz1*(wy1^2 + wz1^2)^2*(q21*wz1 - q31*wy1)*c2 + 2*wz1*c1^3*(q21*wz1 - q31*wy1)*c3 + c1^5*(h*q11*wz1 - 2*q21)*c3) + lambda7*(-h*wz1*(wy1^2 + wz1^2)^2*(q01*wy1 - q11*wz1)*c2 + 2*wz1*c1^3*(q01*wy1 - q11*wz1)*c3 + c1^5*(h*q21*wz1 + 2*q11)*c3) + lambda8*(-h*wz1*(wy1^2 + wz1^2)^2*(q01*wz1 + q11*wy1)*c2 + 2*wz1*c1^3*(q01*wz1 + q11*wy1)*c3 + c1^5*(h*q31*wz1 - 2*q01)*c3)) + (wy1^2 + wz1^2)*(24*lambda10*wz1*(wy1^2 + wz1^2)^2 + lambda5*(-h^2*wz1^2*c1^3*(q21*wy1 + q31*wz1)*c3 + 2*h*q01*c1^5*c3 + 2*h*q31*wz1*(wy1^2 + wz1^2)^2*c2 + 6*h*wz1^2*(wy1^2 + wz1^2)*(q21*wy1 + q31*wz1)*c2 + h*wz1*(wy1^2 + wz1^2)^2*(h*q01*wz1 + 2*q31)*c2 + 2*h*(wy1^2 + wz1^2)^2*(q21*wy1 + q31*wz1)*c2 - 4*q31*wz1*c1^3*c3 - 12*wz1^2*c1*(q21*wy1 + q31*wz1)*c3 + 10*wz1*c1^3*(h*q01*wz1 + 2*q31)*c3 - 4*c1^3*(q21*wy1 + q31*wz1)*c3) + lambda6*(h^2*wz1^2*c1^3*(q21*wz1 - q31*wy1)*c3 + 2*h*q11*c1^5*c3 - 2*h*q21*wz1*(wy1^2 + wz1^2)^2*c2 - 6*h*wz1^2*(wy1^2 + wz1^2)*(q21*wz1 - q31*wy1)*c2 + h*wz1*(wy1^2 + wz1^2)^2*(h*q11*wz1 - 2*q21)*c2 - 2*h*(wy1^2 + wz1^2)^2*(q21*wz1 - q31*wy1)*c2 + 4*q21*wz1*c1^3*c3 + 12*wz1^2*c1*(q21*wz1 - q31*wy1)*c3 + 10*wz1*c1^3*(h*q11*wz1 - 2*q21)*c3 + 4*c1^3*(q21*wz1 - q31*wy1)*c3) + lambda7*(h^2*wz1^2*c1^3*(q01*wy1 - q11*wz1)*c3 + 2*h*q11*wz1*(wy1^2 + wz1^2)^2*c2 + 2*h*q21*c1^5*c3 - 6*h*wz1^2*(wy1^2 + wz1^2)*(q01*wy1 - q11*wz1)*c2 + h*wz1*(wy1^2 + wz1^2)^2*(h*q21*wz1 + 2*q11)*c2 - 2*h*(wy1^2 + wz1^2)^2*(q01*wy1 - q11*wz1)*c2 - 4*q11*wz1*c1^3*c3 + 12*wz1^2*c1*(q01*wy1 - q11*wz1)*c3 + 10*wz1*c1^3*(h*q21*wz1 + 2*q11)*c3 + 4*c1^3*(q01*wy1 - q11*wz1)*c3) + lambda8*(h^2*wz1^2*c1^3*(q01*wz1 + q11*wy1)*c3 - 2*h*q01*wz1*(wy1^2 + wz1^2)^2*c2 + 2*h*q31*c1^5*c3 - 6*h*wz1^2*(wy1^2 + wz1^2)*(q01*wz1 +q11*wy1)*c2 + h*wz1*(wy1^2 + wz1^2)^2*(h*q31*wz1 - 2*q01)*c2 - 2*h*(wy1^2 + wz1^2)^2*(q01*wz1 + q11*wy1)*c2 + 4*q01*wz1*c1^3*c3 + 12*wz1^2*c1*(q01*wz1 + q11*wy1)*c3 + 10*wz1*c1^3*(h*q31*wz1 - 2*q01)*c3 + 4*c1^3*(q01*wz1 + q11*wy1)*c3)))/(4*(wy1^2 + wz1^2)^4);

  hes_block(11,3) = hes_block(3,11);
  hes_block(11,4) = hes_block(4,11);
  hes_block(11,11) = (A*c_L1^4*h*lambda4*s1^2*rho1/2 - A*c_L1^3*h*lambda4*lb_c_L*s1^2*rho1 - A*c_L1^3*h*lambda4*s1^2*ub_c_L*rho1 + A*c_L1^2*h*lambda4*lb_c_L^2*s1^2*rho1/2 + 2*A*c_L1^2*h*lambda4*lb_c_L*s1^2*ub_c_L*rho1 + A*c_L1^2*h*lambda4*s1^2*ub_c_L^2*rho1/2 - A*c_L1*h*lambda4*lb_c_L^2*s1^2*ub_c_L*rho1 - A*c_L1*h*lambda4*lb_c_L*s1^2*ub_c_L^2*rho1 + A*h*lambda4*lb_c_L^2*s1^2*ub_c_L^2*rho1/2 + 2*c_L1^2*m1*mu*pi_AR_e - 2*c_L1*lb_c_L*m1*mu*pi_AR_e - 2*c_L1*m1*mu*pi_AR_e*ub_c_L + lb_c_L^2*m1*mu*pi_AR_e + m1*mu*pi_AR_e*ub_c_L^2)/(m1*pi_AR_e*(c_L1 - lb_c_L)^2*(c_L1 - ub_c_L)^2);
  hes_block(11,12) = -A*s1*(lambda10*cos(sig1) + lambda9*sin(sig1))*rho1/(4*m1);
  hes_block(12,3) = hes_block(3,12);
  hes_block(12,4) = hes_block(4,12);
  hes_block(12,11) = hes_block(11,12);
  hes_block(12,12) = ((lb_sig - sig1)*(sig1 - ub_sig)*(A*c_L1*s1*(lb_sig - sig1)*(sig1 - ub_sig)*(lambda10*sin(sig1) - lambda9*cos(sig1))*rho1 - A*c_L1*s1*(lb_sig - sig1)*(lambda10*cos(sig1) + lambda9*sin(sig1))*rho1 + A*c_L1*s1*(sig1 - ub_sig)*(lambda10*cos(sig1) + lambda9*sin(sig1))*rho1 + 8*m1*mu) + (lb_sig - sig1)*(A*c_L1*s1*(lb_sig - sig1)*(sig1 - ub_sig)*(lambda10*cos(sig1) + lambda9*sin(sig1))*rho1 + 4*m1*mu*(lb_sig - 2*sig1 + ub_sig)) - (sig1 - ub_sig)*(A*c_L1*s1*(lb_sig - sig1)*(sig1 - ub_sig)*(lambda10*cos(sig1) + lambda9*sin(sig1))*rho1 + 4*m1*mu*(lb_sig - 2*sig1 + ub_sig)))/(4*m1*(lb_sig - sig1)^2*(sig1 - ub_sig)^2);
end

function hes_block = get_middle_knot_hes(z_knot, t_knot, knot_params, lambda_knot, z_next, t_next, lambda_last, lb, ub)
  hes_block = zeros(12, 12);
  h = t_next - t_knot;
  pn1 = z_knot(1);
  pe1 = z_knot(2);
  pd1 = z_knot(3);
  s1 = z_knot(4);
  q01 = z_knot(5);
  q11 = z_knot(6);
  q21 = z_knot(7);
  q31 = z_knot(8);
  wy1 = z_knot(9);
  wz1 = z_knot(10);
  c_L1 = z_knot(11);
  sig1 = z_knot(12);

  wy2 = z_next(9);
  wz2 = z_next(10);

  lb_c_L = lb(1);
  lb_sig = lb(2);
  ub_pd = ub(1);
  ub_c_L = ub(2);
  ub_sig = ub(3);

  g = knot_params.g;
  m1 = knot_params.m;
  c_D1 = knot_params.c_D;
  rho1 = knot_params.rho;
  rho_part_h = knot_params.rho_part_h;
  c_D_part_s = knot_params.c_D_part_s;
  rho_part_h2 = knot_params.rho_part_h2;
  c_D_part_s2 = knot_params.c_D_part_s2;
  lambda1 = lambda_knot(1) + lambda_last(1);
  lambda2 = lambda_knot(2) + lambda_last(2);
  lambda3 = lambda_knot(3) + lambda_last(3);
  lambda4 = lambda_knot(4) + lambda_last(4);
  lambda5 = lambda_knot(5);
  lambda6 = lambda_knot(6);
  lambda7 = lambda_knot(7);
  lambda8 = lambda_knot(8);
  lambda9 = lambda_knot(9) + lambda_last(9);
  lambda10 = lambda_knot(10) + lambda_last(10);

  hes_block(3,3) = A*s1*(-c_L1*pi_AR_e*(lambda10*sin(sig1) - lambda9*cos(sig1)) + h*lambda4*s1*(c_L1^2 + pi_AR_e*c_D1))*rho_part_h2/(4*m1*pi_AR_e) + (mu/(pd1 - ub_pd)^2);
  hes_block(3,4) = A*(-c_L1*pi_AR_e*(lambda10*sin(sig1) - lambda9*cos(sig1)) + h*lambda4*pi_AR_e*s1^2*c_D_part_s + 2*h*lambda4*s1*(c_L1^2+ pi_AR_e*c_D1))*rho_part_h/(4*m1*pi_AR_e);
  hes_block(3,11) = A*s1*(2*c_L1*h*lambda4*s1 - pi_AR_e*(lambda10*sin(sig1) - lambda9*cos(sig1)))*rho_part_h/(4*m1*pi_AR_e);
  hes_block(3,12) = -A*c_L1*s1*(lambda10*cos(sig1) + lambda9*sin(sig1))*rho_part_h/(4*m1);
  hes_block(4,3) = hes_block(3,4);
  hes_block(4,4) = A*c_L1^2*h*lambda4*rho1/(2*m1*pi_AR_e) + A*h*lambda4*s1^2*rho1*c_D_part_s2/(4*m1) + A*h*lambda4*s1*rho1*c_D_part_s/m1 + A*h*lambda4*c_D1*rho1/(2*m1) - 2*g*lambda10*q01*q11/s1^3 - 2*g*lambda10*q21*q31/s1^3 - 2*g*lambda9*q11^2/s1^3- 2*g*lambda9*q21^2/s1^3 + g*lambda9/s1^3;
  hes_block(4,5) = g*lambda10*q11/s1^2 - h*lambda2*q31 + h*lambda3*q21;
  hes_block(4,6) = (g*(lambda10*q01 + 2*lambda9*q11) - h*s1^2*(lambda2*q21 + lambda3*q31))/s1^2;
  hes_block(4,7) = (g*(lambda10*q31 + 2*lambda9*q21) + h*s1^2*(2*lambda1*q21 - lambda2*q11 + lambda3*q01))/s1^2;
  hes_block(4,8) = (g*lambda10*q21 - h*s1^2*(-2*lambda1*q31 + lambda2*q01 + lambda3*q11))/s1^2;
  hes_block(4,11) = A*(4*c_L1*h*lambda4*s1 - pi_AR_e*(lambda10*sin(sig1) - lambda9*cos(sig1)))*rho1/(4*m1*pi_AR_e);
  hes_block(4,12) = -A*c_L1*(lambda10*cos(sig1) + lambda9*sin(sig1))*rho1/(4*m1);

  c1 = sqrt(wy1^2 + wz1^2);
  c2 = cos(h*c1/2);
  c3 = sin(h*c1/2);

  hes_block(5,4) = hes_block(4,5);
  hes_block(5,6) = -g*lambda10/s1;
  hes_block(5,7) = h*(g*lambda4 + lambda3*s1);
  hes_block(5,8) = -h*lambda2*s1;
  hes_block(5,9) = (h*lambda5*wy1^3*c3 + h*lambda5*wy1*wz1^2*c3 - h*lambda7*wy1^2*c1*c2 - h*lambda8*wy1*wz1*c1*c2 - 2*lambda7*wz1^2*c3 + 2*lambda8*wy1*wz1*c3)/(2*c1^3);
  hes_block(5,10) = (h*lambda5*wy1^2*wz1*c3 + h*lambda5*wz1^3*c3 - h*lambda7*wy1*wz1*c1*c2 - h*lambda8*wz1^2*c1*c2 + 2*lambda7*wy1*wz1*c3 - 2*lambda8*wy1^2*c3)/(2*c1^3);
  hes_block(6,4) = hes_block(4,6);
  hes_block(6,5) = hes_block(5,6);
  hes_block(6,6) = -2*g*lambda9/s1;
  hes_block(6,7) = -h*lambda2*s1;
  hes_block(6,8) = h*(-g*lambda4 - lambda3*s1);
  hes_block(6,9) = (h*lambda6*wy1^3*c3 + h*lambda6*wy1*wz1^2*c3 + h*lambda7*wy1*wz1*c1*c2 - h*lambda8*wy1^2*c1*c2 - 2*lambda7*wy1*wz1*c3 - 2*lambda8*wz1^2*c3)/(2*c1^3);
  hes_block(6,10) = (h*lambda6*wy1^2*wz1*c3 + h*lambda6*wz1^3*c3 + h*lambda7*wz1^2*c1*c2 - h*lambda8*wy1*wz1*c1*c2 + 2*lambda7*wy1^2*c3 + 2*lambda8*wy1*wz1*c3)/(2*c1^3);
  hes_block(7,4) = hes_block(4,7);
  hes_block(7,5) = hes_block(5,7);
  hes_block(7,6) = hes_block(6,7);
  hes_block(7,7) = -2*g*lambda9/s1 + 2*h*lambda1*s1;
  hes_block(7,8) = -g*lambda10/s1;
  hes_block(7,9) = (h*lambda5*wy1^2*c1*c2 - h*lambda6*wy1*wz1*c1*c2 + h*lambda7*wy1^3*c3 + h*lambda7*wy1*wz1^2*c3 + 2*lambda5*wz1^2*c3 + 2*lambda6*wy1*wz1*c3)/(2*c1^3);
  hes_block(7,10) = (h*lambda5*wy1*wz1*c1*c2 - h*lambda6*wz1^2*c1*c2 + h*lambda7*wy1^2*wz1*c3 + h*lambda7*wz1^3*c3 - 2*lambda5*wy1*wz1*c3 - 2*lambda6*wy1^2*c3)/(2*c1^3);
  hes_block(8,4) = hes_block(4,8);
  hes_block(8,5) = hes_block(5,8);
  hes_block(8,6) = hes_block(6,8);
  hes_block(8,7) = hes_block(7,8);
  hes_block(8,8) = 2*h*lambda1*s1;
  hes_block(8,9) = (h*lambda5*wy1*wz1*c1*c2 + h*lambda6*wy1^2*c1*c2 + h*lambda8*wy1^3*c3 + h*lambda8*wy1*wz1^2*c3 - 2*lambda5*wy1*wz1*c3 + 2*lambda6*wz1^2*c3)/(2*c1^3);
  hes_block(8,10) = (h*lambda5*wz1^2*c1*c2 + h*lambda6*wy1*wz1*c1*c2 + h*lambda8*wy1^2*wz1*c3 + h*lambda8*wz1^3*c3 + 2*lambda5*wy1^2*c3 - 2*lambda6*wy1*wz1*c3)/(2*c1^3);

  hes_block(9,5) = hes_block(5,9);
  hes_block(9,6) = hes_block(6,9);
  hes_block(9,7) = hes_block(7,9);
  hes_block(9,8) = hes_block(8,9);
  hes_block(9,9) = (-12*wy1*(lambda5*(h*wy1*(wy1^2 + wz1^2)^2*(q21*wy1 + q31*wz1)*c2 - 2*wy1*c1^3*(q21*wy1 + q31*wz1)*c3 + c1^5*(h*q01*wy1 + 2*q21)*c3) + lambda6*(-h*wy1*(wy1^2 + wz1^2)^2*(q21*wz1 - q31*wy1)*c2 + 2*wy1*c1^3*(q21*wz1 - q31*wy1)*c3 +c1^5*(h*q11*wy1 + 2*q31)*c3) + lambda7*(-h*wy1*(wy1^2 + wz1^2)^2*(q01*wy1 - q11*wz1)*c2 + 2*wy1*c1^3*(q01*wy1 - q11*wz1)*c3 + c1^5*(h*q21*wy1 - 2*q01)*c3) + lambda8*(-h*wy1*(wy1^2 + wz1^2)^2*(q01*wz1 + q11*wy1)*c2 + 2*wy1*c1^3*(q01*wz1 + q11*wy1)*c3 + c1^5*(h*q31*wy1 - 2*q11)*c3) + 2*lambda9*(wy1^2 + wz1^2)^3) + (wy1^2 + wz1^2)*(lambda5*(-h^2*wy1^2*c1^3*(q21*wy1 + q31*wz1)*c3 + 2*h*q01*c1^5*c3 + 2*h*q21*wy1*(wy1^2 + wz1^2)^2*c2 + 6*h*wy1^2*(wy1^2 + wz1^2)*(q21*wy1 + q31*wz1)*c2 + h*wy1*(wy1^2 + wz1^2)^2*(h*q01*wy1 + 2*q21)*c2 + 2*h*(wy1^2 + wz1^2)^2*(q21*wy1 + q31*wz1)*c2 - 4*q21*wy1*c1^3*c3 - 12*wy1^2*c1*(q21*wy1 + q31*wz1)*c3 + 10*wy1*c1^3*(h*q01*wy1 + 2*q21)*c3 - 4*c1^3*(q21*wy1 + q31*wz1)*c3) + lambda6*(h^2*wy1^2*c1^3*(q21*wz1 - q31*wy1)*c3 + 2*h*q11*c1^5*c3 + 2*h*q31*wy1*(wy1^2 + wz1^2)^2*c2 - 6*h*wy1^2*(wy1^2 + wz1^2)*(q21*wz1 - q31*wy1)*c2 + h*wy1*(wy1^2 + wz1^2)^2*(h*q11*wy1 + 2*q31)*c2 - 2*h*(wy1^2 + wz1^2)^2*(q21*wz1 - q31*wy1)*c2 - 4*q31*wy1*c1^3*c3 + 12*wy1^2*c1*(q21*wz1 - q31*wy1)*c3 + 10*wy1*c1^3*(h*q11*wy1 + 2*q31)*c3 + 4*c1^3*(q21*wz1 - q31*wy1)*c3) + lambda7*(h^2*wy1^2*c1^3*(q01*wy1 - q11*wz1)*c3 - 2*h*q01*wy1*(wy1^2 + wz1^2)^2*c2 + 2*h*q21*c1^5*c3 - 6*h*wy1^2*(wy1^2 + wz1^2)*(q01*wy1 - q11*wz1)*c2 + h*wy1*(wy1^2 + wz1^2)^2*(h*q21*wy1 - 2*q01)*c2 - 2*h*(wy1^2 + wz1^2)^2*(q01*wy1 - q11*wz1)*c2 + 4*q01*wy1*c1^3*c3 + 12*wy1^2*c1*(q01*wy1 - q11*wz1)*c3 + 10*wy1*c1^3*(h*q21*wy1 - 2*q01)*c3 + 4*c1^3*(q01*wy1 - q11*wz1)*c3) + lambda8*(h^2*wy1^2*c1^3*(q01*wz1 + q11*wy1)*c3 - 2*h*q11*wy1*(wy1^2 + wz1^2)^2*c2 + 2*h*q31*c1^5*c3 - 6*h*wy1^2*(wy1^2 + wz1^2)*(q01*wz1 + q11*wy1)*c2 +h*wy1*(wy1^2 + wz1^2)^2*(h*q31*wy1 - 2*q11)*c2 - 2*h*(wy1^2 + wz1^2)^2*(q01*wz1 + q11*wy1)*c2 + 4*q11*wy1*c1^3*c3 + 12*wy1^2*c1*(q01*wz1 + q11*wy1)*c3 + 10*wy1*c1^3*(h*q31*wy1 - 2*q11)*c3 + 4*c1^3*(q01*wz1 + q11*wy1)*c3) + 24*lambda9*wy1*(wy1^2 + wz1^2)^2))/(4*(wy1^2 + wz1^2)^4);
  hes_block(9,10) = (-12*wz1*(lambda5*(h*wy1*(wy1^2 + wz1^2)^2*(q21*wy1 + q31*wz1)*c2 - 2*wy1*c1^3*(q21*wy1 + q31*wz1)*c3 + c1^5*(h*q01*wy1 + 2*q21)*c3) + lambda6*(-h*wy1*(wy1^2 + wz1^2)^2*(q21*wz1 - q31*wy1)*c2 + 2*wy1*c1^3*(q21*wz1 - q31*wy1)*c3 + c1^5*(h*q11*wy1 + 2*q31)*c3) + lambda7*(-h*wy1*(wy1^2 + wz1^2)^2*(q01*wy1 - q11*wz1)*c2 + 2*wy1*c1^3*(q01*wy1 - q11*wz1)*c3 + c1^5*(h*q21*wy1 - 2*q01)*c3) + lambda8*(-h*wy1*(wy1^2 + wz1^2)^2*(q01*wz1 + q11*wy1)*c2 + 2*wy1*c1^3*(q01*wz1 + q11*wy1)*c3 + c1^5*(h*q31*wy1 - 2*q11)*c3) + 2*lambda9*(wy1^2+ wz1^2)^3) + (wy1^2 + wz1^2)*(lambda5*(-h^2*wy1*wz1*c1^3*(q21*wy1 + q31*wz1)*c3 + 2*h*q31*wy1*(wy1^2 + wz1^2)^2*c2 + 6*h*wy1*wz1*(wy1^2 + wz1^2)*(q21*wy1 + q31*wz1)*c2 + h*wz1*(wy1^2 + wz1^2)^2*(h*q01*wy1 + 2*q21)*c2 - 4*q31*wy1*c1^3*c3 - 12*wy1*wz1*c1*(q21*wy1 + q31*wz1)*c3 + 10*wz1*c1^3*(h*q01*wy1 + 2*q21)*c3) + lambda6*(h^2*wy1*wz1*c1^3*(q21*wz1 - q31*wy1)*c3 - 2*h*q21*wy1*(wy1^2 + wz1^2)^2*c2 - 6*h*wy1*wz1*(wy1^2 + wz1^2)*(q21*wz1 - q31*wy1)*c2 + h*wz1*(wy1^2 + wz1^2)^2*(h*q11*wy1 + 2*q31)*c2 + 4*q21*wy1*c1^3*c3 + 12*wy1*wz1*c1*(q21*wz1 - q31*wy1)*c3 + 10*wz1*c1^3*(h*q11*wy1 + 2*q31)*c3) + lambda7*(h^2*wy1*wz1*c1^3*(q01*wy1 - q11*wz1)*c3 + 2*h*q11*wy1*(wy1^2 + wz1^2)^2*c2 - 6*h*wy1*wz1*(wy1^2 + wz1^2)*(q01*wy1 - q11*wz1)*c2 + h*wz1*(wy1^2 + wz1^2)^2*(h*q21*wy1 - 2*q01)*c2 - 4*q11*wy1*c1^3*c3 + 12*wy1*wz1*c1*(q01*wy1 - q11*wz1)*c3 + 10*wz1*c1^3*(h*q21*wy1 - 2*q01)*c3) + lambda8*(h^2*wy1*wz1*c1^3*(q01*wz1 + q11*wy1)*c3 - 2*h*q01*wy1*(wy1^2 + wz1^2)^2*c2 - 6*h*wy1*wz1*(wy1^2 + wz1^2)*(q01*wz1 + q11*wy1)*c2 + h*wz1*(wy1^2 + wz1^2)^2*(h*q31*wy1 - 2*q11)*c2 + 4*q01*wy1*c1^3*c3 + 12*wy1*wz1*c1*(q01*wz1 + q11*wy1)*c3 + 10*wz1*c1^3*(h*q31*wy1 - 2*q11)*c3) + 24*lambda9*wz1*(wy1^2 + wz1^2)^2))/(4*(wy1^2 + wz1^2)^4);

  hes_block(10,5) = hes_block(5,10);
  hes_block(10,6) = hes_block(6,10);
  hes_block(10,7) = hes_block(7,10);
  hes_block(10,8) = hes_block(8,10);
  hes_block(10,9) = hes_block(9,10);
  hes_block(10,10) = (-12*wz1*(2*lambda10*(wy1^2 + wz1^2)^3 + lambda5*(h*wz1*(wy1^2 + wz1^2)^2*(q21*wy1 + q31*wz1)*c2 - 2*wz1*c1^3*(q21*wy1 + q31*wz1)*c3 + c1^5*(h*q01*wz1 + 2*q31)*c3) + lambda6*(-h*wz1*(wy1^2 + wz1^2)^2*(q21*wz1 - q31*wy1)*c2 + 2*wz1*c1^3*(q21*wz1 - q31*wy1)*c3 + c1^5*(h*q11*wz1 - 2*q21)*c3) + lambda7*(-h*wz1*(wy1^2 + wz1^2)^2*(q01*wy1 - q11*wz1)*c2 + 2*wz1*c1^3*(q01*wy1 - q11*wz1)*c3 + c1^5*(h*q21*wz1 + 2*q11)*c3) + lambda8*(-h*wz1*(wy1^2 + wz1^2)^2*(q01*wz1 + q11*wy1)*c2 + 2*wz1*c1^3*(q01*wz1 + q11*wy1)*c3 + c1^5*(h*q31*wz1 - 2*q01)*c3)) + (wy1^2 + wz1^2)*(24*lambda10*wz1*(wy1^2 + wz1^2)^2 + lambda5*(-h^2*wz1^2*c1^3*(q21*wy1 + q31*wz1)*c3 + 2*h*q01*c1^5*c3 + 2*h*q31*wz1*(wy1^2 + wz1^2)^2*c2 + 6*h*wz1^2*(wy1^2 + wz1^2)*(q21*wy1 + q31*wz1)*c2 + h*wz1*(wy1^2 + wz1^2)^2*(h*q01*wz1 + 2*q31)*c2 + 2*h*(wy1^2 + wz1^2)^2*(q21*wy1 + q31*wz1)*c2 - 4*q31*wz1*c1^3*c3 - 12*wz1^2*c1*(q21*wy1 + q31*wz1)*c3 + 10*wz1*c1^3*(h*q01*wz1 + 2*q31)*c3 - 4*c1^3*(q21*wy1 + q31*wz1)*c3) + lambda6*(h^2*wz1^2*c1^3*(q21*wz1 - q31*wy1)*c3 + 2*h*q11*c1^5*c3 - 2*h*q21*wz1*(wy1^2 + wz1^2)^2*c2 - 6*h*wz1^2*(wy1^2 + wz1^2)*(q21*wz1 - q31*wy1)*c2 + h*wz1*(wy1^2 + wz1^2)^2*(h*q11*wz1 - 2*q21)*c2 - 2*h*(wy1^2 + wz1^2)^2*(q21*wz1 - q31*wy1)*c2 + 4*q21*wz1*c1^3*c3 + 12*wz1^2*c1*(q21*wz1 - q31*wy1)*c3 + 10*wz1*c1^3*(h*q11*wz1 - 2*q21)*c3 + 4*c1^3*(q21*wz1 - q31*wy1)*c3) + lambda7*(h^2*wz1^2*c1^3*(q01*wy1 - q11*wz1)*c3 + 2*h*q11*wz1*(wy1^2 + wz1^2)^2*c2 + 2*h*q21*c1^5*c3 - 6*h*wz1^2*(wy1^2 + wz1^2)*(q01*wy1 - q11*wz1)*c2 + h*wz1*(wy1^2 + wz1^2)^2*(h*q21*wz1 + 2*q11)*c2 - 2*h*(wy1^2 + wz1^2)^2*(q01*wy1 - q11*wz1)*c2 - 4*q11*wz1*c1^3*c3 + 12*wz1^2*c1*(q01*wy1 - q11*wz1)*c3 + 10*wz1*c1^3*(h*q21*wz1 + 2*q11)*c3 + 4*c1^3*(q01*wy1 - q11*wz1)*c3) + lambda8*(h^2*wz1^2*c1^3*(q01*wz1 + q11*wy1)*c3 - 2*h*q01*wz1*(wy1^2 + wz1^2)^2*c2 + 2*h*q31*c1^5*c3 - 6*h*wz1^2*(wy1^2 + wz1^2)*(q01*wz1 +q11*wy1)*c2 + h*wz1*(wy1^2 + wz1^2)^2*(h*q31*wz1 - 2*q01)*c2 - 2*h*(wy1^2 + wz1^2)^2*(q01*wz1 + q11*wy1)*c2 + 4*q01*wz1*c1^3*c3 + 12*wz1^2*c1*(q01*wz1 + q11*wy1)*c3 + 10*wz1*c1^3*(h*q31*wz1 - 2*q01)*c3 + 4*c1^3*(q01*wz1 + q11*wy1)*c3)))/(4*(wy1^2 + wz1^2)^4);

  hes_block(11,3) = hes_block(3,11);
  hes_block(11,4) = hes_block(4,11);
  hes_block(11,11) = (A*c_L1^4*h*lambda4*s1^2*rho1/2 - A*c_L1^3*h*lambda4*lb_c_L*s1^2*rho1 - A*c_L1^3*h*lambda4*s1^2*ub_c_L*rho1 + A*c_L1^2*h*lambda4*lb_c_L^2*s1^2*rho1/2 + 2*A*c_L1^2*h*lambda4*lb_c_L*s1^2*ub_c_L*rho1 + A*c_L1^2*h*lambda4*s1^2*ub_c_L^2*rho1/2 - A*c_L1*h*lambda4*lb_c_L^2*s1^2*ub_c_L*rho1 - A*c_L1*h*lambda4*lb_c_L*s1^2*ub_c_L^2*rho1 + A*h*lambda4*lb_c_L^2*s1^2*ub_c_L^2*rho1/2 + 2*c_L1^2*m1*mu*pi_AR_e - 2*c_L1*lb_c_L*m1*mu*pi_AR_e - 2*c_L1*m1*mu*pi_AR_e*ub_c_L + lb_c_L^2*m1*mu*pi_AR_e + m1*mu*pi_AR_e*ub_c_L^2)/(m1*pi_AR_e*(c_L1 - lb_c_L)^2*(c_L1 - ub_c_L)^2);
  hes_block(11,12) = -A*s1*(lambda10*cos(sig1) + lambda9*sin(sig1))*rho1/(4*m1);
  hes_block(12,3) = hes_block(3,12);
  hes_block(12,4) = hes_block(4,12);
  hes_block(12,11) = hes_block(11,12);
  hes_block(12,12) = ((lb_sig - sig1)*(sig1 - ub_sig)*(A*c_L1*s1*(lb_sig - sig1)*(sig1 - ub_sig)*(lambda10*sin(sig1) - lambda9*cos(sig1))*rho1 - A*c_L1*s1*(lb_sig - sig1)*(lambda10*cos(sig1) + lambda9*sin(sig1))*rho1 + A*c_L1*s1*(sig1 - ub_sig)*(lambda10*cos(sig1) + lambda9*sin(sig1))*rho1 + 8*m1*mu) + (lb_sig - sig1)*(A*c_L1*s1*(lb_sig - sig1)*(sig1 - ub_sig)*(lambda10*cos(sig1) + lambda9*sin(sig1))*rho1 + 4*m1*mu*(lb_sig - 2*sig1 + ub_sig)) - (sig1 - ub_sig)*(A*c_L1*s1*(lb_sig - sig1)*(sig1 - ub_sig)*(lambda10*cos(sig1) + lambda9*sin(sig1))*rho1 + 4*m1*mu*(lb_sig - 2*sig1 + ub_sig)))/(4*m1*(lb_sig - sig1)^2*(sig1 - ub_sig)^2);

  hes_block(5,5) = hes_block(5,5) + 2.0*mu_q;
  hes_block(6,6) = hes_block(6,6) + 2.0*mu_q;
  hes_block(7,7) = hes_block(7,7) + 2.0*mu_q;
  hes_block(8,8) = hes_block(8,8) + 2.0*mu_q;

end

function hes_block = get_end_knot_hes(z_knot, t_knot, knot_params, lambda_knot, t_next, lambda_last, lb, ub)
  hes_block = zeros(12, 12);
  h = t_next - t_knot;
  pn1 = z_knot(1);
  pe1 = z_knot(2);
  pd1 = z_knot(3);
  s1 = z_knot(4);
  q01 = z_knot(5);
  q11 = z_knot(6);
  q21 = z_knot(7);
  q31 = z_knot(8);
  wy1 = z_knot(9);
  wz1 = z_knot(10);
  c_L1 = z_knot(11);
  sig1 = z_knot(12);

  lb_c_L = lb(1);
  lb_sig = lb(2);
  ub_pd = ub(1);
  ub_c_L = ub(2);
  ub_sig = ub(3);

  g = knot_params.g;
  m1 = knot_params.m;
  c_D1 = knot_params.c_D;
  rho1 = knot_params.rho;
  rho_part_h = knot_params.rho_part_h;
  c_D_part_s = knot_params.c_D_part_s;
  rho_part_h2 = knot_params.rho_part_h2;
  c_D_part_s2 = knot_params.c_D_part_s2;
  lambda1 = lambda_last(1);
  lambda2 = lambda_last(2);
  lambda3 = lambda_last(3);
  lambda4 = lambda_last(4);
  lambda5 = 0.0;
  lambda6 = 0.0;
  lambda7 = 0.0;
  lambda8 = 0.0;
  lambda9 = lambda_last(9);
  lambda10 = lambda_last(10);

  hes_block(1,1) = 2;

  hes_block(2,2) = 2;

  hes_block(3,3) = 2.0 + A*s1*(-c_L1*pi_AR_e*(lambda10*sin(sig1) - lambda9*cos(sig1)) + h*lambda4*s1*(c_L1^2 + pi_AR_e*c_D1))*rho_part_h2/(4*m1*pi_AR_e) + (mu/(pd1 - ub_pd)^2);
  hes_block(3,4) = A*(-c_L1*pi_AR_e*(lambda10*sin(sig1) - lambda9*cos(sig1)) + h*lambda4*pi_AR_e*s1^2*c_D_part_s + 2*h*lambda4*s1*(c_L1^2+ pi_AR_e*c_D1))*rho_part_h/(4*m1*pi_AR_e);
  hes_block(3,11) = A*s1*(2*c_L1*h*lambda4*s1 - pi_AR_e*(lambda10*sin(sig1) - lambda9*cos(sig1)))*rho_part_h/(4*m1*pi_AR_e);
  hes_block(3,12) = -A*c_L1*s1*(lambda10*cos(sig1) + lambda9*sin(sig1))*rho_part_h/(4*m1);
  hes_block(4,3) = hes_block(3,4);
  hes_block(4,4) = -2.0*q_speed + A*c_L1^2*h*lambda4*rho1/(2*m1*pi_AR_e) + A*h*lambda4*s1^2*rho1*c_D_part_s2/(4*m1) + A*h*lambda4*s1*rho1*c_D_part_s/m1 + A*h*lambda4*c_D1*rho1/(2*m1) - 2*g*lambda10*q01*q11/s1^3 - 2*g*lambda10*q21*q31/s1^3 - 2*g*lambda9*q11^2/s1^3- 2*g*lambda9*q21^2/s1^3 + g*lambda9/s1^3;
  hes_block(4,5) = g*lambda10*q11/s1^2 - h*lambda2*q31 + h*lambda3*q21;
  hes_block(4,6) = (g*(lambda10*q01 + 2*lambda9*q11) - h*s1^2*(lambda2*q21 + lambda3*q31))/s1^2;
  hes_block(4,7) = (g*(lambda10*q31 + 2*lambda9*q21) + h*s1^2*(2*lambda1*q21 - lambda2*q11 + lambda3*q01))/s1^2;
  hes_block(4,8) = (g*lambda10*q21 - h*s1^2*(-2*lambda1*q31 + lambda2*q01 + lambda3*q11))/s1^2;
  hes_block(4,11) = A*(4*c_L1*h*lambda4*s1 - pi_AR_e*(lambda10*sin(sig1) - lambda9*cos(sig1)))*rho1/(4*m1*pi_AR_e);
  hes_block(4,12) = -A*c_L1*(lambda10*cos(sig1) + lambda9*sin(sig1))*rho1/(4*m1);

  c1 = sqrt(wy1^2 + wz1^2);
  c2 = cos(h*c1/2);
  c3 = sin(h*c1/2);

  hes_block(5,4) = hes_block(4,5);
  hes_block(5,6) = -g*lambda10/s1;
  hes_block(5,7) = h*(g*lambda4 + lambda3*s1);
  hes_block(5,8) = -h*lambda2*s1;
  hes_block(6,4) = hes_block(4,6);
  hes_block(6,5) = hes_block(5,6);
  hes_block(6,6) = -2*g*lambda9/s1;
  hes_block(6,7) = -h*lambda2*s1;
  hes_block(6,8) = h*(-g*lambda4 - lambda3*s1);
  hes_block(7,4) = hes_block(4,7);
  hes_block(7,5) = hes_block(5,7);
  hes_block(7,6) = hes_block(6,7);
  hes_block(7,7) = -2*g*lambda9/s1 + 2*h*lambda1*s1;
  hes_block(7,8) = -g*lambda10/s1;
  hes_block(8,4) = hes_block(4,8);
  hes_block(8,5) = hes_block(5,8);
  hes_block(8,6) = hes_block(6,8);
  hes_block(8,7) = hes_block(7,8);
  hes_block(8,8) = 2*h*lambda1*s1;

  hes_block(11,3) = hes_block(3,11);
  hes_block(11,4) = hes_block(4,11);
  hes_block(11,11) = (A*c_L1^4*h*lambda4*s1^2*rho1/2 - A*c_L1^3*h*lambda4*lb_c_L*s1^2*rho1 - A*c_L1^3*h*lambda4*s1^2*ub_c_L*rho1 + A*c_L1^2*h*lambda4*lb_c_L^2*s1^2*rho1/2 + 2*A*c_L1^2*h*lambda4*lb_c_L*s1^2*ub_c_L*rho1 + A*c_L1^2*h*lambda4*s1^2*ub_c_L^2*rho1/2 - A*c_L1*h*lambda4*lb_c_L^2*s1^2*ub_c_L*rho1 - A*c_L1*h*lambda4*lb_c_L*s1^2*ub_c_L^2*rho1 + A*h*lambda4*lb_c_L^2*s1^2*ub_c_L^2*rho1/2 + 2*c_L1^2*m1*mu*pi_AR_e - 2*c_L1*lb_c_L*m1*mu*pi_AR_e - 2*c_L1*m1*mu*pi_AR_e*ub_c_L + lb_c_L^2*m1*mu*pi_AR_e + m1*mu*pi_AR_e*ub_c_L^2)/(m1*pi_AR_e*(c_L1 - lb_c_L)^2*(c_L1 - ub_c_L)^2);
  hes_block(11,12) = -A*s1*(lambda10*cos(sig1) + lambda9*sin(sig1))*rho1/(4*m1);
  hes_block(12,3) = hes_block(3,12);
  hes_block(12,4) = hes_block(4,12);
  hes_block(12,11) = hes_block(11,12);
  hes_block(12,12) = ((lb_sig - sig1)*(sig1 - ub_sig)*(A*c_L1*s1*(lb_sig - sig1)*(sig1 - ub_sig)*(lambda10*sin(sig1) - lambda9*cos(sig1))*rho1 - A*c_L1*s1*(lb_sig - sig1)*(lambda10*cos(sig1) + lambda9*sin(sig1))*rho1 + A*c_L1*s1*(sig1 - ub_sig)*(lambda10*cos(sig1) + lambda9*sin(sig1))*rho1 + 8*m1*mu) + (lb_sig - sig1)*(A*c_L1*s1*(lb_sig - sig1)*(sig1 - ub_sig)*(lambda10*cos(sig1) + lambda9*sin(sig1))*rho1 + 4*m1*mu*(lb_sig - 2*sig1 + ub_sig)) - (sig1 - ub_sig)*(A*c_L1*s1*(lb_sig - sig1)*(sig1 - ub_sig)*(lambda10*cos(sig1) + lambda9*sin(sig1))*rho1 + 4*m1*mu*(lb_sig - 2*sig1 + ub_sig)))/(4*m1*(lb_sig - sig1)^2*(sig1 - ub_sig)^2);

  hes_block(5,5) = hes_block(5,5) + 2.0*mu_q;
  hes_block(6,6) = hes_block(6,6) + 2.0*mu_q;
  hes_block(7,7) = hes_block(7,7) + 2.0*mu_q;
  hes_block(8,8) = hes_block(8,8) + 2.0*mu_q;
end

if option.eval_hes
  hessian = zeros(O.N_decision_variables, O.N_decision_variables);
  con_start = 1;
  for i = 1:O.N_knots
    z_start = O.knot_size*(i-1) + 1;
    z_end = O.knot_size*i;
    z_knot = O.z(z_start:z_end);

    if i == 1
      z_next = O.z(z_end+1:O.knot_size*(i+1));
      con_end = con_start + O.First_knot_constraints - 1;
      lambda_knot = O.lambda(con_start:con_end);
      hes_block = get_first_knot_hes(z_knot, O.t(i), params(i), lambda_knot, z_next, O.t(i+1), O.lb, O.ub);
      con_start = con_end + 1;
    elseif i == O.N_knots
      con_end = con_start + O.End_knot_constraints - 1;
      lambda_knot = O.lambda(con_start:con_end);
      hes_block = get_end_knot_hes(z_knot, O.t(i-1), params(i), lambda_knot, O.t(i), lambda_last, O.lb, O.ub);
    else
      z_next = O.z(z_end+1:O.knot_size*(i+1));
      con_end = con_start + O.Middle_knot_constraints - 1;
      lambda_knot = O.lambda(con_start:con_end);
      hes_block = get_middle_knot_hes(z_knot, O.t(i), params(i), lambda_knot, z_next, O.t(i+1), lambda_last, O.lb, O.ub);
      con_start = con_end + 1;
    end

    hessian(z_start:z_end,z_start:z_end) = hes_block;
    lambda_last = lambda_knot;
  end
end

end
