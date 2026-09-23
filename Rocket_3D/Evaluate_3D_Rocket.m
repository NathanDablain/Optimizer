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
% 3 - p_d,
% 4 - s,
% 5 - q0
% 6 - q1
% 7 - q2
% 8 - q3
% 9 - Cly
% 10 - Clz
% 11 - s_Cly_lb
% 12 - s_Clz_lb
% 13 - s_pd_ub
% 14 - s_Cly_ub
% 15 - s_Clz_ub
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                       Evaluate cost                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  mu = 1.0e-6;
  function contribution = get_first_knot_cost(z_knot)
    contribution = -mu*(log(z_knot(11)) + log(z_knot(12)) + log(z_knot(13)) +...
                        log(z_knot(14)) + log(z_knot(15)));
  end

  function contribution = get_middle_knot_cost(z_knot)
    contribution = -mu*(log(z_knot(11)) + log(z_knot(12)) + log(z_knot(13)) +...
                        log(z_knot(14)) + log(z_knot(15)));
  end

  function contribution = get_end_knot_cost(z_knot, xd)
    contribution = (z_knot(1) - xd(1))^2 +...
                   (z_knot(2) - xd(2))^2 +...
                   (z_knot(3) - xd(3))^2 -...
                    mu*(log(z_knot(11)) + log(z_knot(12)) + log(z_knot(13)) +...
                        log(z_knot(14)) + log(z_knot(15)));
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
  grad_vec(11) = -mu / z_knot(11);
  grad_vec(12) = -mu / z_knot(12);
  grad_vec(13) = -mu / z_knot(13);
  grad_vec(14) = -mu / z_knot(14);
  grad_vec(15) = -mu / z_knot(15);
end

function grad_vec = get_middle_knot_grad(z_knot)
  grad_vec = zeros(length(z_knot), 1);
  grad_vec(11) = -mu / z_knot(11);
  grad_vec(12) = -mu / z_knot(12);
  grad_vec(13) = -mu / z_knot(13);
  grad_vec(14) = -mu / z_knot(14);
  grad_vec(15) = -mu / z_knot(15);
end

function grad_vec = get_end_knot_grad(z_knot, xd)
  grad_vec = zeros(length(z_knot), 1);
  grad_vec(1) = 2.0 * z_knot(1) - 2.0 * xd(1);
  grad_vec(2) = 2.0 * z_knot(2) - 2.0 * xd(2);
  grad_vec(3) = 2.0 * z_knot(3) - 2.0 * xd(3);
  grad_vec(11) = -mu / z_knot(11);
  grad_vec(12) = -mu / z_knot(12);
  grad_vec(13) = -mu / z_knot(13);
  grad_vec(14) = -mu / z_knot(14);
  grad_vec(15) = -mu / z_knot(15);
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
  h = -z_knot(3);
  s = z_knot(4);

  knot_params.g = gravity0*((Re/(Re + h))^2);
  knot_params.g_part_h = -2*Re*gravity0/((Re + h)^3);
  knot_params.g_part_h2 = (6*Re*gravity0)/((Re + h)^4);

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
  s = z_knot(4);
  q0 = z_knot(5);
  q1 = z_knot(6);
  q2 = z_knot(7);
  q3 = z_knot(8);
  Cly = z_knot(9);
  Clz = z_knot(10);

  q = [q0 q1 q2 q3];
  D = 0.5*knot_params.rho*A*knot_params.c_D*s*s;
  az = 0.5*knot_params.rho*A*Clz*s*s/knot_params.m;
  ay = 0.5*knot_params.rho*A*Cly*s*s/knot_params.m;

  wy = (-az - knot_params.g*(-2*q1^2 - 2*q2^2 + 1))/s;
  wz = (ay + 2*knot_params.g*(q0*q1 + q2*q3))/s;
  v_NED = [s*(-2*q2^2 - 2*q3^2 + 1); 2*s*(q0*q3 + q1*q2); 2*s*(-q0*q2 + q1*q3)];
  omega = [0 0 -wy -wz; 0 0 -wz wy; wy wz 0 0; wz -wy 0 0];

  dx(1:3) = v_NED;
  dx(4) = ((knot_params.T - D)/knot_params.m) + 2*knot_params.g*(-q0*q2 + q1*q3);
  dx(5:8) = 0.5*omega*q';
end


function eq_vec = get_first_knot_eq(z_knot, t_knot, knot_params,...
                                    z_next, t_next, next_params, lb, ub, ic)

  h = t_next - t_knot;
  dx1 = get_dx(z_knot(1:10), knot_params);
  dx2 = get_dx(z_next(1:10), next_params);

  eq_vec = [...
            % the trapezoidal constraints for the knot
            z_next(1) - z_knot(1) - 0.5*h*(dx1(1) + dx2(1));...
            z_next(2) - z_knot(2) - 0.5*h*(dx1(2) + dx2(2));...
            z_next(3) - z_knot(3) - 0.5*h*(dx1(3) + dx2(3));...
            z_next(4) - z_knot(4) - 0.5*h*(dx1(4) + dx2(4));...
            z_next(5) - z_knot(5) - 0.5*h*(dx1(5) + dx2(5));...
            z_next(6) - z_knot(6) - 0.5*h*(dx1(6) + dx2(6));...
            z_next(7) - z_knot(7) - 0.5*h*(dx1(7) + dx2(7));...
            z_next(8) - z_knot(8) - 0.5*h*(dx1(8) + dx2(8));...
            % the unit quaternion constraint
            z_knot(5)^2 + z_knot(6)^2 + z_knot(7)^2 + z_knot(8)^2 - 1;...
            % the slack variable constraints for the knot
            % first lower bound
            lb(1) - z_knot(9) + z_knot(11);...
            lb(2) - z_knot(10) + z_knot(12);...
            % then upper bound
            z_knot(3) - ub(1) + z_knot(13);...
            z_knot(9) - ub(2) + z_knot(14);...
            z_knot(10) - ub(3) + z_knot(15);...
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
  dx1 = get_dx(z_knot(1:10), knot_params);
  dx2 = get_dx(z_next(1:10), next_params);

  eq_vec = [...
            % the trapezoidal constraints for the knot
            z_next(1) - z_knot(1) - 0.5*h*(dx1(1) + dx2(1));...
            z_next(2) - z_knot(2) - 0.5*h*(dx1(2) + dx2(2));...
            z_next(3) - z_knot(3) - 0.5*h*(dx1(3) + dx2(3));...
            z_next(4) - z_knot(4) - 0.5*h*(dx1(4) + dx2(4));...
            z_next(5) - z_knot(5) - 0.5*h*(dx1(5) + dx2(5));...
            z_next(6) - z_knot(6) - 0.5*h*(dx1(6) + dx2(6));...
            z_next(7) - z_knot(7) - 0.5*h*(dx1(7) + dx2(7));...
            z_next(8) - z_knot(8) - 0.5*h*(dx1(8) + dx2(8));...
            % the unit quaternion constraint
            z_knot(5)^2 + z_knot(6)^2 + z_knot(7)^2 + z_knot(8)^2 - 1;...
            % the slack variable constraints for the knot
            % first lower bound
            lb(1) - z_knot(9) + z_knot(11);...
            lb(2) - z_knot(10) + z_knot(12);...
            % then upper bound
            z_knot(3) - ub(1) + z_knot(13);...
            z_knot(9) - ub(2) + z_knot(14);...
            z_knot(10) - ub(3) + z_knot(15)];
end

function eq_vec = get_end_knot_eq(z_knot, t_knot, lb, ub)

  eq_vec = [...
            % the unit quaternion constraint
            z_knot(5)^2 + z_knot(6)^2 + z_knot(7)^2 + z_knot(8)^2 - 1;...
            % the slack variable constraints for the knot
            % first lower bound
            lb(1) - z_knot(9) + z_knot(11);...
            lb(2) - z_knot(10) + z_knot(12);...
            % then upper bound
            z_knot(3) - ub(1) + z_knot(13);...
            z_knot(9) - ub(2) + z_knot(14);...
            z_knot(10) - ub(3) + z_knot(15)];
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
  jac_block = zeros(22, 25);
  h = t_next - t_knot;
  pn1 = z_knot(1);
  pe1 = z_knot(2);
  pd1 = z_knot(3);
  s1 = z_knot(4);
  q01 = z_knot(5);
  q11 = z_knot(6);
  q21 = z_knot(7);
  q31 = z_knot(8);
  Cly1 = z_knot(9);
  Clz1 = z_knot(10);
  g1 = knot_params.g;
  m1 = knot_params.m;
  c_D1 = knot_params.c_D;
  rho1 = knot_params.rho;
  g_part_h1 = knot_params.g_part_h;
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
  Cly2 = z_next(9);
  Clz2 = z_next(10);
  g2 = next_params.g;
  m2 = next_params.m;
  c_D2 = next_params.c_D;
  rho2 = next_params.rho;
  g_part_h2 = next_params.g_part_h;
  rho_part_h2 = next_params.rho_part_h;
  c_D_part_s2 = next_params.c_D_part_s;

  jac_block(1,1) = -1;
  jac_block(1,4) = h*(q21^2 + q31^2 - 1/2);
  jac_block(1,7) = 2*h*q21*s1;
  jac_block(1,8) = 2*h*q31*s1;
  jac_block(1,16) = 1;
  jac_block(1,19) = h*(q22^2 + q32^2 - 1/2);
  jac_block(1,22) = 2*h*q22*s2;
  jac_block(1,23) = 2*h*q32*s2;
  jac_block(2,2) = -1;
  jac_block(2,4) = -h*(q01*q31 + q11*q21);
  jac_block(2,5) = -h*q31*s1;
  jac_block(2,6) = -h*q21*s1;
  jac_block(2,7) = -h*q11*s1;
  jac_block(2,8) = -h*q01*s1;
  jac_block(2,17) = 1;
  jac_block(2,19) = -h*(q02*q32 + q12*q22);
  jac_block(2,20) = -h*q32*s2;
  jac_block(2,21) = -h*q22*s2;
  jac_block(2,22) = -h*q12*s2;
  jac_block(2,23) = -h*q02*s2;
  jac_block(3,3) = -1;
  jac_block(3,4) = h*(q01*q21 - q11*q31);
  jac_block(3,5) = h*q21*s1;
  jac_block(3,6) = -h*q31*s1;
  jac_block(3,7) = h*q01*s1;
  jac_block(3,8) = -h*q11*s1;
  jac_block(3,18) = 1;
  jac_block(3,19) = h*(q02*q22 - q12*q32);
  jac_block(3,20) = h*q22*s2;
  jac_block(3,21) = -h*q32*s2;
  jac_block(3,22) = h*q02*s2;
  jac_block(3,23) = -h*q12*s2;
  jac_block(4,3) = h*(A*s1^2*c_D1*rho_part_h1 + 4*m1*(q01*q21 - q11*q31)*g_part_h1)/(4*m1);
  jac_block(4,4) = (A*h*s1*(s1*c_D_part_s1 + 2*c_D1)*rho1/4 - m1)/m1;
  jac_block(4,5) = h*q21*g1;
  jac_block(4,6) = -h*q31*g1;
  jac_block(4,7) = h*q01*g1;
  jac_block(4,8) = -h*q11*g1;
  jac_block(4,18) = h*(A*s2^2*c_D2*rho_part_h2 + 4*m2*(q02*q22 - q12*q32)*g_part_h2)/(4*m2);
  jac_block(4,19) = (A*h*s2*(s2*c_D_part_s2 + 2*c_D2)*rho2/4 + m2)/m2;
  jac_block(4,20) = h*q22*g2;
  jac_block(4,21) = -h*q32*g2;
  jac_block(4,22) = h*q02*g2;
  jac_block(4,23) = -h*q12*g2;
  jac_block(5,3) = -h*(q21*(A*Clz1*s1^2*rho_part_h1 - 2*m1*(2*q11^2 + 2*q21^2 - 1)*g_part_h1) - q31*(A*Cly1*s1^2*rho_part_h1 + 4*m1*(q01*q11 + q21*q31)*g_part_h1))/(8*m1*s1);
  jac_block(5,4) = h*(2*A*s1^2*(Cly1*q31 - Clz1*q21)*rho1 + q21*(A*Clz1*s1^2*rho1 - 2*m1*(2*q11^2 + 2*q21^2 - 1)*g1) - q31*(A*Cly1*s1^2*rho1 + 4*m1*(q01*q11 + q21*q31)*g1))/(8*m1*s1^2);
  jac_block(5,5) = (h*q11*q31*g1/2 - s1)/s1;
  jac_block(5,6) = h*(q01*q31 + 2*q11*q21)*g1/(2*s1);
  jac_block(5,7) = h*(-A*Clz1*s1^2*rho1 + 4*m1*q11^2*g1 + 12*m1*q21^2*g1 + 4*m1*q31^2*g1 - 2*m1*g1)/(8*m1*s1);
  jac_block(5,8) = A*Cly1*h*s1*rho1/(8*m1) + h*q01*q11*g1/(2*s1) + h*q21*q31*g1/s1;
  jac_block(5,9) = A*h*q31*s1*rho1/(8*m1);
  jac_block(5,10) = -A*h*q21*s1*rho1/(8*m1);
  jac_block(5,18) = -h*(q22*(A*Clz2*s2^2*rho_part_h2 - 2*m2*(2*q12^2 + 2*q22^2 -1)*g_part_h2) - q32*(A*Cly2*s2^2*rho_part_h2 + 4*m2*(q02*q12 + q22*q32)*g_part_h2))/(8*m2*s2);
  jac_block(5,19) = h*(2*A*s2^2*(Cly2*q32 - Clz2*q22)*rho2 + q22*(A*Clz2*s2^2*rho2 - 2*m2*(2*q12^2 + 2*q22^2 - 1)*g2) - q32*(A*Cly2*s2^2*rho2 + 4*m2*(q02*q12 + q22*q32)*g2))/(8*m2*s2^2);
  jac_block(5,20) = (h*q12*q32*g2/2 + s2)/s2;
  jac_block(5,21) = h*(q02*q32 + 2*q12*q22)*g2/(2*s2);
  jac_block(5,22) = h*(-A*Clz2*s2^2*rho2 + 4*m2*q12^2*g2 + 12*m2*q22^2*g2 + 4*m2*q32^2*g2 - 2*m2*g2)/(8*m2*s2);
  jac_block(5,23) = A*Cly2*h*s2*rho2/(8*m2) + h*q02*q12*g2/(2*s2) + h*q22*q32*g2/s2;
  jac_block(5,24) = A*h*q32*s2*rho2/(8*m2);
  jac_block(5,25) = -A*h*q22*s2*rho2/(8*m2);
  jac_block(6,3) = h*(A*s1^2*(Cly1*q21 + Clz1*q31)*rho_part_h1 + 2*m1*(2*q01*q11*q21 - 2*q11^2*q31 + q31)*g_part_h1)/(8*m1*s1);
  jac_block(6,4) = h*(A*s1^2*(Cly1*q21 + Clz1*q31)*rho1 + 2*m1*(-2*q01*q11*q21 + 2*q11^2*q31 - q31)*g1)/(8*m1*s1^2);
  jac_block(6,5) = h*q11*q21*g1/(2*s1);
  jac_block(6,6) = (h*(q01*q21 - 2*q11*q31)*g1/2 - s1)/s1;
  jac_block(6,7) = A*Cly1*h*s1*rho1/(8*m1) + h*q01*q11*g1/(2*s1);
  jac_block(6,8) = h*(A*Clz1*s1^2*rho1 + 2*m1*(1 - 2*q11^2)*g1)/(8*m1*s1);
  jac_block(6,9) = A*h*q21*s1*rho1/(8*m1);
  jac_block(6,10) = A*h*q31*s1*rho1/(8*m1);
  jac_block(6,18) = h*(A*s2^2*(Cly2*q22 + Clz2*q32)*rho_part_h2 + 2*m2*(2*q02*q12*q22 - 2*q12^2*q32 + q32)*g_part_h2)/(8*m2*s2);
  jac_block(6,19) = h*(A*s2^2*(Cly2*q22 + Clz2*q32)*rho2 + 2*m2*(-2*q02*q12*q22 + 2*q12^2*q32 - q32)*g2)/(8*m2*s2^2);
  jac_block(6,20) = h*q12*q22*g2/(2*s2);
  jac_block(6,21) = (h*(q02*q22 - 2*q12*q32)*g2/2 + s2)/s2;
  jac_block(6,22) = A*Cly2*h*s2*rho2/(8*m2) + h*q02*q12*g2/(2*s2);
  jac_block(6,23) = h*(A*Clz2*s2^2*rho2 + 2*m2*(1 - 2*q12^2)*g2)/(8*m2*s2);
  jac_block(6,24) = A*h*q22*s2*rho2/(8*m2);
  jac_block(6,25) = A*h*q32*s2*rho2/(8*m2);
  jac_block(7,3) = h*(q01*(A*Clz1*s1^2*rho_part_h1 - 2*m1*(2*q11^2 + 2*q21^2 - 1)*g_part_h1) - q11*(A*Cly1*s1^2*rho_part_h1 + 4*m1*(q01*q11 + q21*q31)*g_part_h1))/(8*m1*s1);
  jac_block(7,4) = h*(-A*Cly1*q11*s1^2*rho1 + A*Clz1*q01*s1^2*rho1 + 8*m1*q01*q11^2*g1 + 4*m1*q01*q21^2*g1 - 2*m1*q01*g1 + 4*m1*q11*q21*q31*g1)/(8*m1*s1^2);
  jac_block(7,5) = h*(A*Clz1*s1^2*rho1 - 8*m1*q11^2*g1 - 4*m1*q21^2*g1 + 2*m1*g1)/(8*m1*s1);
  jac_block(7,6) = h*(-A*Cly1*s1^2*rho1 - 16*m1*q01*q11*g1 - 4*m1*q21*q31*g1)/(8*m1*s1);
  jac_block(7,7) = (-h*(2*q01*q21 + q11*q31)*g1/2 - s1)/s1;
  jac_block(7,8) = -h*q11*q21*g1/(2*s1);
  jac_block(7,9) = -A*h*q11*s1*rho1/(8*m1);
  jac_block(7,10) = A*h*q01*s1*rho1/(8*m1);
  jac_block(7,18) = h*(q02*(A*Clz2*s2^2*rho_part_h2 - 2*m2*(2*q12^2 + 2*q22^2 - 1)*g_part_h2) - q12*(A*Cly2*s2^2*rho_part_h2 + 4*m2*(q02*q12 + q22*q32)*g_part_h2))/(8*m2*s2);
  jac_block(7,19) = h*(-A*Cly2*q12*s2^2*rho2 + A*Clz2*q02*s2^2*rho2 + 8*m2*q02*q12^2*g2 + 4*m2*q02*q22^2*g2 - 2*m2*q02*g2 + 4*m2*q12*q22*q32*g2)/(8*m2*s2^2);
  jac_block(7,20) = h*(A*Clz2*s2^2*rho2 - 8*m2*q12^2*g2 - 4*m2*q22^2*g2 + 2*m2*g2)/(8*m2*s2);
  jac_block(7,21) = h*(-A*Cly2*s2^2*rho2 - 16*m2*q02*q12*g2 - 4*m2*q22*q32*g2)/(8*m2*s2);
  jac_block(7,22) = (-h*(2*q02*q22 + q12*q32)*g2/2 + s2)/s2;
  jac_block(7,23) = -h*q12*q22*g2/(2*s2);
  jac_block(7,24) = -A*h*q12*s2*rho2/(8*m2);
  jac_block(7,25) = A*h*q02*s2*rho2/(8*m2);
  jac_block(8,3) = -h*(q01*(A*Cly1*s1^2*rho_part_h1 + 4*m1*(q01*q11 + q21*q31)*g_part_h1) + q11*(A*Clz1*s1^2*rho_part_h1 - 2*m1*(2*q11^2 + 2*q21^2 - 1)*g_part_h1))/(8*m1*s1);
  jac_block(8,4) = h*(-2*A*s1^2*(Cly1*q01 + Clz1*q11)*rho1 + q01*(A*Cly1*s1^2*rho1 + 4*m1*(q01*q11 + q21*q31)*g1) + q11*(A*Clz1*s1^2*rho1 - 2*m1*(2*q11^2 + 2*q21^2 - 1)*g1))/(8*m1*s1^2);
  jac_block(8,5) = -A*Cly1*h*s1*rho1/(8*m1) - h*q01*q11*g1/s1 - h*q21*q31*g1/(2*s1);
  jac_block(8,6) = h*(-A*Clz1*s1^2*rho1 - 4*m1*q01^2*g1 + 12*m1*q11^2*g1 + 4*m1*q21^2*g1 - 2*m1*g1)/(8*m1*s1);
  jac_block(8,7) = h*(-q01*q31 + 2*q11*q21)*g1/(2*s1);
  jac_block(8,8) = (-h*q01*q21*g1/2 - s1)/s1;
  jac_block(8,9) = -A*h*q01*s1*rho1/(8*m1);
  jac_block(8,10) = -A*h*q11*s1*rho1/(8*m1);
  jac_block(8,18) = -h*(q02*(A*Cly2*s2^2*rho_part_h2 + 4*m2*(q02*q12 + q22*q32)*g_part_h2) + q12*(A*Clz2*s2^2*rho_part_h2 - 2*m2*(2*q12^2 + 2*q22^2- 1)*g_part_h2))/(8*m2*s2);
  jac_block(8,19) = h*(-2*A*s2^2*(Cly2*q02 + Clz2*q12)*rho2 + q02*(A*Cly2*s2^2*rho2 + 4*m2*(q02*q12 + q22*q32)*g2) + q12*(A*Clz2*s2^2*rho2 - 2*m2*(2*q12^2 + 2*q22^2 - 1)*g2))/(8*m2*s2^2);
  jac_block(8,20) = -A*Cly2*h*s2*rho2/(8*m2) - h*q02*q12*g2/s2 - h*q22*q32*g2/(2*s2);
  jac_block(8,21) = h*(-A*Clz2*s2^2*rho2 - 4*m2*q02^2*g2 + 12*m2*q12^2*g2 + 4*m2*q22^2*g2 - 2*m2*g2)/(8*m2*s2);
  jac_block(8,22) = h*(-q02*q32 + 2*q12*q22)*g2/(2*s2);
  jac_block(8,23) = (-h*q02*q22*g2/2 + s2)/s2;
  jac_block(8,24) = -A*h*q02*s2*rho2/(8*m2);
  jac_block(8,25) = -A*h*q12*s2*rho2/(8*m2);
  jac_block(9,5) = 2*q01;
  jac_block(9,6) = 2*q11;
  jac_block(9,7) = 2*q21;
  jac_block(9,8) = 2*q31;
  jac_block(10,9) = -1;
  jac_block(10,11) = 1;
  jac_block(11,10) = -1;
  jac_block(11,12) = 1;
  jac_block(12,3) = 1;
  jac_block(12,13) = 1;
  jac_block(13,9) = 1;
  jac_block(13,14) = 1;
  jac_block(14,10) = 1;
  jac_block(14,15) = 1;
  jac_block(15,1) = -1;
  jac_block(16,2) = -1;
  jac_block(17,3) = -1;
  jac_block(18,4) = -1;
  jac_block(19,5) = -1;
  jac_block(20,6) = -1;
  jac_block(21,7) = -1;
  jac_block(22,8) = -1;


end

function jac_block = get_middle_knot_jac(z_knot, t_knot, knot_params, z_next, t_next, next_params)
  jac_block = zeros(14, 25);
  h = t_next - t_knot;
  pn1 = z_knot(1);
  pe1 = z_knot(2);
  pd1 = z_knot(3);
  s1 = z_knot(4);
  q01 = z_knot(5);
  q11 = z_knot(6);
  q21 = z_knot(7);
  q31 = z_knot(8);
  Cly1 = z_knot(9);
  Clz1 = z_knot(10);
  g1 = knot_params.g;
  m1 = knot_params.m;
  c_D1 = knot_params.c_D;
  rho1 = knot_params.rho;
  g_part_h1 = knot_params.g_part_h;
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
  Cly2 = z_next(9);
  Clz2 = z_next(10);
  g2 = next_params.g;
  m2 = next_params.m;
  c_D2 = next_params.c_D;
  rho2 = next_params.rho;
  g_part_h2 = next_params.g_part_h;
  rho_part_h2 = next_params.rho_part_h;
  c_D_part_s2 = next_params.c_D_part_s;

  jac_block(1,1) = -1;
  jac_block(1,4) = h*(q21^2 + q31^2 - 1/2);
  jac_block(1,7) = 2*h*q21*s1;
  jac_block(1,8) = 2*h*q31*s1;
  jac_block(1,16) = 1;
  jac_block(1,19) = h*(q22^2 + q32^2 - 1/2);
  jac_block(1,22) = 2*h*q22*s2;
  jac_block(1,23) = 2*h*q32*s2;
  jac_block(2,2) = -1;
  jac_block(2,4) = -h*(q01*q31 + q11*q21);
  jac_block(2,5) = -h*q31*s1;
  jac_block(2,6) = -h*q21*s1;
  jac_block(2,7) = -h*q11*s1;
  jac_block(2,8) = -h*q01*s1;
  jac_block(2,17) = 1;
  jac_block(2,19) = -h*(q02*q32 + q12*q22);
  jac_block(2,20) = -h*q32*s2;
  jac_block(2,21) = -h*q22*s2;
  jac_block(2,22) = -h*q12*s2;
  jac_block(2,23) = -h*q02*s2;
  jac_block(3,3) = -1;
  jac_block(3,4) = h*(q01*q21 - q11*q31);
  jac_block(3,5) = h*q21*s1;
  jac_block(3,6) = -h*q31*s1;
  jac_block(3,7) = h*q01*s1;
  jac_block(3,8) = -h*q11*s1;
  jac_block(3,18) = 1;
  jac_block(3,19) = h*(q02*q22 - q12*q32);
  jac_block(3,20) = h*q22*s2;
  jac_block(3,21) = -h*q32*s2;
  jac_block(3,22) = h*q02*s2;
  jac_block(3,23) = -h*q12*s2;
  jac_block(4,3) = h*(A*s1^2*c_D1*rho_part_h1 + 4*m1*(q01*q21 - q11*q31)*g_part_h1)/(4*m1);
  jac_block(4,4) = (A*h*s1*(s1*c_D_part_s1 + 2*c_D1)*rho1/4 - m1)/m1;
  jac_block(4,5) = h*q21*g1;
  jac_block(4,6) = -h*q31*g1;
  jac_block(4,7) = h*q01*g1;
  jac_block(4,8) = -h*q11*g1;
  jac_block(4,18) = h*(A*s2^2*c_D2*rho_part_h2 + 4*m2*(q02*q22 - q12*q32)*g_part_h2)/(4*m2);
  jac_block(4,19) = (A*h*s2*(s2*c_D_part_s2 + 2*c_D2)*rho2/4 + m2)/m2;
  jac_block(4,20) = h*q22*g2;
  jac_block(4,21) = -h*q32*g2;
  jac_block(4,22) = h*q02*g2;
  jac_block(4,23) = -h*q12*g2;
  jac_block(5,3) = -h*(q21*(A*Clz1*s1^2*rho_part_h1 - 2*m1*(2*q11^2 + 2*q21^2 - 1)*g_part_h1) - q31*(A*Cly1*s1^2*rho_part_h1 + 4*m1*(q01*q11 + q21*q31)*g_part_h1))/(8*m1*s1);
  jac_block(5,4) = h*(2*A*s1^2*(Cly1*q31 - Clz1*q21)*rho1 + q21*(A*Clz1*s1^2*rho1 - 2*m1*(2*q11^2 + 2*q21^2 - 1)*g1) - q31*(A*Cly1*s1^2*rho1 + 4*m1*(q01*q11 + q21*q31)*g1))/(8*m1*s1^2);
  jac_block(5,5) = (h*q11*q31*g1/2 - s1)/s1;
  jac_block(5,6) = h*(q01*q31 + 2*q11*q21)*g1/(2*s1);
  jac_block(5,7) = h*(-A*Clz1*s1^2*rho1 + 4*m1*q11^2*g1 + 12*m1*q21^2*g1 + 4*m1*q31^2*g1 - 2*m1*g1)/(8*m1*s1);
  jac_block(5,8) = A*Cly1*h*s1*rho1/(8*m1) + h*q01*q11*g1/(2*s1) + h*q21*q31*g1/s1;
  jac_block(5,9) = A*h*q31*s1*rho1/(8*m1);
  jac_block(5,10) = -A*h*q21*s1*rho1/(8*m1);
  jac_block(5,18) = -h*(q22*(A*Clz2*s2^2*rho_part_h2 - 2*m2*(2*q12^2 + 2*q22^2 -1)*g_part_h2) - q32*(A*Cly2*s2^2*rho_part_h2 + 4*m2*(q02*q12 + q22*q32)*g_part_h2))/(8*m2*s2);
  jac_block(5,19) = h*(2*A*s2^2*(Cly2*q32 - Clz2*q22)*rho2 + q22*(A*Clz2*s2^2*rho2 - 2*m2*(2*q12^2 + 2*q22^2 - 1)*g2) - q32*(A*Cly2*s2^2*rho2 + 4*m2*(q02*q12 + q22*q32)*g2))/(8*m2*s2^2);
  jac_block(5,20) = (h*q12*q32*g2/2 + s2)/s2;
  jac_block(5,21) = h*(q02*q32 + 2*q12*q22)*g2/(2*s2);
  jac_block(5,22) = h*(-A*Clz2*s2^2*rho2 + 4*m2*q12^2*g2 + 12*m2*q22^2*g2 + 4*m2*q32^2*g2 - 2*m2*g2)/(8*m2*s2);
  jac_block(5,23) = A*Cly2*h*s2*rho2/(8*m2) + h*q02*q12*g2/(2*s2) + h*q22*q32*g2/s2;
  jac_block(5,24) = A*h*q32*s2*rho2/(8*m2);
  jac_block(5,25) = -A*h*q22*s2*rho2/(8*m2);
  jac_block(6,3) = h*(A*s1^2*(Cly1*q21 + Clz1*q31)*rho_part_h1 + 2*m1*(2*q01*q11*q21 - 2*q11^2*q31 + q31)*g_part_h1)/(8*m1*s1);
  jac_block(6,4) = h*(A*s1^2*(Cly1*q21 + Clz1*q31)*rho1 + 2*m1*(-2*q01*q11*q21 + 2*q11^2*q31 - q31)*g1)/(8*m1*s1^2);
  jac_block(6,5) = h*q11*q21*g1/(2*s1);
  jac_block(6,6) = (h*(q01*q21 - 2*q11*q31)*g1/2 - s1)/s1;
  jac_block(6,7) = A*Cly1*h*s1*rho1/(8*m1) + h*q01*q11*g1/(2*s1);
  jac_block(6,8) = h*(A*Clz1*s1^2*rho1 + 2*m1*(1 - 2*q11^2)*g1)/(8*m1*s1);
  jac_block(6,9) = A*h*q21*s1*rho1/(8*m1);
  jac_block(6,10) = A*h*q31*s1*rho1/(8*m1);
  jac_block(6,18) = h*(A*s2^2*(Cly2*q22 + Clz2*q32)*rho_part_h2 + 2*m2*(2*q02*q12*q22 - 2*q12^2*q32 + q32)*g_part_h2)/(8*m2*s2);
  jac_block(6,19) = h*(A*s2^2*(Cly2*q22 + Clz2*q32)*rho2 + 2*m2*(-2*q02*q12*q22 + 2*q12^2*q32 - q32)*g2)/(8*m2*s2^2);
  jac_block(6,20) = h*q12*q22*g2/(2*s2);
  jac_block(6,21) = (h*(q02*q22 - 2*q12*q32)*g2/2 + s2)/s2;
  jac_block(6,22) = A*Cly2*h*s2*rho2/(8*m2) + h*q02*q12*g2/(2*s2);
  jac_block(6,23) = h*(A*Clz2*s2^2*rho2 + 2*m2*(1 - 2*q12^2)*g2)/(8*m2*s2);
  jac_block(6,24) = A*h*q22*s2*rho2/(8*m2);
  jac_block(6,25) = A*h*q32*s2*rho2/(8*m2);
  jac_block(7,3) = h*(q01*(A*Clz1*s1^2*rho_part_h1 - 2*m1*(2*q11^2 + 2*q21^2 - 1)*g_part_h1) - q11*(A*Cly1*s1^2*rho_part_h1 + 4*m1*(q01*q11 + q21*q31)*g_part_h1))/(8*m1*s1);
  jac_block(7,4) = h*(-A*Cly1*q11*s1^2*rho1 + A*Clz1*q01*s1^2*rho1 + 8*m1*q01*q11^2*g1 + 4*m1*q01*q21^2*g1 - 2*m1*q01*g1 + 4*m1*q11*q21*q31*g1)/(8*m1*s1^2);
  jac_block(7,5) = h*(A*Clz1*s1^2*rho1 - 8*m1*q11^2*g1 - 4*m1*q21^2*g1 + 2*m1*g1)/(8*m1*s1);
  jac_block(7,6) = h*(-A*Cly1*s1^2*rho1 - 16*m1*q01*q11*g1 - 4*m1*q21*q31*g1)/(8*m1*s1);
  jac_block(7,7) = (-h*(2*q01*q21 + q11*q31)*g1/2 - s1)/s1;
  jac_block(7,8) = -h*q11*q21*g1/(2*s1);
  jac_block(7,9) = -A*h*q11*s1*rho1/(8*m1);
  jac_block(7,10) = A*h*q01*s1*rho1/(8*m1);
  jac_block(7,18) = h*(q02*(A*Clz2*s2^2*rho_part_h2 - 2*m2*(2*q12^2 + 2*q22^2 - 1)*g_part_h2) - q12*(A*Cly2*s2^2*rho_part_h2 + 4*m2*(q02*q12 + q22*q32)*g_part_h2))/(8*m2*s2);
  jac_block(7,19) = h*(-A*Cly2*q12*s2^2*rho2 + A*Clz2*q02*s2^2*rho2 + 8*m2*q02*q12^2*g2 + 4*m2*q02*q22^2*g2 - 2*m2*q02*g2 + 4*m2*q12*q22*q32*g2)/(8*m2*s2^2);
  jac_block(7,20) = h*(A*Clz2*s2^2*rho2 - 8*m2*q12^2*g2 - 4*m2*q22^2*g2 + 2*m2*g2)/(8*m2*s2);
  jac_block(7,21) = h*(-A*Cly2*s2^2*rho2 - 16*m2*q02*q12*g2 - 4*m2*q22*q32*g2)/(8*m2*s2);
  jac_block(7,22) = (-h*(2*q02*q22 + q12*q32)*g2/2 + s2)/s2;
  jac_block(7,23) = -h*q12*q22*g2/(2*s2);
  jac_block(7,24) = -A*h*q12*s2*rho2/(8*m2);
  jac_block(7,25) = A*h*q02*s2*rho2/(8*m2);
  jac_block(8,3) = -h*(q01*(A*Cly1*s1^2*rho_part_h1 + 4*m1*(q01*q11 + q21*q31)*g_part_h1) + q11*(A*Clz1*s1^2*rho_part_h1 - 2*m1*(2*q11^2 + 2*q21^2 - 1)*g_part_h1))/(8*m1*s1);
  jac_block(8,4) = h*(-2*A*s1^2*(Cly1*q01 + Clz1*q11)*rho1 + q01*(A*Cly1*s1^2*rho1 + 4*m1*(q01*q11 + q21*q31)*g1) + q11*(A*Clz1*s1^2*rho1 - 2*m1*(2*q11^2 + 2*q21^2 - 1)*g1))/(8*m1*s1^2);
  jac_block(8,5) = -A*Cly1*h*s1*rho1/(8*m1) - h*q01*q11*g1/s1 - h*q21*q31*g1/(2*s1);
  jac_block(8,6) = h*(-A*Clz1*s1^2*rho1 - 4*m1*q01^2*g1 + 12*m1*q11^2*g1 + 4*m1*q21^2*g1 - 2*m1*g1)/(8*m1*s1);
  jac_block(8,7) = h*(-q01*q31 + 2*q11*q21)*g1/(2*s1);
  jac_block(8,8) = (-h*q01*q21*g1/2 - s1)/s1;
  jac_block(8,9) = -A*h*q01*s1*rho1/(8*m1);
  jac_block(8,10) = -A*h*q11*s1*rho1/(8*m1);
  jac_block(8,18) = -h*(q02*(A*Cly2*s2^2*rho_part_h2 + 4*m2*(q02*q12 + q22*q32)*g_part_h2) + q12*(A*Clz2*s2^2*rho_part_h2 - 2*m2*(2*q12^2 + 2*q22^2- 1)*g_part_h2))/(8*m2*s2);
  jac_block(8,19) = h*(-2*A*s2^2*(Cly2*q02 + Clz2*q12)*rho2 + q02*(A*Cly2*s2^2*rho2 + 4*m2*(q02*q12 + q22*q32)*g2) + q12*(A*Clz2*s2^2*rho2 - 2*m2*(2*q12^2 + 2*q22^2 - 1)*g2))/(8*m2*s2^2);
  jac_block(8,20) = -A*Cly2*h*s2*rho2/(8*m2) - h*q02*q12*g2/s2 - h*q22*q32*g2/(2*s2);
  jac_block(8,21) = h*(-A*Clz2*s2^2*rho2 - 4*m2*q02^2*g2 + 12*m2*q12^2*g2 + 4*m2*q22^2*g2 - 2*m2*g2)/(8*m2*s2);
  jac_block(8,22) = h*(-q02*q32 + 2*q12*q22)*g2/(2*s2);
  jac_block(8,23) = (-h*q02*q22*g2/2 + s2)/s2;
  jac_block(8,24) = -A*h*q02*s2*rho2/(8*m2);
  jac_block(8,25) = -A*h*q12*s2*rho2/(8*m2);
  jac_block(9,5) = 2*q01;
  jac_block(9,6) = 2*q11;
  jac_block(9,7) = 2*q21;
  jac_block(9,8) = 2*q31;
  jac_block(10,9) = -1;
  jac_block(10,11) = 1;
  jac_block(11,10) = -1;
  jac_block(11,12) = 1;
  jac_block(12,3) = 1;
  jac_block(12,13) = 1;
  jac_block(13,9) = 1;
  jac_block(13,14) = 1;
  jac_block(14,10) = 1;
  jac_block(14,15) = 1;
end

function jac_block = get_end_knot_jac(z_knot)
  jac_block = zeros(6, 15);
  q01 = z_knot(5);
  q11 = z_knot(6);
  q21 = z_knot(7);
  q31 = z_knot(8);

  jac_block(1,5) = 2*q01;
  jac_block(1,6) = 2*q11;
  jac_block(1,7) = 2*q21;
  jac_block(1,8) = 2*q31;
  jac_block(2,9) = -1;
  jac_block(2,11) = 1;
  jac_block(3,10) = -1;
  jac_block(3,12) = 1;
  jac_block(4,3) = 1;
  jac_block(4,13) = 1;
  jac_block(5,9) = 1;
  jac_block(5,14) = 1;
  jac_block(6,10) = 1;
  jac_block(6,15) = 1;
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
  hes_block = zeros(15, 15);
  h = t_next - t_knot;
  pn1 = z_knot(1);
  pe1 = z_knot(2);
  pd1 = z_knot(3);
  s1 = z_knot(4);
  q01 = z_knot(5);
  q11 = z_knot(6);
  q21 = z_knot(7);
  q31 = z_knot(8);
  Cly1 = z_knot(9);
  Clz1 = z_knot(10);
  s_Cly_lb = z_knot(11);
  s_Clz_lb = z_knot(12);
  s_pd_ub = z_knot(13);
  s_Cly_ub = z_knot(14);
  s_Clz_ub = z_knot(15);
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

end

function hes_block = get_middle_knot_hes(z_knot, t_knot, knot_params, lambda_knot, t_next, lambda_last)
  hes_block = zeros(15, 15);
  h = t_next - t_knot;
  pn1 = z_knot(1);
  pe1 = z_knot(2);
  pd1 = z_knot(3);
  s1 = z_knot(4);
  q01 = z_knot(5);
  q11 = z_knot(6);
  q21 = z_knot(7);
  q31 = z_knot(8);
  Cly1 = z_knot(9);
  Clz1 = z_knot(10);
  s_Cly_lb = z_knot(11);
  s_Clz_lb = z_knot(12);
  s_pd_ub = z_knot(13);
  s_Cly_ub = z_knot(14);
  s_Clz_ub = z_knot(15);
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

end

function hes_block = get_end_knot_hes(z_knot, t_knot, knot_params, lambda_knot, t_next, lambda_last)
  hes_block = zeros(15, 15);
  h = t_next - t_knot;
  pn1 = z_knot(1);
  pe1 = z_knot(2);
  pd1 = z_knot(3);
  s1 = z_knot(4);
  q01 = z_knot(5);
  q11 = z_knot(6);
  q21 = z_knot(7);
  q31 = z_knot(8);
  Cly1 = z_knot(9);
  Clz1 = z_knot(10);
  s_Cly_lb = z_knot(11);
  s_Clz_lb = z_knot(12);
  s_pd_ub = z_knot(13);
  s_Cly_ub = z_knot(14);
  s_Clz_ub = z_knot(15);
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
