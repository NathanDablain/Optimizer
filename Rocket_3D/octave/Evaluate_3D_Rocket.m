function [cost, grad, eq, jacobian, hessian] = Evaluate_3D_Rocket(O, option)
% Initialize outputs
cost = [];
grad = [];
eq = [];
jacobian = [];
hessian = [];
% knot vars are :
% 1 - pn,
% 2 - pe,
% 3 - pd,
% 4 - vn,
% 5 - ve,
% 6 - vd,
% 7 - an,
% 8 - ae,
% 9 - ad,
% 10 - c_L
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                       Evaluate cost                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  mu = 1.0e-3;
  function contribution = get_first_knot_cost(z_knot, lb, ub)

    contribution = -mu*(log(z_knot(10) - lb(1)) +...
                        log(ub(1) - z_knot(3)) +...
                        log(ub(2) - z_knot(10)) +...
                        log(z_knot(7) + 1.1) +...
                        log(z_knot(8) + 1.1) +...
                        log(z_knot(9) + 1.1) +...
                        log(1.1 - z_knot(7)) +...
                        log(1.1 - z_knot(8)) +...
                        log(1.1 - z_knot(9))...
                        );
  end

  function contribution = get_middle_knot_cost(z_knot, lb, ub)

    contribution = -mu*(log(z_knot(10) - lb(1)) +...
                        log(ub(1) - z_knot(3)) +...
                        log(ub(2) - z_knot(10)) +...
                        log(z_knot(7) + 1.1) +...
                        log(z_knot(8) + 1.1) +...
                        log(z_knot(9) + 1.1) +...
                        log(1.1 - z_knot(7)) +...
                        log(1.1 - z_knot(8)) +...
                        log(1.1 - z_knot(9))...
                        );
  end

  function contribution = get_end_knot_cost(z_knot, xd, lb, ub)
    contribution = (z_knot(1) - xd(1))^2 +...
                   (z_knot(2) - xd(2))^2 +...
                   (z_knot(3) - xd(3))^2 -...
                    mu*(log(z_knot(10) - lb(1)) +...
                        log(ub(1) - z_knot(3)) +...
                        log(ub(2) - z_knot(10)) +...
                        log(z_knot(7) + 1.1) +...
                        log(z_knot(8) + 1.1) +...
                        log(z_knot(9) + 1.1) +...
                        log(1.1 - z_knot(7)) +...
                        log(1.1 - z_knot(8)) +...
                        log(1.1 - z_knot(9))...
                        );;

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
  grad_vec(3) = -mu/(z_knot(3) - ub(1));
  grad_vec(7) = -mu/(z_knot(7) - 1.1) + mu/(-1.1 - z_knot(7));
  grad_vec(8) = -mu/(z_knot(8) - 1.1) + mu/(-1.1 - z_knot(8));
  grad_vec(9) = -mu/(z_knot(9) - 1.1) + mu/(-1.1 - z_knot(9));
  grad_vec(10) = -mu/(z_knot(10) - ub(2)) + mu/(lb(1) - z_knot(10));
end

function grad_vec = get_middle_knot_grad(z_knot, lb, ub)
  grad_vec = zeros(length(z_knot), 1);
  grad_vec(3) = -mu/(z_knot(3) - ub(1));
  grad_vec(7) = -mu/(z_knot(7) - 1.1) + mu/(-1.1 - z_knot(7));
  grad_vec(8) = -mu/(z_knot(8) - 1.1) + mu/(-1.1 - z_knot(8));
  grad_vec(9) = -mu/(z_knot(9) - 1.1) + mu/(-1.1 - z_knot(9));
  grad_vec(10) = -mu/(z_knot(10) - ub(2)) + mu/(lb(1) - z_knot(10));
end

function grad_vec = get_end_knot_grad(z_knot, xd, lb, ub)
  grad_vec = zeros(length(z_knot), 1);
  grad_vec(1) = 2.0 * z_knot(1) - 2.0 * xd(1);
  grad_vec(2) = 2.0 * z_knot(2) - 2.0 * xd(2);
  grad_vec(3) = 2.0 * z_knot(3) - 2.0 * xd(3) - mu/(z_knot(3) - ub(1));
  grad_vec(7) = -mu/(z_knot(7) - 1.1) + mu/(-1.1 - z_knot(7));
  grad_vec(8) = -mu/(z_knot(8) - 1.1) + mu/(-1.1 - z_knot(8));
  grad_vec(9) = -mu/(z_knot(9) - 1.1) + mu/(-1.1 - z_knot(9));
  grad_vec(10) = -mu/(z_knot(10) - ub(2)) + mu/(lb(1) - z_knot(10));
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
  s = norm(z_knot(4:6));

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
  vn = z_knot(4);
  ve = z_knot(5);
  vd = z_knot(6);
  an = z_knot(7);
  ae = z_knot(8);
  ad = z_knot(9);
  c_L = z_knot(10);

  v_vec = [vn; ve; vd];
  v_uv = v_vec ./ norm(v_vec);
  a_uv = [an; ae; ad];
  g_vec = [0;0;knot_params.g];

  D = 0.5*knot_params.rho*A*knot_params.c_D*norm(v_vec)^2;
  L = 0.5*knot_params.rho*A*c_L*norm(v_vec)^2;
  dx(1:3) = v_vec;
  dx(4:6) = ((knot_params.T - D)*v_uv + L*a_uv)/knot_params.m + g_vec;
end

function eq_vec = get_first_knot_eq(z_knot, t_knot, knot_params,...
                                    z_next, t_next, next_params, lb, ub, ic)

  h = t_next - t_knot;
  dx1 = get_dx(z_knot, knot_params);
  dx2 = get_dx(z_next, next_params);

  eq_vec = [...
            % the trapezoidal constraints for the knot
            z_next(1) - z_knot(1) - 0.5*h*(dx1(1) + dx2(1));...
            z_next(2) - z_knot(2) - 0.5*h*(dx1(2) + dx2(2));...
            z_next(3) - z_knot(3) - 0.5*h*(dx1(3) + dx2(3));...
            z_next(4) - z_knot(4) - 0.5*h*(dx1(4) + dx2(4));...
            z_next(5) - z_knot(5) - 0.5*h*(dx1(5) + dx2(5));...
            z_next(6) - z_knot(6) - 0.5*h*(dx1(6) + dx2(6));...
             % the orthogonality constraint
            dot(z_knot(4:6),z_knot(7:9));...
            % the unit length constraint
            norm(z_knot(7:9)) - 1.0;...
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
  dx1 = get_dx(z_knot, knot_params);
  dx2 = get_dx(z_next, next_params);

  eq_vec = [...
            % the trapezoidal constraints for the knot
            z_next(1) - z_knot(1) - 0.5*h*(dx1(1) + dx2(1));...
            z_next(2) - z_knot(2) - 0.5*h*(dx1(2) + dx2(2));...
            z_next(3) - z_knot(3) - 0.5*h*(dx1(3) + dx2(3));...
            z_next(4) - z_knot(4) - 0.5*h*(dx1(4) + dx2(4));...
            z_next(5) - z_knot(5) - 0.5*h*(dx1(5) + dx2(5));...
            z_next(6) - z_knot(6) - 0.5*h*(dx1(6) + dx2(6));...
             % the orthogonality constraint
            dot(z_knot(4:6),z_knot(7:9));...
            % the unit length constraint
            norm(z_knot(7:9)) - 1.0];

end

function eq_vec = get_end_knot_eq(z_knot, t_knot, knot_params, lb, ub)

  eq_vec = [...
            % the orthogonality constraint
            dot(z_knot(4:6),z_knot(7:9));...
            % the unit length constraint
            norm(z_knot(7:9)) - 1.0];
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
  jac_block = zeros(14, 20);
  h = t_next - t_knot;
  pn1 = z_knot(1);
  pe1 = z_knot(2);
  pd1 = z_knot(3);
  vn1 = z_knot(4);
  ve1 = z_knot(5);
  vd1 = z_knot(6);
  an1 = z_knot(7);
  ae1 = z_knot(8);
  ad1 = z_knot(9);
  c_L1 = z_knot(10);
  g1 = knot_params.g;
  T1 = knot_params.T;
  m1 = knot_params.m;
  c_D1 = knot_params.c_D;
  rho1 = knot_params.rho;
  s1 = norm([vn1 ve1 vd1]);

  pn2 = z_next(1);
  pe2 = z_next(2);
  pd2 = z_next(3);
  vn2 = z_next(4);
  ve2 = z_next(5);
  vd2 = z_next(6);
  an2 = z_next(7);
  ae2 = z_next(8);
  ad2 = z_next(9);
  c_L2 = z_next(10);
  g2 = next_params.g;
  T2 = next_params.T;
  m2 = next_params.m;
  c_D2 = next_params.c_D;
  rho2 = next_params.rho;
  s2 = norm([vn2 ve2 vd2]);

  jac_block(1,1) = -1;
  jac_block(1,4) = -h/2;
  jac_block(1,11) = 1;
  jac_block(1,14) = -h/2;
  jac_block(2,2) = -1;
  jac_block(2,5) = -h/2;
  jac_block(2,12) = 1;
  jac_block(2,15) = -h/2;
  jac_block(3,3) = -1;
  jac_block(3,6) = -h/2;
  jac_block(3,13) = 1;
  jac_block(3,16) = -h/2;
  jac_block(4,4) = (-2*A*an1*c_L1*h*rho1*vd1^2*vn1*s1 - 2*A*an1*c_L1*h*rho1*ve1^2*vn1*s1 - 2*A*an1*c_L1*h*rho1*vn1^3*s1 + A*c_D1*h*rho1*vd1^4 + 2*A*c_D1*h*rho1*vd1^2*ve1^2 + 3*A*c_D1*h*rho1*vd1^2*vn1^2 + A*c_D1*h*rho1*ve1^4 + 3*A*c_D1*h*rho1*ve1^2*vn1^2 + 2*A*c_D1*h*rho1*vn1^4 - 2*T1*h*vd1^2 - 2*T1*h*ve1^2 - 4*m1*vd1^2*s1 - 4*m1*ve1^2*s1 - 4*m1*vn1^2*s1)/(4*m1*s1^2^(3/2));
  jac_block(4,5) = h*ve1*(-2*A*an1*c_L1*rho1*vd1^2*s1 - 2*A*an1*c_L1*rho1*ve1^2*s1 - 2*A*an1*c_L1*rho1*vn1^2*s1 + A*c_D1*rho1*vd1^2*vn1 + A*c_D1*rho1*ve1^2*vn1 + A*c_D1*rho1*vn1^3 + 2*T1*vn1)/(4*m1*s1^2^(3/2));
  jac_block(4,6) = h*vd1*(-2*A*an1*c_L1*rho1*vd1^2*s1 - 2*A*an1*c_L1*rho1*ve1^2*s1 - 2*A*an1*c_L1*rho1*vn1^2*s1 + A*c_D1*rho1*vd1^2*vn1 + A*c_D1*rho1*ve1^2*vn1 + A*c_D1*rho1*vn1^3 + 2*T1*vn1)/(4*m1*s1^2^(3/2));
  jac_block(4,7) = -A*c_L1*h*rho1*s1^2/(4*m1);
  jac_block(4,10) = -A*an1*h*rho1*s1^2/(4*m1);
  jac_block(4,14) = (-2*A*an2*c_L2*h*rho2*vd2^2*vn2*s2 - 2*A*an2*c_L2*h*rho2*ve2^2*vn2*s2 - 2*A*an2*c_L2*h*rho2*vn2^3*s2 + A*c_D2*h*rho2*vd2^4 + 2*A*c_D2*h*rho2*vd2^2*ve2^2 + 3*A*c_D2*h*rho2*vd2^2*vn2^2 + A*c_D2*h*rho2*ve2^4 + 3*A*c_D2*h*rho2*ve2^2*vn2^2 + 2*A*c_D2*h*rho2*vn2^4 - 2*T2*h*vd2^2 - 2*T2*h*ve2^2 + 4*m2*vd2^2*s2 + 4*m2*ve2^2*s2 + 4*m2*vn2^2*s2)/(4*m2*s2^2^(3/2));
  jac_block(4,15) = h*ve2*(-2*A*an2*c_L2*rho2*vd2^2*s2 - 2*A*an2*c_L2*rho2*ve2^2*s2 - 2*A*an2*c_L2*rho2*vn2^2*s2 + A*c_D2*rho2*vd2^2*vn2 + A*c_D2*rho2*ve2^2*vn2 + A*c_D2*rho2*vn2^3 + 2*T2*vn2)/(4*m2*s2^2^(3/2));
  jac_block(4,16) = h*vd2*(-2*A*an2*c_L2*rho2*vd2^2*s2 - 2*A*an2*c_L2*rho2*ve2^2*s2 - 2*A*an2*c_L2*rho2*vn2^2*s2 + A*c_D2*rho2*vd2^2*vn2 + A*c_D2*rho2*ve2^2*vn2 + A*c_D2*rho2*vn2^3 + 2*T2*vn2)/(4*m2*s2^2^(3/2));
  jac_block(4,17) = -A*c_L2*h*rho2*s2^2/(4*m2);
  jac_block(4,20) = -A*an2*h*rho2*s2^2/(4*m2);
  jac_block(5,4) = h*vn1*(-2*A*ae1*c_L1*rho1*vd1^2*s1 - 2*A*ae1*c_L1*rho1*ve1^2*s1 - 2*A*ae1*c_L1*rho1*vn1^2*s1 + A*c_D1*rho1*vd1^2*ve1 + A*c_D1*rho1*ve1^3 + A*c_D1*rho1*ve1*vn1^2 + 2*T1*ve1)/(4*m1*s1^2^(3/2));
  jac_block(5,5) = (-2*A*ae1*c_L1*h*rho1*vd1^2*ve1*s1 - 2*A*ae1*c_L1*h*rho1*ve1^3*s1 - 2*A*ae1*c_L1*h*rho1*ve1*vn1^2*s1 + A*c_D1*h*rho1*vd1^4 + 3*A*c_D1*h*rho1*vd1^2*ve1^2 + 2*A*c_D1*h*rho1*vd1^2*vn1^2 + 2*A*c_D1*h*rho1*ve1^4 + 3*A*c_D1*h*rho1*ve1^2*vn1^2 + A*c_D1*h*rho1*vn1^4 - 2*T1*h*vd1^2 - 2*T1*h*vn1^2 - 4*m1*vd1^2*s1 - 4*m1*ve1^2*s1 - 4*m1*vn1^2*s1)/(4*m1*s1^2^(3/2));
  jac_block(5,6) = h*vd1*(-2*A*ae1*c_L1*rho1*vd1^2*s1 - 2*A*ae1*c_L1*rho1*ve1^2*s1 - 2*A*ae1*c_L1*rho1*vn1^2*s1 + A*c_D1*rho1*vd1^2*ve1 + A*c_D1*rho1*ve1^3 + A*c_D1*rho1*ve1*vn1^2 + 2*T1*ve1)/(4*m1*s1^2^(3/2));
  jac_block(5,8) = -A*c_L1*h*rho1*s1^2/(4*m1);
  jac_block(5,10) = -A*ae1*h*rho1*s1^2/(4*m1);
  jac_block(5,14) = h*vn2*(-2*A*ae2*c_L2*rho2*vd2^2*s2 - 2*A*ae2*c_L2*rho2*ve2^2*s2 - 2*A*ae2*c_L2*rho2*vn2^2*s2 + A*c_D2*rho2*vd2^2*ve2 + A*c_D2*rho2*ve2^3 + A*c_D2*rho2*ve2*vn2^2 + 2*T2*ve2)/(4*m2*s2^2^(3/2));
  jac_block(5,15) = (-2*A*ae2*c_L2*h*rho2*vd2^2*ve2*s2 - 2*A*ae2*c_L2*h*rho2*ve2^3*s2 - 2*A*ae2*c_L2*h*rho2*ve2*vn2^2*s2 + A*c_D2*h*rho2*vd2^4 + 3*A*c_D2*h*rho2*vd2^2*ve2^2 + 2*A*c_D2*h*rho2*vd2^2*vn2^2 + 2*A*c_D2*h*rho2*ve2^4 + 3*A*c_D2*h*rho2*ve2^2*vn2^2 + A*c_D2*h*rho2*vn2^4 - 2*T2*h*vd2^2 - 2*T2*h*vn2^2 + 4*m2*vd2^2*s2 + 4*m2*ve2^2*s2 + 4*m2*vn2^2*s2)/(4*m2*s2^2^(3/2));
  jac_block(5,16) = h*vd2*(-2*A*ae2*c_L2*rho2*vd2^2*s2 - 2*A*ae2*c_L2*rho2*ve2^2*s2 - 2*A*ae2*c_L2*rho2*vn2^2*s2 + A*c_D2*rho2*vd2^2*ve2 + A*c_D2*rho2*ve2^3 + A*c_D2*rho2*ve2*vn2^2 + 2*T2*ve2)/(4*m2*s2^2^(3/2));
  jac_block(5,18) = -A*c_L2*h*rho2*s2^2/(4*m2);
  jac_block(5,20) = -A*ae2*h*rho2*s2^2/(4*m2);
  jac_block(6,4) = h*vn1*(-2*A*ad1*c_L1*rho1*vd1^2*s1 - 2*A*ad1*c_L1*rho1*ve1^2*s1 - 2*A*ad1*c_L1*rho1*vn1^2*s1 + A*c_D1*rho1*vd1^3 + A*c_D1*rho1*vd1*ve1^2 + A*c_D1*rho1*vd1*vn1^2 + 2*T1*vd1)/(4*m1*s1^2^(3/2));
  jac_block(6,5) = h*ve1*(-2*A*ad1*c_L1*rho1*vd1^2*s1 - 2*A*ad1*c_L1*rho1*ve1^2*s1 - 2*A*ad1*c_L1*rho1*vn1^2*s1 + A*c_D1*rho1*vd1^3 + A*c_D1*rho1*vd1*ve1^2 + A*c_D1*rho1*vd1*vn1^2 + 2*T1*vd1)/(4*m1*s1^2^(3/2));
  jac_block(6,6) = (-2*A*ad1*c_L1*h*rho1*vd1^3*s1 - 2*A*ad1*c_L1*h*rho1*vd1*ve1^2*s1 - 2*A*ad1*c_L1*h*rho1*vd1*vn1^2*s1 + 2*A*c_D1*h*rho1*vd1^4 + 3*A*c_D1*h*rho1*vd1^2*ve1^2 + 3*A*c_D1*h*rho1*vd1^2*vn1^2 + A*c_D1*h*rho1*ve1^4 + 2*A*c_D1*h*rho1*ve1^2*vn1^2 + A*c_D1*h*rho1*vn1^4 - 2*T1*h*ve1^2 - 2*T1*h*vn1^2 - 4*m1*vd1^2*s1 - 4*m1*ve1^2*s1 - 4*m1*vn1^2*s1)/(4*m1*s1^2^(3/2));
  jac_block(6,9) = -A*c_L1*h*rho1*s1^2/(4*m1);
  jac_block(6,10) = -A*ad1*h*rho1*s1^2/(4*m1);
  jac_block(6,14) = h*vn2*(-2*A*ad2*c_L2*rho2*vd2^2*s2 - 2*A*ad2*c_L2*rho2*ve2^2*s2 - 2*A*ad2*c_L2*rho2*vn2^2*s2 + A*c_D2*rho2*vd2^3 + A*c_D2*rho2*vd2*ve2^2 + A*c_D2*rho2*vd2*vn2^2 + 2*T2*vd2)/(4*m2*s2^2^(3/2));
  jac_block(6,15) = h*ve2*(-2*A*ad2*c_L2*rho2*vd2^2*s2 - 2*A*ad2*c_L2*rho2*ve2^2*s2 - 2*A*ad2*c_L2*rho2*vn2^2*s2 + A*c_D2*rho2*vd2^3 + A*c_D2*rho2*vd2*ve2^2 + A*c_D2*rho2*vd2*vn2^2 + 2*T2*vd2)/(4*m2*s2^2^(3/2));
  jac_block(6,16) = (-2*A*ad2*c_L2*h*rho2*vd2^3*s2 - 2*A*ad2*c_L2*h*rho2*vd2*ve2^2*s2 - 2*A*ad2*c_L2*h*rho2*vd2*vn2^2*s2 + 2*A*c_D2*h*rho2*vd2^4 + 3*A*c_D2*h*rho2*vd2^2*ve2^2 + 3*A*c_D2*h*rho2*vd2^2*vn2^2 + A*c_D2*h*rho2*ve2^4 + 2*A*c_D2*h*rho2*ve2^2*vn2^2 + A*c_D2*h*rho2*vn2^4 - 2*T2*h*ve2^2 - 2*T2*h*vn2^2 + 4*m2*vd2^2*s2 + 4*m2*ve2^2*s2 + 4*m2*vn2^2*s2)/(4*m2*s2^2^(3/2));
  jac_block(6,19) = -A*c_L2*h*rho2*s2^2/(4*m2);
  jac_block(6,20) = -A*ad2*h*rho2*s2^2/(4*m2);
  jac_block(7,4) = an1;
  jac_block(7,5) = ae1;
  jac_block(7,6) = ad1;
  jac_block(7,7) = vn1;
  jac_block(7,8) = ve1;
  jac_block(7,9) = vd1;
  jac_block(8,7) = an1/sqrt(ad1^2 + ae1^2 + an1^2);
  jac_block(8,8) = ae1/sqrt(ad1^2 + ae1^2 + an1^2);
  jac_block(8,9) = ad1/sqrt(ad1^2 + ae1^2 + an1^2);
  jac_block(9,1) = -1;
  jac_block(10,2) = -1;
  jac_block(11,3) = -1;
  jac_block(12,4) = -1;
  jac_block(13,5) = -1;
  jac_block(14,6) = -1;

end

function jac_block = get_middle_knot_jac(z_knot, t_knot, knot_params, z_next, t_next, next_params)
  jac_block = zeros(8, 20);
  h = t_next - t_knot;
  pn1 = z_knot(1);
  pe1 = z_knot(2);
  pd1 = z_knot(3);
  vn1 = z_knot(4);
  ve1 = z_knot(5);
  vd1 = z_knot(6);
  an1 = z_knot(7);
  ae1 = z_knot(8);
  ad1 = z_knot(9);
  c_L1 = z_knot(10);
  g1 = knot_params.g;
  T1 = knot_params.T;
  m1 = knot_params.m;
  c_D1 = knot_params.c_D;
  rho1 = knot_params.rho;
  s1 = norm([vn1 ve1 vd1]);

  pn2 = z_next(1);
  pe2 = z_next(2);
  pd2 = z_next(3);
  vn2 = z_next(4);
  ve2 = z_next(5);
  vd2 = z_next(6);
  an2 = z_next(7);
  ae2 = z_next(8);
  ad2 = z_next(9);
  c_L2 = z_next(10);
  g2 = next_params.g;
  T2 = next_params.T;
  m2 = next_params.m;
  c_D2 = next_params.c_D;
  rho2 = next_params.rho;
  s2 = norm([vn2 ve2 vd2]);

  jac_block(1,1) = -1;
  jac_block(1,4) = -h/2;
  jac_block(1,11) = 1;
  jac_block(1,14) = -h/2;
  jac_block(2,2) = -1;
  jac_block(2,5) = -h/2;
  jac_block(2,12) = 1;
  jac_block(2,15) = -h/2;
  jac_block(3,3) = -1;
  jac_block(3,6) = -h/2;
  jac_block(3,13) = 1;
  jac_block(3,16) = -h/2;
  jac_block(4,4) = (-2*A*an1*c_L1*h*rho1*vd1^2*vn1*s1 - 2*A*an1*c_L1*h*rho1*ve1^2*vn1*s1 - 2*A*an1*c_L1*h*rho1*vn1^3*s1 + A*c_D1*h*rho1*vd1^4 + 2*A*c_D1*h*rho1*vd1^2*ve1^2 + 3*A*c_D1*h*rho1*vd1^2*vn1^2 + A*c_D1*h*rho1*ve1^4 + 3*A*c_D1*h*rho1*ve1^2*vn1^2 + 2*A*c_D1*h*rho1*vn1^4 - 2*T1*h*vd1^2 - 2*T1*h*ve1^2 - 4*m1*vd1^2*s1 - 4*m1*ve1^2*s1 - 4*m1*vn1^2*s1)/(4*m1*s1^2^(3/2));
  jac_block(4,5) = h*ve1*(-2*A*an1*c_L1*rho1*vd1^2*s1 - 2*A*an1*c_L1*rho1*ve1^2*s1 - 2*A*an1*c_L1*rho1*vn1^2*s1 + A*c_D1*rho1*vd1^2*vn1 + A*c_D1*rho1*ve1^2*vn1 + A*c_D1*rho1*vn1^3 + 2*T1*vn1)/(4*m1*s1^2^(3/2));
  jac_block(4,6) = h*vd1*(-2*A*an1*c_L1*rho1*vd1^2*s1 - 2*A*an1*c_L1*rho1*ve1^2*s1 - 2*A*an1*c_L1*rho1*vn1^2*s1 + A*c_D1*rho1*vd1^2*vn1 + A*c_D1*rho1*ve1^2*vn1 + A*c_D1*rho1*vn1^3 + 2*T1*vn1)/(4*m1*s1^2^(3/2));
  jac_block(4,7) = -A*c_L1*h*rho1*s1^2/(4*m1);
  jac_block(4,10) = -A*an1*h*rho1*s1^2/(4*m1);
  jac_block(4,14) = (-2*A*an2*c_L2*h*rho2*vd2^2*vn2*s2 - 2*A*an2*c_L2*h*rho2*ve2^2*vn2*s2 - 2*A*an2*c_L2*h*rho2*vn2^3*s2 + A*c_D2*h*rho2*vd2^4 + 2*A*c_D2*h*rho2*vd2^2*ve2^2 + 3*A*c_D2*h*rho2*vd2^2*vn2^2 + A*c_D2*h*rho2*ve2^4 + 3*A*c_D2*h*rho2*ve2^2*vn2^2 + 2*A*c_D2*h*rho2*vn2^4 - 2*T2*h*vd2^2 - 2*T2*h*ve2^2 + 4*m2*vd2^2*s2 + 4*m2*ve2^2*s2 + 4*m2*vn2^2*s2)/(4*m2*s2^2^(3/2));
  jac_block(4,15) = h*ve2*(-2*A*an2*c_L2*rho2*vd2^2*s2 - 2*A*an2*c_L2*rho2*ve2^2*s2 - 2*A*an2*c_L2*rho2*vn2^2*s2 + A*c_D2*rho2*vd2^2*vn2 + A*c_D2*rho2*ve2^2*vn2 + A*c_D2*rho2*vn2^3 + 2*T2*vn2)/(4*m2*s2^2^(3/2));
  jac_block(4,16) = h*vd2*(-2*A*an2*c_L2*rho2*vd2^2*s2 - 2*A*an2*c_L2*rho2*ve2^2*s2 - 2*A*an2*c_L2*rho2*vn2^2*s2 + A*c_D2*rho2*vd2^2*vn2 + A*c_D2*rho2*ve2^2*vn2 + A*c_D2*rho2*vn2^3 + 2*T2*vn2)/(4*m2*s2^2^(3/2));
  jac_block(4,17) = -A*c_L2*h*rho2*s2^2/(4*m2);
  jac_block(4,20) = -A*an2*h*rho2*s2^2/(4*m2);
  jac_block(5,4) = h*vn1*(-2*A*ae1*c_L1*rho1*vd1^2*s1 - 2*A*ae1*c_L1*rho1*ve1^2*s1 - 2*A*ae1*c_L1*rho1*vn1^2*s1 + A*c_D1*rho1*vd1^2*ve1 + A*c_D1*rho1*ve1^3 + A*c_D1*rho1*ve1*vn1^2 + 2*T1*ve1)/(4*m1*s1^2^(3/2));
  jac_block(5,5) = (-2*A*ae1*c_L1*h*rho1*vd1^2*ve1*s1 - 2*A*ae1*c_L1*h*rho1*ve1^3*s1 - 2*A*ae1*c_L1*h*rho1*ve1*vn1^2*s1 + A*c_D1*h*rho1*vd1^4 + 3*A*c_D1*h*rho1*vd1^2*ve1^2 + 2*A*c_D1*h*rho1*vd1^2*vn1^2 + 2*A*c_D1*h*rho1*ve1^4 + 3*A*c_D1*h*rho1*ve1^2*vn1^2 + A*c_D1*h*rho1*vn1^4 - 2*T1*h*vd1^2 - 2*T1*h*vn1^2 - 4*m1*vd1^2*s1 - 4*m1*ve1^2*s1 - 4*m1*vn1^2*s1)/(4*m1*s1^2^(3/2));
  jac_block(5,6) = h*vd1*(-2*A*ae1*c_L1*rho1*vd1^2*s1 - 2*A*ae1*c_L1*rho1*ve1^2*s1 - 2*A*ae1*c_L1*rho1*vn1^2*s1 + A*c_D1*rho1*vd1^2*ve1 + A*c_D1*rho1*ve1^3 + A*c_D1*rho1*ve1*vn1^2 + 2*T1*ve1)/(4*m1*s1^2^(3/2));
  jac_block(5,8) = -A*c_L1*h*rho1*s1^2/(4*m1);
  jac_block(5,10) = -A*ae1*h*rho1*s1^2/(4*m1);
  jac_block(5,14) = h*vn2*(-2*A*ae2*c_L2*rho2*vd2^2*s2 - 2*A*ae2*c_L2*rho2*ve2^2*s2 - 2*A*ae2*c_L2*rho2*vn2^2*s2 + A*c_D2*rho2*vd2^2*ve2 + A*c_D2*rho2*ve2^3 + A*c_D2*rho2*ve2*vn2^2 + 2*T2*ve2)/(4*m2*s2^2^(3/2));
  jac_block(5,15) = (-2*A*ae2*c_L2*h*rho2*vd2^2*ve2*s2 - 2*A*ae2*c_L2*h*rho2*ve2^3*s2 - 2*A*ae2*c_L2*h*rho2*ve2*vn2^2*s2 + A*c_D2*h*rho2*vd2^4 + 3*A*c_D2*h*rho2*vd2^2*ve2^2 + 2*A*c_D2*h*rho2*vd2^2*vn2^2 + 2*A*c_D2*h*rho2*ve2^4 + 3*A*c_D2*h*rho2*ve2^2*vn2^2 + A*c_D2*h*rho2*vn2^4 - 2*T2*h*vd2^2 - 2*T2*h*vn2^2 + 4*m2*vd2^2*s2 + 4*m2*ve2^2*s2 + 4*m2*vn2^2*s2)/(4*m2*s2^2^(3/2));
  jac_block(5,16) = h*vd2*(-2*A*ae2*c_L2*rho2*vd2^2*s2 - 2*A*ae2*c_L2*rho2*ve2^2*s2 - 2*A*ae2*c_L2*rho2*vn2^2*s2 + A*c_D2*rho2*vd2^2*ve2 + A*c_D2*rho2*ve2^3 + A*c_D2*rho2*ve2*vn2^2 + 2*T2*ve2)/(4*m2*s2^2^(3/2));
  jac_block(5,18) = -A*c_L2*h*rho2*s2^2/(4*m2);
  jac_block(5,20) = -A*ae2*h*rho2*s2^2/(4*m2);
  jac_block(6,4) = h*vn1*(-2*A*ad1*c_L1*rho1*vd1^2*s1 - 2*A*ad1*c_L1*rho1*ve1^2*s1 - 2*A*ad1*c_L1*rho1*vn1^2*s1 + A*c_D1*rho1*vd1^3 + A*c_D1*rho1*vd1*ve1^2 + A*c_D1*rho1*vd1*vn1^2 + 2*T1*vd1)/(4*m1*s1^2^(3/2));
  jac_block(6,5) = h*ve1*(-2*A*ad1*c_L1*rho1*vd1^2*s1 - 2*A*ad1*c_L1*rho1*ve1^2*s1 - 2*A*ad1*c_L1*rho1*vn1^2*s1 + A*c_D1*rho1*vd1^3 + A*c_D1*rho1*vd1*ve1^2 + A*c_D1*rho1*vd1*vn1^2 + 2*T1*vd1)/(4*m1*s1^2^(3/2));
  jac_block(6,6) = (-2*A*ad1*c_L1*h*rho1*vd1^3*s1 - 2*A*ad1*c_L1*h*rho1*vd1*ve1^2*s1 - 2*A*ad1*c_L1*h*rho1*vd1*vn1^2*s1 + 2*A*c_D1*h*rho1*vd1^4 + 3*A*c_D1*h*rho1*vd1^2*ve1^2 + 3*A*c_D1*h*rho1*vd1^2*vn1^2 + A*c_D1*h*rho1*ve1^4 + 2*A*c_D1*h*rho1*ve1^2*vn1^2 + A*c_D1*h*rho1*vn1^4 - 2*T1*h*ve1^2 - 2*T1*h*vn1^2 - 4*m1*vd1^2*s1 - 4*m1*ve1^2*s1 - 4*m1*vn1^2*s1)/(4*m1*s1^2^(3/2));
  jac_block(6,9) = -A*c_L1*h*rho1*s1^2/(4*m1);
  jac_block(6,10) = -A*ad1*h*rho1*s1^2/(4*m1);
  jac_block(6,14) = h*vn2*(-2*A*ad2*c_L2*rho2*vd2^2*s2 - 2*A*ad2*c_L2*rho2*ve2^2*s2 - 2*A*ad2*c_L2*rho2*vn2^2*s2 + A*c_D2*rho2*vd2^3 + A*c_D2*rho2*vd2*ve2^2 + A*c_D2*rho2*vd2*vn2^2 + 2*T2*vd2)/(4*m2*s2^2^(3/2));
  jac_block(6,15) = h*ve2*(-2*A*ad2*c_L2*rho2*vd2^2*s2 - 2*A*ad2*c_L2*rho2*ve2^2*s2 - 2*A*ad2*c_L2*rho2*vn2^2*s2 + A*c_D2*rho2*vd2^3 + A*c_D2*rho2*vd2*ve2^2 + A*c_D2*rho2*vd2*vn2^2 + 2*T2*vd2)/(4*m2*s2^2^(3/2));
  jac_block(6,16) = (-2*A*ad2*c_L2*h*rho2*vd2^3*s2 - 2*A*ad2*c_L2*h*rho2*vd2*ve2^2*s2 - 2*A*ad2*c_L2*h*rho2*vd2*vn2^2*s2 + 2*A*c_D2*h*rho2*vd2^4 + 3*A*c_D2*h*rho2*vd2^2*ve2^2 + 3*A*c_D2*h*rho2*vd2^2*vn2^2 + A*c_D2*h*rho2*ve2^4 + 2*A*c_D2*h*rho2*ve2^2*vn2^2 + A*c_D2*h*rho2*vn2^4 - 2*T2*h*ve2^2 - 2*T2*h*vn2^2 + 4*m2*vd2^2*s2 + 4*m2*ve2^2*s2 + 4*m2*vn2^2*s2)/(4*m2*s2^2^(3/2));
  jac_block(6,19) = -A*c_L2*h*rho2*s2^2/(4*m2);
  jac_block(6,20) = -A*ad2*h*rho2*s2^2/(4*m2);
  jac_block(7,4) = an1;
  jac_block(7,5) = ae1;
  jac_block(7,6) = ad1;
  jac_block(7,7) = vn1;
  jac_block(7,8) = ve1;
  jac_block(7,9) = vd1;
  jac_block(8,7) = an1/sqrt(ad1^2 + ae1^2 + an1^2);
  jac_block(8,8) = ae1/sqrt(ad1^2 + ae1^2 + an1^2);
  jac_block(8,9) = ad1/sqrt(ad1^2 + ae1^2 + an1^2);

end

function jac_block = get_end_knot_jac(z_knot, knot_params)
  jac_block = zeros(2, 10);
  pn1 = z_knot(1);
  pe1 = z_knot(2);
  pd1 = z_knot(3);
  vn1 = z_knot(4);
  ve1 = z_knot(5);
  vd1 = z_knot(6);
  an1 = z_knot(7);
  ae1 = z_knot(8);
  ad1 = z_knot(9);
  c_L1 = z_knot(10);

  jac_block(1,4) = an1;
  jac_block(1,5) = ae1;
  jac_block(1,6) = ad1;
  jac_block(1,7) = vn1;
  jac_block(1,8) = ve1;
  jac_block(1,9) = vd1;
  jac_block(2,7) = an1/sqrt(ad1^2 + ae1^2 + an1^2);
  jac_block(2,8) = ae1/sqrt(ad1^2 + ae1^2 + an1^2);
  jac_block(2,9) = ad1/sqrt(ad1^2 + ae1^2 + an1^2);

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
  hes_block = zeros(10, 10);
  h = t_next - t_knot;
  pn1 = z_knot(1);
  pe1 = z_knot(2);
  pd1 = z_knot(3);
  vn1 = z_knot(4);
  ve1 = z_knot(5);
  vd1 = z_knot(6);
  an1 = z_knot(7);
  ae1 = z_knot(8);
  ad1 = z_knot(9);
  c_L1 = z_knot(10);

  lb_c_L = lb(1);
  ub_pd = ub(1);
  ub_c_L = ub(2);

  g1 = knot_params.g;
  m1 = knot_params.m;
  c_D1 = knot_params.c_D;
  rho1 = knot_params.rho;

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

  hes_block(3,3) = (mu/(z_knot(3) - ub(1))^2);
  hes_block(7,7) = (mu/(z_knot(7) - 1.1)^2) + (mu/(-1.1 - z_knot(7))^2);
  hes_block(8,8) = (mu/(z_knot(8) - 1.1)^2) + (mu/(-1.1 - z_knot(8))^2);
  hes_block(9,9) = (mu/(z_knot(9) - 1.1)^2) + (mu/(-1.1 - z_knot(9))^2);
  hes_block(10,10) = (mu/(z_knot(10) - ub(2))^2) + (mu/(lb(1) - z_knot(10))^2);

end

function hes_block = get_middle_knot_hes(z_knot, t_knot, knot_params, lambda_knot, z_next, t_next, lambda_last, lb, ub)
  hes_block = zeros(10, 10);
  h = t_next - t_knot;
  pn1 = z_knot(1);
  pe1 = z_knot(2);
  pd1 = z_knot(3);
  vn1 = z_knot(4);
  ve1 = z_knot(5);
  vd1 = z_knot(6);
  an1 = z_knot(7);
  ae1 = z_knot(8);
  ad1 = z_knot(9);
  c_L1 = z_knot(10);

  lb_c_L = lb(1);
  ub_pd = ub(1);
  ub_c_L = ub(2);

  g1 = knot_params.g;
  m1 = knot_params.m;
  c_D1 = knot_params.c_D;
  rho1 = knot_params.rho;

  lambda1 = lambda_knot(1);
  lambda2 = lambda_knot(2);
  lambda3 = lambda_knot(3);
  lambda4 = lambda_knot(4);
  lambda5 = lambda_knot(5);
  lambda6 = lambda_knot(6);
  lambda7 = lambda_knot(7);
  lambda8 = lambda_knot(8);

  hes_block(3,3) = (mu/(z_knot(3) - ub(1))^2);
  hes_block(7,7) = (mu/(z_knot(7) - 1.1)^2) + (mu/(-1.1 - z_knot(7))^2);
  hes_block(8,8) = (mu/(z_knot(8) - 1.1)^2) + (mu/(-1.1 - z_knot(8))^2);
  hes_block(9,9) = (mu/(z_knot(9) - 1.1)^2) + (mu/(-1.1 - z_knot(9))^2);
  hes_block(10,10) = (mu/(z_knot(10) - ub(2))^2) + (mu/(lb(1) - z_knot(10))^2);
end

function hes_block = get_end_knot_hes(z_knot, t_knot, knot_params, lambda_knot, t_next, lambda_last, lb, ub)
  hes_block = zeros(10, 10);
  h = t_next - t_knot;
  pn1 = z_knot(1);
  pe1 = z_knot(2);
  pd1 = z_knot(3);
  vn1 = z_knot(4);
  ve1 = z_knot(5);
  vd1 = z_knot(6);
  an1 = z_knot(7);
  ae1 = z_knot(8);
  ad1 = z_knot(9);
  c_L1 = z_knot(10);

  lb_c_L = lb(1);
  ub_pd = ub(1);
  ub_c_L = ub(2);

  g1 = knot_params.g;
  m1 = knot_params.m;
  c_D1 = knot_params.c_D;
  rho1 = knot_params.rho;

  lambda1 = lambda_knot(1);
  lambda2 = lambda_knot(2);

  hes_block(1,1) = 2.0;
  hes_block(2,2) = 2.0;
  hes_block(3,3) = 2.0 + (mu/(z_knot(3) - ub(1))^2);
  hes_block(7,7) = (mu/(z_knot(7) - 1.1)^2) + (mu/(-1.1 - z_knot(7))^2);
  hes_block(8,8) = (mu/(z_knot(8) - 1.1)^2) + (mu/(-1.1 - z_knot(8))^2);
  hes_block(9,9) = (mu/(z_knot(9) - 1.1)^2) + (mu/(-1.1 - z_knot(9))^2);
  hes_block(10,10) = (mu/(z_knot(10) - ub(2))^2) + (mu/(lb(1) - z_knot(10))^2);

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
