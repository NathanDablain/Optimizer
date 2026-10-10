function O = Simulate_Rocket(O)

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

function jac = get_jacobian(z_next, next_params, h)
  jac = zeros(6, 6);

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
  s2 = sqrt(vd2^2 + ve2^2 + vn2^2);

  jac(1,1) = 1;
  jac(1,4) = -h/2;
  jac(2,2) = 1;
  jac(2,5) = -h/2;
  jac(3,3) = 1;
  jac(3,6) = -h/2;
  jac(4,4) = (h*(vn2*(A*an2*c_L2*rho2*s2^2^(3/2) - vn2*(A*c_D2*rho2*s2^2 - 2*T2)) + s2^2*(-3*A*an2*c_L2*rho2*vn2*s2 + 2*A*c_D2*rho2*vn2^2 + A*c_D2*rho2*s2^2 - 2*T2))*s2/4 + m2*s2^2^2)/(m2*s2^2^2);
  jac(4,5) = h*ve2*(A*an2*c_L2*rho2*s2^2^(3/2) + A*rho2*(-3*an2*c_L2*s2 + 2*c_D2*vn2)*s2^2 - vn2*(A*c_D2*rho2*s2^2 - 2*T2))/(4*m2*s2^2^(3/2));
  jac(4,6) = h*vd2*(A*an2*c_L2*rho2*s2^2^(3/2) + A*rho2*(-3*an2*c_L2*s2 + 2*c_D2*vn2)*s2^2 - vn2*(A*c_D2*rho2*s2^2 - 2*T2))/(4*m2*s2^2^(3/2));
  jac(5,4) = h*vn2*(A*ae2*c_L2*rho2*s2^2^(3/2) + A*rho2*(-3*ae2*c_L2*s2 + 2*c_D2*ve2)*s2^2 - ve2*(A*c_D2*rho2*s2^2 - 2*T2))/(4*m2*s2^2^(3/2));
  jac(5,5) = (h*(ve2*(A*ae2*c_L2*rho2*s2^2^(3/2) - ve2*(A*c_D2*rho2*s2^2 - 2*T2)) + s2^2*(-3*A*ae2*c_L2*rho2*ve2*s2 + 2*A*c_D2*rho2*ve2^2 + A*c_D2*rho2*s2^2 - 2*T2))*s2/4 + m2*s2^2^2)/(m2*s2^2^2);
  jac(5,6) = h*vd2*(A*ae2*c_L2*rho2*s2^2^(3/2) + A*rho2*(-3*ae2*c_L2*s2 + 2*c_D2*ve2)*s2^2 - ve2*(A*c_D2*rho2*s2^2 - 2*T2))/(4*m2*s2^2^(3/2));
  jac(6,4) = h*vn2*(-2*A*ad2*c_L2*rho2*vd2^2*s2 - 2*A*ad2*c_L2*rho2*ve2^2*s2 - 2*A*ad2*c_L2*rho2*vn2^2*s2 + A*c_D2*rho2*vd2^3 + A*c_D2*rho2*vd2*ve2^2 + A*c_D2*rho2*vd2*vn2^2 + 2*T2*vd2)/(4*m2*s2^2^(3/2));
  jac(6,5) = h*ve2*(-2*A*ad2*c_L2*rho2*vd2^2*s2 - 2*A*ad2*c_L2*rho2*ve2^2*s2 - 2*A*ad2*c_L2*rho2*vn2^2*s2 + A*c_D2*rho2*vd2^3 + A*c_D2*rho2*vd2*ve2^2 + A*c_D2*rho2*vd2*vn2^2 + 2*T2*vd2)/(4*m2*s2^2^(3/2));
  jac(6,6) = (h*(vd2*s2^2*(A*ad2*c_L2*rho2*s2^2^(3/2) + 2*g2*m2*s2 - vd2*(A*c_D2*rho2*s2^2 - 2*T2)) + (-2*g2*m2*vd2 + s2*(-3*A*ad2*c_L2*rho2*vd2*s2 + 2*A*c_D2*rho2*vd2^2 + A*c_D2*rho2*s2^2 - 2*T2))*s2^2^(3/2))/4 + m2*s2^2^(5/2))/(m2*s2^2^(5/2));

end

function func = get_function(z_knot, h, knot_params, z_next, next_params)

  dx1 = get_dx(z_knot, knot_params);
  dx2 = get_dx(z_next, next_params);

  func = z_next(1:6) - z_knot(1:6) - 0.5*h*(dx1' + dx2');

end

params_knot = struct('T', 0, 'm', 0, 'g', 0, 'rho', 0, 'c_D', 0,...
                 'rho_part_h', 0, 'c_D_part_s', 0,...
                 'rho_part_h2', 0, 'c_D_part_s2', 0);
params_next = struct('T', 0, 'm', 0, 'g', 0, 'rho', 0, 'c_D', 0,...
                 'rho_part_h', 0, 'c_D_part_s', 0,...
                 'rho_part_h2', 0, 'c_D_part_s2', 0);

for i = 1:O.N_knots-1
  knot_start = O.knot_size*(i-1) + 1;
  knot_end = O.knot_size*i;
  z_knot = O.z(knot_start:knot_end);
  params_knot = get_knot_params(z_knot, O.t(i));
  z_next = O.z(knot_end+1:O.knot_size*(i+1));
  z_next(1:6) = z_knot(1:6);
  params_next = get_knot_params(z_next, O.t(i+1));
  dt = O.t(i+1) - O.t(i);
  for j = 1:5
    jac = get_jacobian(z_next, params_next, dt);
    func = get_function(z_knot, dt, params_knot, z_next, params_next);
    delta_z = jac \ -func;
    z_next(1:6) = z_next(1:6) + delta_z;
    params_next = get_knot_params(z_next, O.t(i+1));
  end
  temp = cross(z_next(4:6), [0; 0; -1]);
  temp = temp ./ norm(temp);
  z_next(7:9) = temp;
  z_next(10) = 0.01;
  O.z(knot_end+1:O.knot_size*(i+1)) = z_next;
end
end
