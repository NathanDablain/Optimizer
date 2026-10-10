function O = Simulate_Rocket(O)

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

function [dx, w] = get_dx(z_knot, knot_params)
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

  wy = (-az - knot_params.g*(-2*q1^2 - 2*q2^2 + 1))/s;
  wz = (ay + 2*knot_params.g*(q0*q1 + q2*q3))/s;
  v_NED = [s*(-2*q2^2 - 2*q3^2 + 1); 2*s*(q0*q3 + q1*q2); 2*s*(-q0*q2 + q1*q3)];
  omega = [0 0 -wy -wz; 0 0 wz -wy; wy -wz 0 0; wz wy 0 0];

  dx(1:3) = v_NED;
  dx(4) = ((knot_params.T - D)/knot_params.m) + 2*knot_params.g*(-q0*q2 + q1*q3);
##  dx(5:8) = 0.5*omega*q';
  w = [wy; wz];
end

function jac = get_jacobian(z_next, next_params, h)
  jac = zeros(8, 8);

  pd2 = z_next(3);
  s2 = z_next(4);
  q02 = z_next(5);
  q12 = z_next(6);
  q22 = z_next(7);
  q32 = z_next(8);

  g = next_params.g;
  m = next_params.m;
  c_D = next_params.c_D;
  rho = next_params.rho;
  rho_part_h = next_params.rho_part_h;
  c_D_part_s = next_params.c_D_part_s;

  jac(1,1) = 1;
  jac(1,4) = h*(q22^2 + q32^2 - 1/2);
  jac(1,7) = 2*h*q22*s2;
  jac(1,8) = 2*h*q32*s2;
  jac(2,2) = 1;
  jac(2,4) = -h*(q02*q32 + q12*q22);
  jac(2,5) = -h*q32*s2;
  jac(2,6) = -h*q22*s2;
  jac(2,7) = -h*q12*s2;
  jac(2,8) = -h*q02*s2;
  jac(3,3) = 1;
  jac(3,4) = h*(q02*q22 - q12*q32);
  jac(3,5) = h*q22*s2;
  jac(3,6) = -h*q32*s2;
  jac(3,7) = h*q02*s2;
  jac(3,8) = -h*q12*s2;
  jac(4,3) = A*h*s2^2*c_D*rho_part_h/(4*m);
  jac(4,4) = (A*h*s2*(s2*c_D_part_s + 2*c_D)*rho/4 + m)/m;
  jac(4,5) = g*h*q22;
  jac(4,6) = -g*h*q32;
  jac(4,7) = g*h*q02;
  jac(4,8) = -g*h*q12;
  jac(5,5) = 1;
  jac(6,6) = 1;
  jac(7,7) = 1;
  jac(8,8) = 1;

end

function func = get_function(z_knot, h, knot_params, z_next, next_params)

  [dx1, w1] = get_dx(z_knot, knot_params);
  [dx2, w2] = get_dx(z_next, next_params);

  w = 0.5*(w1 + w2);

  if norm(w) > 1.0e-5
    theta = norm(w)*h;
    delta_q = [cos(theta/2);0;(w(1)/norm(w))*sin(theta/2);(w(2)/norm(w))*sin(theta/2)];
  else
    delta_q = [1;0;0;0];
  end
  q_new = quaternion_multiply(z_knot(5:8), delta_q, 'right');

  func = [z_next(1:4) - z_knot(1:4) - 0.5*h*(dx1(1:4)' + dx2(1:4)');...
          z_next(5:8) - q_new];
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
  z_next(1:8) = z_knot(1:8);
  z_next(11) = 0.1;
  z_next(12) = pi;
  params_next = get_knot_params(z_next, O.t(i+1));
  dt = O.t(i+1) - O.t(i);
  for j = 1:5
    jac = get_jacobian(z_next, params_next, dt);
    func = get_function(z_knot, dt, params_knot, z_next, params_next);
    delta_z = jac \ -func;
    z_next(1:8) = z_next(1:8) + delta_z;
    params_next = get_knot_params(z_next, O.t(i+1));
  end
  [~, w1] = get_dx(z_knot, params_knot);
  [~, w2] = get_dx(z_next, params_next);

  w = 0.5*(w1 + w2);
  O.z(knot_start+8:knot_start+9) = w;
  O.z(knot_end+1:O.knot_size*(i+1)) = z_next;
end
end
