function O = Simulate_Rocket(O)

A = 0.25;
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

function [dx, w] = get_dx(z_knot, knot_params)
  s = z_knot(4);
  q0 = z_knot(5);
  q1 = z_knot(6);
  q2 = z_knot(7);
  q3 = z_knot(8);
  wy = z_knot(9);
  wz = z_knot(10);
  Cly = z_knot(11);
  Clz = z_knot(12);

  q = [q0 q1 q2 q3];
##  q0 = q0 / norm(q);
##  q1 = q1 / norm(q);
##  q2 = q2 / norm(q);
##  q3 = q3 / norm(q);

  D = 0.5*knot_params.rho*A*knot_params.c_D*s*s;
  az = 0.5*knot_params.rho*A*Clz*s*s/knot_params.m;
  ay = 0.5*knot_params.rho*A*Cly*s*s/knot_params.m;

  wy = (-az - knot_params.g*(-2*q1^2 - 2*q2^2 + 1))/s;
  wz = (ay + 2*knot_params.g*(q0*q1 + q2*q3))/s;
  v_NED = [s*(-2*q2^2 - 2*q3^2 + 1); 2*s*(q0*q3 + q1*q2); 2*s*(-q0*q2 + q1*q3)];
  omega = [0 0 -wy -wz; 0 0 wz -wy; wy -wz 0 0; wz wy 0 0];

  dx(1:3) = v_NED;
  dx(4) = ((knot_params.T - D)/knot_params.m) + 2*knot_params.g*(-q0*q2 + q1*q3);
##  dx(5:8) = 0.5*omega*q';
  w = [wy; wz];
end


params(1:O.N_knots) = struct('T', 0, 'm', 0, 'g', 0, 'rho', 0, 'c_D', 0,...
                             'rho_part_h', 0, 'c_D_part_s', 0,...
                             'rho_part_h2', 0, 'c_D_part_s2', 0);

for i = 1:O.N_knots
  knot_start = O.knot_size*(i-1) + 1;
  knot_end = O.knot_size*i;
  O.z(knot_start+10) = 0.0;
  O.z(knot_start+11) = -0.1;
  if i > 1
##    O.z(knot_start:knot_start+7) = z_knot(1:8) + dx'*(O.t(i) - O.t(i-1));
    O.z(knot_start:knot_start+3) = z_knot(1:4) + dx(1:4)'*(O.t(i) - O.t(i-1));
    theta = norm(w)*(O.t(i) - O.t(i-1));
    delta_q = [cos(theta/2);0;(w(1)/norm(w))*sin(theta/2);(w(2)/norm(w))*sin(theta/2)];
    q0 = z_knot(5);
    q1 = z_knot(6);
    q2 = z_knot(7);
    q3 = z_knot(8);
    q_new = [q0 -q1 -q2 -q3; q1 q0 -q3 q2; q2 q3 q0 -q1; q3 -q2 q1 q0]*delta_q;
    O.z(knot_start+4:knot_start+7) = q_new;
##    O.z(knot_start+4:knot_start+7) = O.z(knot_start+4:knot_start+7)./norm(O.z(knot_start+4:knot_start+7));
  end
  z_knot = O.z(knot_start:knot_end);
  params(i) = get_knot_params(z_knot, O.t(i));
  [dx, w] = get_dx(z_knot, params(i));
  O.z(knot_start+8:knot_start+9) = w;
end

end

