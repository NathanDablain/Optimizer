function O = Simulate_Rocket(O)

A = 10.52;
C = 3000.0;
function knot_params = get_knot_params(z_knot, t_knot)
  gravity0 = 9.8065;
  Re = 6371e3;
  h = z_knot(1);
  v = z_knot(2);

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
  Mach = v / speed_of_sound;
  if Mach < 0.8
    knot_params.c_D = 0.22;

    knot_params.c_D_part_s = 0;
    knot_params.c_D_part_s2 = 0;
  elseif Mach >= 0.8 && Mach < 1.2
    c1 = 0.22;
    c2 = 0.48;
    c3 = 0.8;
    knot_params.c_D = c1 + c2*sin(pi*(Mach - c3)/c3)^2;

    knot_params.c_D_part_s = (pi*c2*sin(2*pi*v/(speed_of_sound*c3)))/(speed_of_sound*c3);
    knot_params.c_D_part_s2 = (2*pi^2*c2*cos((2*pi*v)/(speed_of_sound*c3)))/(speed_of_sound^2*c3^2);
  elseif Mach >= 1.2
    c1 = 0.25;
    c2 = 0.54;
    c3 = 1.2;
    knot_params.c_D = c1 + c2/(Mach^c3);

    knot_params.c_D_part_s = (-c2*c3*((v/speed_of_sound)^-c3))/v;
    knot_params.c_D_part_s2 = (c2*c3*((v/speed_of_sound)^-c3)*(c3 + 1))/(v^2);
  end


end


function dx = get_dx(z_knot, knot_params)
  h = z_knot(1);
  v = z_knot(2);
  m = z_knot(3);
  T = z_knot(4);

  Drag = 0.5*knot_params.c_D*knot_params.rho*A*v*v;

  dx(1) = v;
  dx(2) = (T - Drag)/m - knot_params.g;
  dx(3) = -T / C;
end

params(1:O.N_knots) = struct('g', 0, 'rho', 0, 'c_D', 0,...
                               'g_part_h', 0, 'rho_part_h', 0, 'c_D_part_s', 0,...
                               'g_part_h2', 0, 'rho_part_h2', 0, 'c_D_part_s2', 0);
O.z(1:length(O.ic)) = O.ic;
for i = 1:O.N_knots
  knot_start = O.knot_size*(i-1) + 1;
  knot_end = O.knot_size*i;
  O.z(knot_start+4) = 6.0e6;
  if i > 1
    O.z(knot_start:knot_start+2) = z_knot(1:3) + dx'*(O.t(i) - O.t(i-1));
  end
  z_knot = O.z(knot_start:knot_end);
  params(i) = get_knot_params(z_knot, O.t(i));
  dx = get_dx(z_knot, params(i));
end

end

