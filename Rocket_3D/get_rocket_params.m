function [gravity, rho, c_D, mass, Thrust] = get_rocket_params(h, s, t)

  gravity0 = 9.8065;
  Re = 6371e3;
  gravity = gravity0*((Re/(Re + h))^2);

  if h < 11000
    Temperature = 15.04 - 0.00649*h;
    Pressure = 101.29*(((Temperature+273.1)/288.08)^5.256);
  elseif h >= 11000 && h < 25000
    Temperature = -56.46;
    Pressure = 22.65*exp(1.73 - 0.000157*h);
  else
    Temperature = -131.21 + 0.00299*h;
    Pressure = 2.488*(((Temperature+273.1)/216.6)^-11.388);
  end

  rho = Pressure/(0.2869*(Temperature+273.1));
  speed_of_sound = sqrt(1.4*287*(Temperature+273.1));
  Mach = s / speed_of_sound;

  if Mach < 0.8
    c_D = 0.22;
  elseif Mach >= 0.8 && Mach < 1.2
    c_D = 0.22 + 0.48*sin(pi*(Mach - 0.8)/0.8)^2;
  elseif Mach >= 1.2
    c_D = 0.25 + 0.54/(Mach^1.2);
  end

  % Made up thrust mass tables to take us supersonic
  t_table = [0 0.2  0.5  2.5 3   3.25 4   6   8  10  11  12   13 13.5];
  T_table = [0 300 1000 1000 800 600 550 525 500 450 350 250 100 0];
  m_table = [15 14.92 14.52 11.8533 11.32 11.12 10.57 9.17 7.8367 6.6367 6.17 5.8367 5.7033 5.7033];
  if t >= t_table(end)
    Thrust = T_table(end);
    mass = m_table(end);
  else
    Thrust = interp1(t_table, T_table, t);
    mass = interp1(t_table, m_table, t);
  end

end

