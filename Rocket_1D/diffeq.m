% Specific function to problem being solved, describes the continuous time differential equations
function dX = diffeq(state, input)
  gravity = 9.81;
  C = 3000.0;
  A = 10.52;

  if state(1) < 11000
    Temperature = 15.04 - 0.00649*state(1);
    Pressure = 101.29*(((Temperature+273.1)/288.08)^5.256);
  elseif state(1) >= 11000 && state(1) < 25000
    Temperature = -56.46;
    Pressure = 22.65*exp(1.73 - 0.000157*state(1));
  else
    Temperature = -131.21 + 0.00299*state(1);
    Pressure = 2.488*(((Temperature+273.1)/216.6)^-11.388);
  endif
  rho = Pressure/(0.2869*(Temperature+273.1));
  speed_of_sound = sqrt(1.4*287*(Temperature+273.1));
  Mach = state(2) / speed_of_sound;
  if Mach < 0.8
    C_D = 0.22;
  elseif Mach >= 0.8 && Mach < 1.2
    C_D = 0.22 + 0.48*sin(pi*(Mach - 0.8)/0.8)^2;
  elseif Mach >= 1.2
    C_D = 0.25 + 0.54/(Mach^1.2);
  endif

  Drag = 0.5*C_D*rho*A*state(2)*state(2);

  dX(1) = state(2);
  dX(2) = (input(1) - Drag)/state(3) - gravity;
  dX(3) = -input(1) / C;
  end
