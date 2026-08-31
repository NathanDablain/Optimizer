function dx = diffeq(state, input, t)
  % We have states of position, speed, heading, fpa
  % We have inputs of hdg lift coefficent and fpa lift coefficient
  dx = zeros(6,1);

  r = state(1:3);
  s = state(4);
  psi = state(5);
  gam = state(6);

  Clp = input(1);
  Clg = input(2);

  [gravity, rho, c_D, m, Thrust] = get_rocket_params(-r(3), s, t);

  A = 1.52;
  Q = 0.5 * rho * s * s;
  Drag = Q * c_D * A;

  dx(1) = s * cos(gam) * cos(psi);
  dx(2) = s * cos(gam) * sin(psi);
  dx(3) = -s * sin(gam);
  dx(4) = (Thrust - Drag)/m - gravity*sin(gam);
  dx(5) = (rho * A * s * Clp) / (m * 2 * cos(gam));
  dx(6) = (0.5 * rho * A * s * Clg)/m - (gravity * cos(gam)) / s;

end


