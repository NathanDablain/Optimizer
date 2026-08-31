function partial_block = evaluate_partial_derivative(states, inputs, t, test_command)
t1 = t(1);
t2 = t(2);
h  = t2 - t1;

rn1  = states(1,1);
re1  = states(1,2);
rd1  = states(1,3);
s1   = states(1,4);
psi1 = states(1,5);
gam1 = states(1,6);
Clp1 = inputs(1,1);
Clg1 = inputs(1,2);

rn2  = states(2,1);
re2  = states(2,2);
rd2  = states(2,3);
s2   = states(2,4);
psi2 = states(2,5);
gam2 = states(2,6);
Clp2 = inputs(2,1);
Clg2 = inputs(2,2);

##[gravity, rho, c_D, m, Thrust] = get_rocket_params(-r(3), s, t);

g = 9.8065;
c_D = 0.3;
A = 1.52;
rho0 = 1.225;
r0   = 8500.0;
rho1 = rho0*exp(rd1/r0);
rho2 = rho0*exp(rd2/r0);
D1 = 0.5 * A * c_D * rho1 * s1 * s1;
D2 = 0.5 * A * c_D * rho2 * s2 * s2;
[greal1, rhoreal1, c_Dreal1, m1, Thrust1] = get_rocket_params(-rd1, s1, t1);
[greal2, rhoreal2, c_Dreal2, m2, Thrust2] = get_rocket_params(-rd2, s2, t2);


if test_command
  partial_block = rand(6,16);
else
  % wrt position - north
  partial_block(1,1) = -1.0;
  partial_block(1,4) = -0.5 * h * cos(gam1) * cos(psi1);
  partial_block(1,5) = 0.5 * h * s1 * sin(psi1) * cos(gam1);
  partial_block(1,6) = 0.5 * h * s1 * sin(gam1) * cos(psi1);
  partial_block(1,9) = 1.0;
  partial_block(1,12) = -0.5 * h * cos(gam2) * cos(psi2);
  partial_block(1,13) = 0.5 * h * s2 * sin(psi2) * cos(gam2);
  partial_block(1,14) = 0.5 * h * s2 * sin(gam2) * cos(psi2);

  % wrt position - east
  partial_block(2,2) = -1.0;
  partial_block(2,4) = -0.5 * h * sin(psi1) * cos(gam1);
  partial_block(2,5) = -0.5 * h * s1 * cos(gam1) * cos(psi1);
  partial_block(2,6) = 0.5 * h * s1 * sin(gam1) * sin(psi1);
  partial_block(2,10) = 1.0;
  partial_block(2,12) = -0.5 * h * sin(psi2) * cos(gam2);
  partial_block(2,13) = -0.5 * h * s2 * cos(gam2) * cos(psi2);
  partial_block(2,14) = 0.5 * h * s2 * sin(gam2) * sin(psi2);

  % wrt position - down
  partial_block(3,3) = -1.0;
  partial_block(3,4) = 0.5 * h * sin(gam1);
  partial_block(3,6) = 0.5 * h * s1 * cos(gam1);
  partial_block(3,11) = 1.0;
  partial_block(3,12) = 0.5 * h * sin(gam2);
  partial_block(3,14) = 0.5 * h * s2 * cos(gam2);

  % wrt speed
  partial_block(4,3) = h*D1/(2*m1*r0);
  partial_block(4,4) = (h*(D1/s1) - m1)/m1;
  partial_block(4,6) = 0.5*h*g*cos(gam1);
  partial_block(4,11) = h*D2/(2*m2*r0);
  partial_block(4,12) = (h*(D2/s2) - m2)/m2;
  partial_block(4,14) = 0.5*h*g*cos(gam2);

  % wrt heading
  partial_block(5,3) = (-A*Clp1*h*rho1*s1)/(4*r0*cos(gam1));
  partial_block(5,4) = (-A*Clp1*h*rho1)/(4*cos(gam1));
  partial_block(5,5) = -1.0;
  partial_block(5,6) = (-A*Clp1*h*rho1*s1*sin(gam1))/(4*(cos(gam1)^2));
  partial_block(5,7) = (-A*h*rho1*s1)/(4*cos(gam1));
  partial_block(5,11) = (-A*Clp2*h*rho2*s2)/(4*r0*cos(gam2));
  partial_block(5,12) = (-A*Clp2*h*rho2)/(4*cos(gam2));
  partial_block(5,13) = 1.0;
  partial_block(5,14) = (-A*Clp2*h*rho2*s2*sin(gam2))/(4*(cos(gam2)^2));
  partial_block(5,15) = (-A*h*rho2*s2)/(4*cos(gam2));

  % wrt flight path angle
  partial_block(6,3) = (-A*Clg1*h*rho1*s1)/(4*r0);
  partial_block(6,4) = ((-A*Clg1*h*rho1)/4) - ((g*h*cos(gam1))/(2*s1*s1));
  partial_block(6,6) = (((-g*h*sin(gam1))/2) - s1)/s1;
  partial_block(6,8) = (-A*h*rho1*s1)/4;
  partial_block(6,11) = (-A*Clg2*h*rho2*s2)/(4*r0);
  partial_block(6,12) = ((-A*Clg2*h*rho2)/4) - ((g*h*cos(gam2))/(2*s2*s2));
  partial_block(6,14) = (((-g*h*sin(gam2))/2) + s2)/s2;
  partial_block(6,16) = (-A*h*rho2*s2)/4;
end


end

