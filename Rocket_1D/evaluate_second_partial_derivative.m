function partial_block = evaluate_second_partial_derivative(states,inputs,lambdas,h)
  partial_block = zeros(4,4);

  % Use simplified dynamics for computing partial derivatives
  rho0 = 1.225;
  r0 = 8500;
##  rho = 1.225*exp(-state(1)/8500);
  C_D = 0.235;
  C   = 3000.0;
  A = 10.52;

  r = states(1);
  v = states(2);
  m = states(3);
  T = inputs(1);
  lambda2 = lambdas(2);

  % Should return a matrix of states+inputs by states+inputs
  % First row is partial of first state followed by partial of each state and input
  % Second row is partial of second state followed by partial of each state and input
  % ...
  % we are left with lambda_1*(-0.5*h*dx_1) + lambda_2*(-0.5*h*dx_2) + lambda_3*(-0.5*h*dx_3)
  % This becomes:
  % lambda_1*(-0.5*h*v) +...
  % lambda_2*(-0.5*h*((T - 0.5*C_D*rho*A*v^2)/m - gravity)) +...
  % lambda_3*(-0.5*h*(-T / C))
  % For smooth functions the off diagonal terms should equal

  partial_block(1,1) = (A*C_D*h*lambda2*rho0*(v^2)*exp(-r/r0))/(4*m*(r0^2));
  partial_block(1,2) = (-A*C_D*h*lambda2*rho0*v*exp(-r/r0))/(2*m*r0);
  partial_block(1,3) = (A*C_D*h*lambda2*rho0*(v^2)*exp(-r/r0))/(4*(m^2)*r0);
  partial_block(1,4) = 0;

  partial_block(2,1) = (-A*C_D*h*lambda2*rho0*v*exp(-r/r0))/(2*m*r0); % Same as 1,2
  partial_block(2,2) = (A*C_D*h*lambda2*rho0*exp(-r/r0))/(2*m);
  partial_block(2,3) = (-A*C_D*h*lambda2*rho0*v*exp(-r/r0))/(2*(m^2));
  partial_block(2,4) = 0;

  partial_block(3,1) = (A*C_D*h*lambda2*rho0*(v^2)*exp(-r/r0))/(4*(m^2)*r0); % same as 1,3
  partial_block(3,2) = (-A*C_D*h*lambda2*rho0*v*exp(-r/r0))/(2*(m^2)); % Same as 2,3
  partial_block(3,3) = (h*lambda2*(A*C_D*rho0*(v^2) - 2*T*exp(r/r0))*exp(-r/r0))/(2*(m^3));
  partial_block(3,4) = (h*lambda2)/(2*(m^2));

  partial_block(4,1) = 0;
  partial_block(4,2) = 0;
  partial_block(4,3) = (h*lambda2)/(2*(m^2));
  partial_block(4,4) = 0;
end
