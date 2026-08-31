function partial_block = evaluate_second_partial_derivative(state, input, lambda, t, h, test_command)
  if test_command
    partial_block = rand(8,8);
  else
    partial_block = zeros(8,8);

    rn  = state(1);
    re  = state(2);
    rd  = state(3);
    s   = state(4);
    psi = state(5);
    gam = state(6);
    Clp = input(1);
    Clg = input(2);
    lambda1 = lambda(1);
    lambda2 = lambda(2);
    lambda3 = lambda(3);
    lambda4 = lambda(4);
    lambda5 = lambda(5);
    lambda6 = lambda(6);

    g = 9.8065;
    c_D = 0.3;
    A = 1.52;
    rho0 = 1.225;
    r0   = 8500.0;
    rho = rho0*exp(rd/r0);
    [greal, rhoreal, c_Dreal, m, Thrust] = get_rocket_params(-rd, s, t);

    % Row 1 is all zeros

    % Row 2 is all zeros

    % Row 3
    partial_block(3,3) = -0.5*h*((A*rho*s*(Clg*lambda6*m + ((Clp*lambda5*m)/cos(gam)) - c_D*lambda4*s))/(2*m*r0*r0));
    partial_block(3,4) = -0.5*h*((A*rho*(Clg*lambda6*m + ((Clp*lambda5*m)/cos(gam)) - 2*c_D*lambda4*s))/(2*m*r0));
    partial_block(3,6) = -0.5*h*((A*Clp*lambda5*rho*s*sin(gam))/(2*r0*(cos(gam)^2)));
    partial_block(3,7) = -0.5*h*((A*lambda5*rho*s)/(2*r0*cos(gam)));
    partial_block(3,8) = -0.5*h*((A*lambda6*rho*s)/(2*r0));

    % Row 4
    partial_block(4,3) = partial_block(3,4);
    partial_block(4,4) = -0.5*h*(((-A*c_D*lambda4*rho)/m) - ((2*g*lambda6*cos(gam))/(s^3)));
    partial_block(4,5) = -0.5*h*(((-lambda1*sin(psi)) + (lambda2*cos(psi)))*cos(gam));
    partial_block(4,6) = -0.5*h*(((A*Clp*lambda5*rho*sin(gam))/(2*(cos(gam)^2))) - ((g*lambda6*sin(gam))/(s^2)) -...
                          lambda1*sin(gam)*cos(psi) - lambda2*sin(gam)*sin(psi) - lambda3*cos(gam));
    partial_block(4,7) = -0.5*h*((A*lambda5*rho)/(2*cos(gam)));
    partial_block(4,8) = -0.5*h*((A*lambda6*rho)/2);

    % Row 5
    partial_block(5,4) = partial_block(4,5);
    partial_block(5,5) = -0.5*h*(-s*(lambda1*cos(psi) + lambda2*sin(psi))*cos(gam));
    partial_block(5,6) = -0.5*h*(s*(lambda1*sin(psi) - lambda2*cos(psi))*sin(gam));

    % Row 6
    partial_block(6,3) = partial_block(3,6);
    partial_block(6,4) = partial_block(4,6);
    partial_block(6,5) = partial_block(5,6);
    partial_block(6,6) = -0.5*h*(((A*Clp*lambda5*rho*s*(sin(gam)^2))/(cos(gam)^3)) + ((A*Clp*lambda5*rho*s)/(2*cos(gam))) +...
                          g*lambda4*sin(gam) + ((g*lambda6*cos(gam))/s) - lambda1*s*cos(gam)*cos(psi) -...
                          lambda2*s*sin(psi)*cos(gam) + lambda3*s*sin(gam));
    partial_block(6,7) = -0.5*h*(((A*lambda5*rho*s*sin(gam))/(2*(cos(gam)^2))));

    % Row 7
    partial_block(7,3) = partial_block(3,7);
    partial_block(7,4) = partial_block(4,7);
    partial_block(7,6) = partial_block(6,7);

    % Row 8
    partial_block(8,3) = partial_block(3,8);
    partial_block(8,4) = partial_block(4,8);

  end
end
