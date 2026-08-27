function partial_block = evaluate_partial_derivative(states, inputs, h)

  % Take partial wrt r_k, v_k, m_k, T_k, r_k+1, v_k+1, m_k+1, T_k+1
  % c = z_k+1 - z_k - 0.5*h*(f_k + f_k+1)
  % c_1 = r_k+1 - r_k - 0.5*h*(v_k + v_k+1)
  % c_2 = v_k+1 - v_k - 0.5*h*((T_k - 0.5*C_D*A*rho_k*v_k^2)/m_k + (T_k+1 - 0.5*C_D*A*rho_k+1*v_k+1^2)/m_k+1)
  % c_3 = m_k+1 - m_k - 0.5*h*((-T_k / C) + (-T_k+1 / C))

  r1 = states(1,1);
  v1 = states(1,2);
  m1 = states(1,3);
  T1 = inputs(1,1);
  r2 = states(2,1);
  v2 = states(2,2);
  m2 = states(2,3);
  T2 = inputs(2,1);

  rho0 = 1.225;
  r0   = 8500.0;
  C_D  = 0.235;
  C    = 3000.0;
  A    = 10.52;

  partial_block = zeros(3, 8);

  partial_block(1,1) = -1;
  partial_block(1,2) = -0.5*h;
  partial_block(1,3) = 0;
  partial_block(1,4) = 0;
  partial_block(1,5) = 1;
  partial_block(1,6) = -0.5*h;
  partial_block(1,7) = 0;
  partial_block(1,8) = 0;

  partial_block(2,1) = (-A*C_D*h*rho0*(v1^2)*exp(-r1/r0))/(4*m1*r0);
  partial_block(2,2) = ((A*C_D*h*rho0*v1*exp(-r1/r0))/(2*m1)) - 1;
  partial_block(2,3) = (h*(-A*C_D*rho0*(v1^2) + 2*T1*exp(r1/r0))*exp(-r1/r0))/(4*(m1^2));
  partial_block(2,4) = -h/(2*m1);
  partial_block(2,5) = (-A*C_D*h*rho0*(v2^2)*exp(-r2/r0))/(4*m2*r0);
  partial_block(2,6) = ((A*C_D*h*rho0*v2*exp(-r2/r0))/(2*m2)) + 1;
  partial_block(2,7) = (h*(-A*C_D*rho0*(v2^2) + 2*T2*exp(r2/r0))*exp(-r2/r0))/(4*(m2^2));
  partial_block(2,8) = -h/(2*m2);

  partial_block(3,1) = 0;
  partial_block(3,2) = 0;
  partial_block(3,3) = -1;
  partial_block(3,4) = 0.5*h/C;
  partial_block(3,5) = 0;
  partial_block(3,6) = 0;
  partial_block(3,7) = 1;
  partial_block(3,8) = 0.5*h/C;

 end