function [states, inputs] = build_traj(t, x0, inputs)

  function f = get_f(x_c, x_f, u_c, u_f, h)
  %   f = x_i+1 - x_i - 0.5*h*(y_i+1 + y_i) = 0
    f = x_f - x_c - 0.5*h*(diffeq(x_f,u_f)' + diffeq(x_c,u_c)');
  endfunction

  % Performs trapezoidal integration using the Broyden method
  N = length(t);
  D = length(x0);
  states = zeros(N, D);
  states(1,:) = x0;

  for i = 2:N
    h = t(i) - t(i-1);

    % Here B is the inverse of our estimate of the jacobian
    % We need an initial guess for B, and y
    % The equation(vector) f we are trying to solve is x_i+1 - x_i - 0.5*h*(y_i+1 + y_i) = 0
    % the variable(vector) x we are solving for is x_i+1
    iterations = 10;
    stop_index = iterations;
    B = zeros(D,D,iterations-1);
    y_guess = zeros(D,iterations);
    f_guess = zeros(D,iterations-1);
    delta_x = zeros(D,iterations-1);
    delta_f = zeros(D,iterations-1);

    B(:,:,1) = eye(D);
    y_guess(:,1) = (states(i-1,:) + h*diffeq(states(i-1,:),inputs(i-1)))';
    f_guess(:,1) = get_f(states(i-1,:)', y_guess(:,1), inputs(i-1), inputs(i), h);
    y_guess(:,2) = y_guess(:,1) - B(:,:,1)*f_guess(:,1);
    for k = 2:iterations-1
      f_guess(:,k) = get_f(states(i-1,:)', y_guess(:,k), inputs(i-1), inputs(i), h);
      delta_x(:,k) = y_guess(:,k) - y_guess(:,k-1);
      if norm(delta_x(:,k)) < 1e-6 && norm(f_guess(:,k)) < 1e-10
        stop_index = k;
        break;
      endif
      delta_f(:,k) = f_guess(:,k) - f_guess(:,k-1);
      B(:,:,k) = B(:,:,k-1) + ((delta_x(:,k) - B(:,:,k-1)*delta_f(:,k)) / ...
                              (delta_x(:,k)'*B(:,:,k-1)*delta_f(:,k)))*delta_x(:,k)'*B(:,:,k-1);
      y_guess(:,k+1) = y_guess(:,k) - B(:,:,k)*f_guess(:,k);
      if k == iterations-1
        disp(['Reached iteration limit on knot index:' num2str(i)]);
      endif
    endfor
    states(i,:) = y_guess(:,stop_index)';
  endfor

endfunction

