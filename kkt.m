function [z_opt, cost_opt] = kkt(cost_func, z0, cons, constraints, analytical_jac, analytical_hes)
  pf_tolerance = 1.0;
  df_tolerance = 1.0;
  max_iterations = 100;
  perturbation = 1.0e-8;
  z_opt = z0;
  z = z0;
  M = length(z0);
  if analytical_jac
    [ineq, eq, J] = constraints(z0);
  else
    [ineq, eq, ~] = constraints(z0);
  end
  [cost_opt, grad, ~] = cost_func(z0);
  N = length(eq);

  lambda = zeros(N,1);
  if analytical_hes
    H = evaluate_hessian(z, lambda, cons);
  else
    H = eye(M);
  end

  if ~analytical_jac
    J = zeros(N,M);
    for i = 1:M
      z_perturbed = z0;
      z_perturbed(i) = z_perturbed(i) + perturbation;
      [~, eq_perturbed, ~] = constraints(z_perturbed);
      J(:,i) = (eq_perturbed - eq) ./ perturbation;
    end
  end
  O = M+N;

  % When taking a full step would push you over bound, calc alpha = (bound - current)/step
  % Apply this alpha ratio to all steps for the iteration, move out of working set

  for i = 1 : max_iterations


    A = zeros(O,O);
    A(1:M,1:M) = H;
    A(1:M,(M+1):end) = J';
    A((M+1):end,1:M) = J;
    b = [-grad;-eq];
    tic
    x = A \ b;
    toc

    delta_z = x(1:M,1);
    lambda_new = x(M+1:end,1);
    alpha = select_step_size(z, delta_z, cons, cost_func);

    delta_lambda = lambda_new - lambda;
    z_new = z + alpha*delta_z;
    lambda_new = lambda + alpha*delta_lambda;

    s = delta_z;
    if analytical_jac
      [~, eq_new, J_new] = constraints(z_new);
    else
      [~, eq_new, ~] = constraints(z_new);
      delta_c = eq_new - eq;
      J_new = J + ((delta_c - J*s)*s')./(s'*s);
    end
    [~, grad_new, ~] = cost_func(z_new);

    if analytical_hes
      H_new = evaluate_hessian(z_new, lambda_new, cons);
    else
      r = (grad_new - grad) - (J_new - J)' * lambda_new;
      H_new = H + (r*r')/(r'*s) - (H*s*s'*H)/(s'*H*s);
    end


    lambda = lambda_new;
    z      = z_new;
    eq     = eq_new;
    grad   = grad_new;
    H      = H_new;
    J      = J_new;

    % Primal feasability checks if we are satisfying our constraints
    pf = max(abs(eq));
    if pf < pf_tolerance
      pf_satisfied = true;
    else
      pf_satisfied = false;
    end
    % Dual feasability checks if we are at a local minimum while balancing our constraints
    lagrange_grad = (grad - J'*lambda);
    df = max(abs(lagrange_grad));
    if df < df_tolerance
      df_satisfied = true;
    else
      df_satisfied = false;
    end

    if pf_satisfied && df_satisfied
      z_opt = z;
      [cost_opt, ~] = cost_func(z_opt);
      disp(['Converged in ' num2str(i) ' iterations'])
      break;
    elseif i == max_iterations
      disp(['Failed to converge with pf of ' num2str(pf) ' and df of ' num2str(df)])
    end
  end
 end

