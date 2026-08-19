##Use broyden for jacobian and broyden-fletcher-goldfarb-shanno for the hessian
% We are solving : |H J'||delta_z|   |f(z)|
%                  |J 0 ||lambda|  = |c(z)|
% f(z) is the cost gradient vector and has dimension Mx1 where M is the number of decision points
% it should be evaluated analytically where you take the partial derivative of the cost function
% wrt each decision point
% c(z) is just the constraints evaluated with the current decision points, it has dimension Nx1
% where N is the number of constraints
##H(:,:,k+1) = H(:,:,k) + (delta_y(:,k)*delta_y(:,k)')/(delta_y(:,k)'*delta_z(:,k)) -
##            (H(:,:,k)*delta_z(:,k)*delta_z(:,k)'*H(:,:,k))/(delta_z(:,k)'*H(:,:,k)*delta_z(:,k));
##delta_z = z(:,k+1) - z(:,k) = s
##delta_y = delta_z*lagrange(z_k+1,lambda_k+1) - delta_z*lagrange(z_k, lambda_k+1) = r
##delta_c = constraints(z(:,k+1)) - constraints(z(:,k));
##J(:,:,k+1) = J(:,:,k) + ((delta_c(:,k) - J(:,:,k)*delta_z(:,k))/(norm(delta_z(:,k))^2))*delta_z(:,k)'

##lagrange(z,lambda) = f(z) + J(z)'*lambda
function [z_opt, cost_opt] = kkt(cost_func, z0, lb, ub, constraints)
  pf_tolerance = 1.0e-6;
  df_tolerance = 1.0e-6;
  max_iterations = 100;
  perturbation = 1.0e-8;
  z_opt = z0;
  z = z0;
  M = length(z0);
  [ineq, eq] = constraints(z0);
  [cost_opt, grad] = cost_func(z0);
  N = length(eq);

  H = eye(M);
  J = zeros(N,M);
  for i = 1:M
    z_perturbed = z0;
    z_perturbed(i) = z_perturbed(i) + perturbation;
    [~, eq_perturbed] = constraints(z_perturbed);
    J(:,i) = (eq_perturbed - eq) ./ perturbation;
  end
  O = M+N;

  for i = 1 : max_iterations

    A = zeros(O,O);
    A(1:M,1:M) = H;
    A(1:M,(M+1):end) = -J';
    A((M+1):end,1:M) = J;
    b = [-grad;-eq];
    x = A \ b;

    delta_z = x(1:M,1);
    lambda_new = x(M+1:end,1);
    z_new = z + delta_z;

    s = delta_z;
    [~, eq_new] = constraints(z_new);
    delta_c = eq_new - eq;
    [~, grad_new] = cost_func(z_new);
    J_new = J + ((delta_c - J*s)*s')./(s'*s);

    r = (grad_new - grad) - (J_new - J)' * lambda_new;

    H_new = H + (r*r')/(r'*s) - (H*s*s'*H)/(s'*H*s);

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
      break;
    end
  end
 end


