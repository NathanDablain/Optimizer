##Use broyden for jacobian and broyden-fletcher-goldfarb-shanno for the hessian
% We are solving : |H J'||delta_z|   |f(z)|
%                  |J 0||lambda|  = |c(z)|
% f(z) is the cost gradient vector and has dimension Mx1 where M is the number of decision points
% it should be evaluated analytically where you take the partial derivative of the cost function
% wrt each decision point
% c(z) is just the constraints evaluated with the current decision points, it has dimension Nx1
% where N is the number of constraints
##H(:,:,k+1) = H(:,:,k) + (delta_y(:,k)*delta_y(:,k)')/(delta_y(:,k)'*delta_z(:,k)) -
##            (H(:,:,k)*delta_z(:,k)*delta_z(:,k)'*H(:,:,k))/(delta_z(:,k)'*H(:,:,k)*delta_z(:,k));
##delta_z = z(:,k+1) - z(:,k)
##delta_y = delta_z*lagrange(z_k+1,lambda_k+1) - delta_z*lagrange(z_k, lambda_k+1)
##delta_c = constraints(z(:,k+1)) - constraints(z(:,k));
##J(:,:,k+1) = J(:,:,k) + ((delta_c(:,k) - J(:,:,k)*delta_z(:,k))/(norm(delta_z(:,k))^2))*delta_z(:,k)'

##lagrange(z,lambda) = f(z) + J(z)'*lambda
