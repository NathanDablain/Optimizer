function [J, grad, double_grad] = cost(z, cons, test_command, iteration)
  grad = zeros(length(z),1);
  double_grad = zeros(length(z),length(z));

  final_knot = (cons.N_states + cons.N_inputs) * cons.N_knots;
  id_s1 = final_knot - 7;
  id_s2 = final_knot - 6;
  id_s3 = final_knot - 5;
  id_s4 = final_knot - 4;
  id_s5 = final_knot - 3;
  id_s6 = final_knot - 2;

  J = (z(id_s1) - cons.xd(1))^2 +...
      (z(id_s2) - cons.xd(2))^2 +...
      (z(id_s3) - cons.xd(3))^2;

  grad(id_s1) = 2.0*z(id_s1) - 2*cons.xd(1);
  grad(id_s2) = 2.0*z(id_s2) - 2*cons.xd(2);
  grad(id_s3) = 2.0*z(id_s3) - 2*cons.xd(3);

  mu = 1 / iteration;
  for i = final_knot+1:length(z)
    J = J - log(z(i));
    grad(i) = -mu / z(i);
    if test_command
      double_grad(i,i) = rand;
    else
      double_grad(i,i) = mu/(z(i)^2);
    end
  end

  double_grad(id_s1, id_s1) = 2.0;
  double_grad(id_s2, id_s2) = 2.0;
  double_grad(id_s3, id_s3) = 2.0;

end
