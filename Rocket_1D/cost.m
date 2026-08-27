function [J, grad, double_grad] = cost(z, cons)
    grad = zeros(length(z),1);
    double_grad = zeros(length(z),length(z));

    final_knot = (cons.N_states + cons.N_inputs) * cons.N_knots;
    id_s1 = final_knot - 3;
    id_s2 = final_knot - 2;

    rf = cons.xd(1);
    J = (z(id_s1) - rf)^2;

    grad(id_s1) = 2.0*z(id_s1) - 2*(rf);

  for i = final_knot+1:length(z)
    J = J - log(z(i));
    grad(i) = -1 / z(i);
    double_grad(i,i) = 1/(z(i)^2);
  end

  double_grad(id_s1, id_s1) = 2.0;
end