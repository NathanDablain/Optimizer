function [J, grad, double_grad] = cost(z, cons)
   %(z - 2)^2 -> z^2 - 4*z + 4
   grad = zeros(length(z),1);
   double_grad = zeros(length(z),length(z));

   final_knot = (cons.N_states + cons.N_inputs) * cons.N_knots;
   id_s1 = final_knot - 2;
   id_s2 = final_knot - 1;

   xd = cons.xd;
   J = (z(id_s1) - xd(1))^2 + (z(id_s2) - xd(2))^2;

   for i = final_knot+1:length(z)
    J = J - log(z(i));
    grad(i) = -1 / z(i);
    double_grad(i,i) = 1/(z(i)^2);
   end

   grad(id_s1) = 2*z(id_s1) - 2*xd(1);
   grad(id_s2) = 2*z(id_s2) - 2*xd(2);

   double_grad(id_s1,id_s1) = 2;
   double_grad(id_s2,id_s2) = 2;
end
