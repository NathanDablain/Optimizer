clear
clc

% Should be able to supply differential equations, break out states and inputs
% as well as what variables will be bounded, and then be able to build the jacobian
% and hessian from there.

% Make a new ordering that groups knot vars together so slack isnt afterwords but
% instead packed with the states and inputs of a knot

% This model is a double integrator with an upper and lower bound constraint imposed
% on the acceleration command


pkg load symbolic

function [dx, states, vars] = get_dx(p, v, a)
  dx(1) = v;
  dx(2) = a;

  states = [p v];
  inputs = a;
  vars = [states inputs];
end

function c = get_constraints_first_knot(states, inputs, slack, lbs, ubs, ic)
  syms h real;

  p1 = states(1,1);
  v1 = states(1,2);
  a1 = inputs(1,1);
  s_a_lb = slack(1,1);
  s_a_ub = slack(1,2);
  lb_a = lbs(1);
  ub_a = ubs(1);
  p_ic = ic(1);
  v_ic = ic(2);

  p2 = states(2,1);
  v2 = states(2,2);
  a2 = inputs(2,1);

  dx1 = get_dx(p1, v1, a1);
  dx2 = get_dx(p2, v2, a2);
  c = [...
    % the trapezoidal constraints for the knot
    p2 - p1 - 0.5*h*(dx1(1) + dx2(1));...
    v2 - v1 - 0.5*h*(dx1(2) + dx2(2));...
    % the slack variable constraints for the knot
    % first lower bound
    lb_a - a1 + s_a_lb;...
    % then upper bound
    a1 - ub_a + s_a_ub
    % the initial condition constraints for the knot
    p_ic - p1;...
    v_ic - v1];
end

function c = get_constraints_middle_knot(states, inputs, slack, lbs, ubs)
  syms h real;

  p1 = states(1,1);
  v1 = states(1,2);
  a1 = inputs(1,1);
  s_a_lb = slack(1,1);
  s_a_ub = slack(1,2);
  lb_a = lbs(1);
  ub_a = ubs(1);

  p2 = states(2,1);
  v2 = states(2,2);
  a2 = inputs(2,1);

  dx1 = get_dx(p1, v1, a1);
  dx2 = get_dx(p2, v2, a2);
  c = [...
    % the trapezoidal constraints for the knot
    p2 - p1 - 0.5*h*(dx1(1) + dx2(1));...
    v2 - v1 - 0.5*h*(dx1(2) + dx2(2));...
    % the slack variable constraints for the knot
    % first lower bound
    lb_a - a1 + s_a_lb;...
    % then upper bound
    a1 - ub_a + s_a_ub];
end

function c = get_constraints_end_knot(states, inputs, slack, lbs, ubs)
  syms h real;

  p1 = states(1,1);
  v1 = states(1,2);
  a1 = inputs(1,1);
  s_a_lb = slack(1,1);
  s_a_ub = slack(1,2);
  lb_a = lbs(1);
  ub_a = ubs(1);

  c = [...
    % the slack variable constraints for the knot
    % first lower bound
    lb_a - a1 + s_a_lb;...
    % then upper bound
    a1 - ub_a + s_a_ub];
end

function lagrangian = get_lagrangian_first_knot(states, inputs, slack, lbs, ubs, ic)
  syms mu lambda1 lambda2 lambda3 lambda4 lambda5 lambda6 real;

  lambda = [lambda1 lambda2 lambda3 lambda4 lambda5 lambda6];
  cost = -mu*(log(slack(1)) + log(slack(2)));

  lagrangian = cost + lambda*get_constraints_first_knot(states, inputs, slack, lbs, ubs, ic);
end

function lagrangian = get_lagrangian_middle_knot(states, inputs, slack, lbs, ubs)
  syms mu lambda1 lambda2 lambda3 lambda4 real;

  lambda = [lambda1 lambda2 lambda3 lambda4];
  cost = -mu*(log(slack(1)) + log(slack(2)));

  lagrangian = cost + lambda*get_constraints_middle_knot(states, inputs, slack, lbs, ubs);
end

function lagrangian = get_lagrangian_end_knot(states, inputs, slack, lbs, ubs, xd)
  syms mu lambda1 lambda2 real;

  lambda = [lambda1 lambda2];
  cost = (states(1,1) - xd(1))^2 + (states(1,2) - xd(2))^2-mu*(log(slack(1)) + log(slack(2)));

  lagrangian = cost + lambda*get_constraints_end_knot(states, inputs, slack, lbs, ubs);
end

syms p1 v1 a1 p2 v2 a2 s_a_lb s_a_ub lb_a ub_a p_ic v_ic xd1 xd2 real

states = [p1 v1;...
          p2 v2];
inputs = [a1;...
          a2];
slack = [s_a_lb, s_a_ub];
lbs = lb_a;
ubs = ub_a;
ic = [p_ic, v_ic];
xd = [xd1 xd2];

z = [p1 v1 a1 s_a_lb s_a_ub p2 v2 a2];
% All knots have 2 states, 1 input, 2 slack variables -> 5 decision variables
% So the first knot has 6 constraints
cf = get_constraints_first_knot(states, inputs, slack, lbs, ubs, ic);

% The middle knots have 4 constraints
cm = get_constraints_middle_knot(states, inputs, slack, lbs, ubs);

% The end knot has 2 constraints
ce = get_constraints_end_knot(states, inputs, slack, lbs, ubs);
% So the total number of decision variables is 5 * N_knots
% The total number of constraints is 6 + 4*(N_knots - 2) + 2

Lf = get_lagrangian_first_knot(states, inputs, slack, lbs, ubs, ic);
Lm = get_lagrangian_middle_knot(states, inputs, slack, lbs, ubs);
Le = get_lagrangian_end_knot(states, inputs, slack, lbs, ubs, xd);

clc
disp(cf)
for i = 1:length(cf)
  for j = 1:length(z)
    blockcf(i,j) = simplify(diff(cf(i), z(j)));
    if blockcf(i,j) ~= 0
      fprintf('jac_block(%d,%d) = %s;\n',i,j,char(blockcf(i,j)));
    end
  end
end

fprintf('\n\n\n')
disp(cm)
for i = 1:length(cm)
  for j = 1:length(z)
    blockcm(i,j) = simplify(diff(cm(i), z(j)));
    if blockcm(i,j) ~= 0
      fprintf('jac_block(%d,%d) = %s;\n',i,j,char(blockcm(i,j)));
    end
  end
end

fprintf('\n\n\n')
disp(ce)
z = [p1 v1 a1 s_a_lb s_a_ub];
for i = 1:length(ce)
  for j = 1:length(z)
    blockce(i,j) = simplify(diff(ce(i), z(j)));
    if blockce(i,j) ~= 0
      fprintf('jac_block(%d,%d) = %s;\n',i,j,char(blockce(i,j)));
    end
  end
end

fprintf('\n\n\n')
disp(Lf)
for i = 1:length(z)
  for j = 1:length(z)
    dif1 = diff(Lf, z(i));
    blocklf(i,j) = simplify(diff(dif1, z(j)));
    if blocklf(i,j) ~= 0
      fprintf('hes_block(%d,%d) = %s;\n',i,j,char(blocklf(i,j)));
    end
  end
end

fprintf('\n\n\n')
disp(Lm)
for i = 1:length(z)
  for j = 1:length(z)
    dif1 = diff(Lm, z(i));
    blocklm(i,j) = simplify(diff(dif1, z(j)));
    if blocklm(i,j) ~= 0
      fprintf('hes_block(%d,%d) = %s;\n',i,j,char(blocklm(i,j)));
    end
  end
end

fprintf('\n\n\n')
disp(Le)
for i = 1:length(z)
  for j = 1:length(z)
    dif1 = diff(Le, z(i));
    blockle(i,j) = simplify(diff(dif1, z(j)));
    if blockle(i,j) ~= 0
      fprintf('hes_block(%d,%d) = %s;\n',i,j,char(blockle(i,j)));
    end
  end
end
