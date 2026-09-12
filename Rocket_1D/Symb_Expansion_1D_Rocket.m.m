clear
clc

% Should be able to supply differential equations, break out states and inputs
% as well as what variables will be bounded, and then be able to build the jacobian
% and hessian from there.

% Make a new ordering that groups knot vars together so slack isnt afterwords but
% instead packed with the states and inputs of a knot

% This model is a 1D rocket with position, velocity, and mass as states and thrust
% as an input. There is an upper and lower bound on thrust and a lower bound on mass
% The goal is to hit a target altitude at a target time

pkg load symbolic

function [dx, states, vars] = get_dx1(p1, v1, m1, T1)
  syms c_D(v1) rho(p1) g(p1) A C

  Drag = 0.5*c_D(v1)*rho(p1)*A*v1*v1;

  dx(1) = v1;
  dx(2) = (T1 - Drag)/m1 - g(p1);
  dx(3) = -T1 / C;

  states = [p1 v1 m1];
  inputs = T1;
  vars = [states inputs];
end

function [dx, states, vars] = get_dx2(p2, v2, m2, T2)
  syms c_D(v2) rho(p2) g(p2) A C

  Drag = 0.5*c_D(v2)*rho(p2)*A*v2*v2;

  dx(1) = v2;
  dx(2) = (T2 - Drag)/m2 - g(p2);
  dx(3) = -T2 / C;

  states = [p2 v2 m2];
  inputs = T2;
  vars = [states inputs];
end

function c = get_constraints_first_knot(states, inputs, slack, lbs, ubs, ic)
  syms h;

  p1 = states(1,1);
  v1 = states(1,2);
  m1 = states(1,3);
  T1 = inputs(1,1);
  s_m_lb = slack(1);
  s_T_lb = slack(2);
  s_T_ub = slack(3);
  lb_m = lbs(1);
  lb_T = lbs(2);
  ub_T = ubs(1);
  p_ic = ic(1);
  v_ic = ic(2);
  m_ic = ic(3);

  p2 = states(2,1);
  v2 = states(2,2);
  m2 = states(2,3);
  T2 = inputs(2,1);

  dx1 = get_dx1(p1, v1, m1, T1);
  dx2 = get_dx2(p2, v2, m2, T2);
  c = [...
    % the trapezoidal defect constraints for the knot
    p2 - p1 - 0.5*h*(dx1(1) + dx2(1));...
    v2 - v1 - 0.5*h*(dx1(2) + dx2(2));...
    m2 - m1 - 0.5*h*(dx1(3) + dx2(3));...
    % the slack variable constraints for the knot
    % first lower bound
    lb_m - m1 + s_m_lb;...
    lb_T - T1 + s_T_lb;...
    % then upper bound
    T1 - ub_T + s_T_ub;...
    % the initial condition constraints for the knot
    p_ic - p1;...
    v_ic - v1;...
    m_ic - m1];
end

function c = get_constraints_middle_knot(states, inputs, slack, lbs, ubs)
  syms h;

  p1 = states(1,1);
  v1 = states(1,2);
  m1 = states(1,3);
  T1 = inputs(1,1);
  s_m_lb = slack(1);
  s_T_lb = slack(2);
  s_T_ub = slack(3);
  lb_m = lbs(1);
  lb_T = lbs(2);
  ub_T = ubs(1);

  p2 = states(2,1);
  v2 = states(2,2);
  m2 = states(2,3);
  T2 = inputs(2,1);

  dx1 = get_dx1(p1, v1, m1, T1);
  dx2 = get_dx2(p2, v2, m2, T2);
  c = [...
    % the trapezoidal defect constraints for the knot
    p2 - p1 - 0.5*h*(dx1(1) + dx2(1));...
    v2 - v1 - 0.5*h*(dx1(2) + dx2(2));...
    m2 - m1 - 0.5*h*(dx1(3) + dx2(3));...
    % the slack variable constraints for the knot
    % first lower bound
    lb_m - m1 + s_m_lb;...
    lb_T - T1 + s_T_lb;...
    % then upper bound
    T1 - ub_T + s_T_ub];
end

function c = get_constraints_end_knot(states, inputs, slack, lbs, ubs)
  syms h;

  p1 = states(1,1);
  v1 = states(1,2);
  m1 = states(1,3);
  T1 = inputs(1,1);
  s_m_lb = slack(1);
  s_T_lb = slack(2);
  s_T_ub = slack(3);
  lb_m = lbs(1);
  lb_T = lbs(2);
  ub_T = ubs(1);

  c = [...
    % the slack variable constraints for the knot
    % first lower bound
    lb_m - m1 + s_m_lb;...
    lb_T - T1 + s_T_lb;...
    % then upper bound
    T1 - ub_T + s_T_ub];
end

function lagrangian = get_lagrangian_first_knot(states, inputs, slack, lbs, ubs, ic)
  syms mu lambda1 lambda2 lambda3 lambda4 lambda5 lambda6 lambda7 lambda8 lambda9 real;

  lambda = [lambda1 lambda2 lambda3 lambda4 lambda5 lambda6 lambda7 lambda8 lambda9];
  cost = -mu*(log(slack(1)) + log(slack(2)) + log(slack(3)));

  lagrangian = cost + lambda*get_constraints_first_knot(states, inputs, slack, lbs, ubs, ic);
end

function lagrangian = get_lagrangian_middle_knot(states, inputs, slack, lbs, ubs)
  syms mu lambda1 lambda2 lambda3 lambda4 lambda5 lambda6 real;

  lambda = [lambda1 lambda2 lambda3 lambda4 lambda5 lambda6];
  cost = -mu*(log(slack(1)) + log(slack(2)) + log(slack(3)));

  lagrangian = cost + lambda*get_constraints_middle_knot(states, inputs, slack, lbs, ubs);
end

function lagrangian = get_lagrangian_end_knot(states, inputs, slack, lbs, ubs, xd)
  syms mu lambda1 lambda2 lambda3 real;

  lambda = [lambda1 lambda2 lambda3];
  cost = (states(1,1) - xd(1))^2 - mu*(log(slack(1)) + log(slack(2)) + log(slack(3)));

  lagrangian = cost + lambda*get_constraints_end_knot(states, inputs, slack, lbs, ubs);
end

syms p1 v1 m1 T1 p2 v2 m2 T2 s_m_lb s_T_lb s_T_ub lb_m lb_T ub_T p_ic v_ic m_ic xd1 xd2

states = [p1 v1 m1;...
          p2 v2 m2];
inputs = [T1;...
          T2];
slack = [s_m_lb, s_T_lb, s_T_ub];
lbs = [lb_m, lb_T];
ubs = ub_T;
ic = [p_ic, v_ic m_ic];
xd = xd1;

z = [p1 v1 m1 T1 s_m_lb s_T_lb s_T_ub p2 v2 m2 T2];
% All knots have 3 states, 1 input, 3 slack variables -> 7 decision variables
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
z = [p1 v1 m1 T1 s_m_lb s_T_lb s_T_ub];
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
      if i > j
        fprintf('hes_block(%d,%d) = hes_block(%d,%d);\n',i,j,j,i);
      else
        fprintf('hes_block(%d,%d) = %s;\n',i,j,char(blocklf(i,j)));
      end
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
      if i > j
        fprintf('hes_block(%d,%d) = hes_block(%d,%d);\n',i,j,j,i);
      else
        fprintf('hes_block(%d,%d) = %s;\n',i,j,char(blocklm(i,j)));
      end
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
      if i > j
        fprintf('hes_block(%d,%d) = hes_block(%d,%d);\n',i,j,j,i);
      else
        fprintf('hes_block(%d,%d) = %s;\n',i,j,char(blockle(i,j)));
      end
    end
  end
end
