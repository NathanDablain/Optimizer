clear
clc

% Should be able to supply differential equations, break out states and inputs
% as well as what variables will be bounded, and then be able to build the jacobian
% and hessian from there.

% Make a new ordering that groups knot vars together so slack isnt afterwords but
% instead packed with the states and inputs of a knot

% This model is a 3D rocket with position, speed, heading angle and flight path angle
% as states and lift as an input. There is an upper bound on lift and a lower bound on height
% The goal is to hit a target position at a target time

pkg load symbolic

function [dx, states, vars] = get_dx1(pn1, pe1, pu1, s1, psi1, gam1, Clp1, Clg1)
  syms c_D(s1) rho(pu1) g(pu1) A T1 m1

  Drag = 0.5*c_D(s1)*rho(pu1)*A*s1*s1;

  dx(1) = s1 * cos(gam1) * cos(psi1);
  dx(2) = s1 * cos(gam1) * sin(psi1);
  dx(3) = s1 * sin(gam1);
  dx(4) = ((T1 - Drag)/m1) - (g(pu1)*sin(gam1));
  dx(5) = (rho(pu1) * A * s1 * Clp1) / (m1 * 2 * cos(gam1));
  dx(6) = ((0.5 * rho(pu1) * A * s1 * Clg1)/m1) - ((g(pu1) * cos(gam1)) / s1);

  states = [pn1 pe1 pu1 s1 psi1 gam1];
  inputs = [Clp1 Clg1];
  vars = [states inputs];
end

function [dx, states, vars] = get_dx2(pn2, pe2, pu2, s2, psi2, gam2, Clp2, Clg2)
  syms c_D(s2) rho(pu2) g(pu2) A T2 m2

  Drag = 0.5*c_D(s2)*rho(pu2)*A*s2*s2;

  dx(1) = s2 * cos(gam2) * cos(psi2);
  dx(2) = s2 * cos(gam2) * sin(psi2);
  dx(3) = s2 * sin(gam2);
  dx(4) = ((T2 - Drag)/m2) - (g(pu2)*sin(gam2));
  dx(5) = (rho(pu2) * A * s2 * Clp2) / (m2 * 2 * cos(gam2));
  dx(6) = ((0.5 * rho(pu2) * A * s2 * Clg2)/m2) - ((g(pu2) * cos(gam2)) / s2);

  states = [pn2 pe2 pu2 s2 psi2 gam2];
  inputs = [Clp2 Clg2];
  vars = [states inputs];
end

function c = get_constraints_first_knot(states, inputs, slack, lbs, ubs, ic)
  syms h;

  pn1 = states(1,1);
  pe1 = states(1,2);
  pu1 = states(1,3);
  s1 = states(1,4);
  psi1 = states(1,5);
  gam1 = states(1,6);
  Clp1 = inputs(1,1);
  Clg1 = inputs(1,2);
  s_pu_lb = slack(1);
  s_Clp_lb = slack(2);
  s_Clg_lb = slack(3);
  s_Clp_ub = slack(4);
  s_Clg_ub = slack(5);
  lb_pu = lbs(1);
  lb_Clp = lbs(2);
  lb_Clg = lbs(3);
  ub_Clp = ubs(1);
  ub_Clg = ubs(2);
  pn_ic = ic(1);
  pe_ic = ic(2);
  pu_ic = ic(3);
  s_ic = ic(4);
  psi_ic = ic(5);
  gam_ic = ic(6);

  pn2 = states(2,1);
  pe2 = states(2,2);
  pu2 = states(2,3);
  s2 = states(2,4);
  psi2 = states(2,5);
  gam2 = states(2,6);
  Clp2 = inputs(2,1);
  Clg2 = inputs(2,2);

  dx1 = get_dx1(pn1, pe1, pu1, s1, psi1, gam1, Clp1, Clg1);
  dx2 = get_dx2(pn2, pe2, pu2, s2, psi2, gam2, Clp2, Clg2);
  c = [...
    % the trapezoidal defect constraints for the knot
    pn2 - pn1 - 0.5*h*(dx1(1) + dx2(1));...
    pe2 - pe1 - 0.5*h*(dx1(2) + dx2(2));...
    pu2 - pu1 - 0.5*h*(dx1(3) + dx2(3));...
    s2  - s1  - 0.5*h*(dx1(4) + dx2(4));...
    psi2-psi1 - 0.5*h*(dx1(5) + dx2(5));...
    gam2-gam1 - 0.5*h*(dx1(6) + dx2(6));...
    % the slack variable constraints for the knot
    % first lower bound
    lb_pu - pu1 + s_pu_lb;...
    lb_Clp - Clp1 + s_Clp_lb;...
    lb_Clg - Clg1 + s_Clg_lb;...
    % then upper bound
    Clp1 - ub_Clp + s_Clp_ub;...
    Clg1 - ub_Clg + s_Clg_ub;...
    % the initial condition constraints for the knot
    pn_ic - pn1;...
    pe_ic - pe1;...
    pu_ic - pu1;...
    s_ic - s1;...
    psi_ic - psi1;...
    gam_ic - gam1];
end

function c = get_constraints_middle_knot(states, inputs, slack, lbs, ubs)
  syms h;

  pn1 = states(1,1);
  pe1 = states(1,2);
  pu1 = states(1,3);
  s1 = states(1,4);
  psi1 = states(1,5);
  gam1 = states(1,6);
  Clp1 = inputs(1,1);
  Clg1 = inputs(1,2);
  s_pu_lb = slack(1);
  s_Clp_lb = slack(2);
  s_Clg_lb = slack(3);
  s_Clp_ub = slack(4);
  s_Clg_ub = slack(5);
  lb_pu = lbs(1);
  lb_Clp = lbs(2);
  lb_Clg = lbs(3);
  ub_Clp = ubs(1);
  ub_Clg = ubs(2);

  pn2 = states(2,1);
  pe2 = states(2,2);
  pu2 = states(2,3);
  s2 = states(2,4);
  psi2 = states(2,5);
  gam2 = states(2,6);
  Clp2 = inputs(2,1);
  Clg2 = inputs(2,2);

  dx1 = get_dx1(pn1, pe1, pu1, s1, psi1, gam1, Clp1, Clg1);
  dx2 = get_dx2(pn2, pe2, pu2, s2, psi2, gam2, Clp2, Clg2);
  c = [...
    % the trapezoidal defect constraints for the knot
    pn2 - pn1 - 0.5*h*(dx1(1) + dx2(1));...
    pe2 - pe1 - 0.5*h*(dx1(2) + dx2(2));...
    pu2 - pu1 - 0.5*h*(dx1(3) + dx2(3));...
    s2  - s1  - 0.5*h*(dx1(4) + dx2(4));...
    psi2-psi1 - 0.5*h*(dx1(5) + dx2(5));...
    gam2-gam1 - 0.5*h*(dx1(6) + dx2(6));...
    % the slack variable constraints for the knot
    % first lower bound
    lb_pu - pu1 + s_pu_lb;...
    lb_Clp - Clp1 + s_Clp_lb;...
    lb_Clg - Clg1 + s_Clg_lb;...
    % then upper bound
    Clp1 - ub_Clp + s_Clp_ub;...
    Clg1 - ub_Clg + s_Clg_ub];
end

function c = get_constraints_end_knot(states, inputs, slack, lbs, ubs)

  pn1 = states(1,1);
  pe1 = states(1,2);
  pu1 = states(1,3);
  s1 = states(1,4);
  psi1 = states(1,5);
  gam1 = states(1,6);
  Clp1 = inputs(1,1);
  Clg1 = inputs(1,2);
  s_pu_lb = slack(1);
  s_Clp_lb = slack(2);
  s_Clg_lb = slack(3);
  s_Clp_ub = slack(4);
  s_Clg_ub = slack(5);
  lb_pu = lbs(1);
  lb_Clp = lbs(2);
  lb_Clg = lbs(3);
  ub_Clp = ubs(1);
  ub_Clg = ubs(2);

  c = [...
    % the slack variable constraints for the knot
    % first lower bound
    lb_pu - pu1 + s_pu_lb;...
    lb_Clp - Clp1 + s_Clp_lb;...
    lb_Clg - Clg1 + s_Clg_lb;...
    % then upper bound
    Clp1 - ub_Clp + s_Clp_ub;...
    Clg1 - ub_Clg + s_Clg_ub];
end

function lagrangian = get_lagrangian_first_knot(states, inputs, slack, lbs, ubs, ic)
  syms mu lambda1 lambda2 lambda3 lambda4 lambda5 lambda6 lambda7...
          lambda8 lambda9 lambda10 lambda11 lambda12 lambda13 lambda14...
          lambda15 lambda16 lambda17 real;

  lambda = [lambda1 lambda2 lambda3 lambda4 lambda5 lambda6 lambda7...
            lambda8 lambda9 lambda10 lambda11 lambda12 lambda13 lambda14...
            lambda15 lambda16 lambda17];
  cost = -mu*(log(slack(1)) + log(slack(2)) + log(slack(3)) + log(slack(4)) + log(slack(5)));

  lagrangian = cost + lambda*get_constraints_first_knot(states, inputs, slack, lbs, ubs, ic);
end

function lagrangian = get_lagrangian_middle_knot(states, inputs, slack, lbs, ubs)
  syms mu lambda1 lambda2 lambda3 lambda4 lambda5 lambda6 lambda7 lambda8...
          lambda9 lambda10 lambda11 real;

  lambda = [lambda1 lambda2 lambda3 lambda4 lambda5 lambda6 lambda7 lambda8...
            lambda9 lambda10 lambda11];
  cost = -mu*(log(slack(1)) + log(slack(2)) + log(slack(3)) + log(slack(4)) + log(slack(5)));

  lagrangian = cost + lambda*get_constraints_middle_knot(states, inputs, slack, lbs, ubs);
end

function lagrangian = get_lagrangian_end_knot(states, inputs, slack, lbs, ubs, xd)
  syms mu lambda1 lambda2 lambda3 lambda4 lambda5 real;

  lambda = [lambda1 lambda2 lambda3 lambda4 lambda5];
  cost = (states(1,1) - xd(1))^2 + (states(1,2) - xd(2))^2 + (states(1,3) - xd(3))^2 -...
          mu*(log(slack(1)) + log(slack(2)) + log(slack(3)) + log(slack(4)) + log(slack(5)));

  lagrangian = cost + lambda*get_constraints_end_knot(states, inputs, slack, lbs, ubs);
end

syms pn1 pe1 pu1 s1 psi1 gam1 Clp1 Clg1 pn2 pe2 pu2 s2 psi2 gam2 Clp2 Clg2...
     s_pu_lb s_Clp_lb s_Clg_lb s_Clp_ub s_Clg_ub...
     lb_pu lb_Clp lb_Clg ub_Clp ub_Clg pn_ic pe_ic pu_ic s_ic psi_ic gam_ic xd1 xd2 xd3

states = [pn1 pe1 pu1 s1 psi1 gam1;...
          pn2 pe2 pu2 s2 psi2 gam2];
inputs = [Clp1 Clg1;...
          Clp2 Clg2];
slack = [s_pu_lb s_Clp_lb s_Clg_lb s_Clp_ub s_Clg_ub];
lbs = [lb_pu lb_Clp lb_Clg];
ubs = [ub_Clp ub_Clg];
ic = [pn_ic pe_ic pu_ic s_ic psi_ic gam_ic];
xd = [xd1 xd2 xd3];

z = [pn1 pe1 pu1 s1 psi1 gam1 Clp1 Clg1 s_pu_lb s_Clp_lb s_Clg_lb s_Clp_ub s_Clg_ub pn2 pe2 pu2 s2 psi2 gam2 Clp2 Clg2];
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
##disp(cf)
##for i = 1:length(cf)
##  for j = 1:length(z)
##    blockcf(i,j) = simplify(diff(cf(i), z(j)));
##    if blockcf(i,j) ~= 0
##      fprintf('jac_block(%d,%d) = %s;\n',i,j,char(blockcf(i,j)));
##    end
##  end
##end
##
##fprintf('\n\n\n')
##disp(cm)
##for i = 1:length(cm)
##  for j = 1:length(z)
##    blockcm(i,j) = simplify(diff(cm(i), z(j)));
##    if blockcm(i,j) ~= 0
##      fprintf('jac_block(%d,%d) = %s;\n',i,j,char(blockcm(i,j)));
##    end
##  end
##end
##
##fprintf('\n\n\n')
##disp(ce)
##z = [pn1 pe1 pu1 s1 psi1 gam1 Clp1 Clg1 s_pu_lb s_Clp_lb s_Clg_lb s_Clp_ub s_Clg_ub];
##for i = 1:length(ce)
##  for j = 1:length(z)
##    blockce(i,j) = simplify(diff(ce(i), z(j)));
##    if blockce(i,j) ~= 0
##      fprintf('jac_block(%d,%d) = %s;\n',i,j,char(blockce(i,j)));
##    end
##  end
##end

z = [pn1 pe1 pu1 s1 psi1 gam1 Clp1 Clg1 s_pu_lb s_Clp_lb s_Clg_lb s_Clp_ub s_Clg_ub];
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
