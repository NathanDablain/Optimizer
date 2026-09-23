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

function [dx, states, vars] = get_dx1(pn1, pe1, pd1, s1, q01, q11, q21, q31, Cly1, Clz1)
  syms c_D(s1) rho(pd1) g(pd1) A T1 m1

  dx(1) = s1*(-2*q21^2 - 2*q31^2 + 1);
  dx(2) = 2*s1*(q01*q31 + q11*q21);
  dx(3) = 2*s1*(-q01*q21 + q11*q31);
  dx(4) = (-A*c_D(s1)*rho(pd1)*s1^2/2 + T1 - 2*g(pd1)*m1*(q01*q21 - q11*q31))/m1;
  dx(5) = (q21*(A*Clz1*rho(pd1)*s1^2 - 2*g(pd1)*m1*(2*q11^2 + 2*q21^2 - 1)) - q31*(A*Cly1*rho(pd1)*s1^2 + 4*g(pd1)*m1*(q01*q11 + q21*q31)))/(4*m1*s1);
  dx(6) = -A*Cly1*q21*rho(pd1)*s1/(4*m1) - A*Clz1*q31*rho(pd1)*s1/(4*m1) - g(pd1)*q01*q11*q21/s1 + g(pd1)*q11^2*q31/s1 - g(pd1)*q31/(2*s1);
  dx(7) = (-q01*(A*Clz1*rho(pd1)*s1^2 - 2*g(pd1)*m1*(2*q11^2 + 2*q21^2 - 1)) + q11*(A*Cly1*rho(pd1)*s1^2 + 4*g(pd1)*m1*(q01*q11 + q21*q31)))/(4*m1*s1);
  dx(8) = (q01*(A*Cly1*rho(pd1)*s1^2 + 4*g(pd1)*m1*(q01*q11 + q21*q31)) + q11*(A*Clz1*rho(pd1)*s1^2 - 2*g(pd1)*m1*(2*q11^2 + 2*q21^2 - 1)))/(4*m1*s1);

  states = [pn1, pe1, pd1, s1, q01, q11, q21, q31];
  inputs = [Cly1, Clz1];
  vars = [states inputs];
end

function [dx, states, vars] = get_dx2(pn2, pe2, pd2, s2, q02, q12, q22, q32, Cly2, Clz2)
  syms c_D(s2) rho(pd2) g(pd2) A T2 m2

  dx(1) = s2*(-2*q22^2 - 2*q32^2 + 1);
  dx(2) = 2*s2*(q02*q32 + q12*q22);
  dx(3) = 2*s2*(-q02*q22 + q12*q32);
  dx(4) = (-A*c_D(s2)*rho(pd2)*s2^2/2 + T2 - 2*g(pd2)*m2*(q02*q22 - q12*q32))/m2;
  dx(5) = (q22*(A*Clz2*rho(pd2)*s2^2 - 2*g(pd2)*m2*(2*q12^2 + 2*q22^2 - 1)) - q32*(A*Cly2*rho(pd2)*s2^2 + 4*g(pd2)*m2*(q02*q12 + q22*q32)))/(4*m2*s2);
  dx(6) = -A*Cly2*q22*rho(pd2)*s2/(4*m2) - A*Clz2*q32*rho(pd2)*s2/(4*m2) - g(pd2)*q02*q12*q22/s2 + g(pd2)*q12^2*q32/s2 - g(pd2)*q32/(2*s2);
  dx(7) = (-q02*(A*Clz2*rho(pd2)*s2^2 - 2*g(pd2)*m2*(2*q12^2 + 2*q22^2 - 1)) + q12*(A*Cly2*rho(pd2)*s2^2 + 4*g(pd2)*m2*(q02*q12 + q22*q32)))/(4*m2*s2);
  dx(8) = (q02*(A*Cly2*rho(pd2)*s2^2 + 4*g(pd2)*m2*(q02*q12 + q22*q32)) + q12*(A*Clz2*rho(pd2)*s2^2 - 2*g(pd2)*m2*(2*q12^2 + 2*q22^2 - 1)))/(4*m2*s2);

  states = [pn2, pe2, pd2, s2, q02, q12, q22, q32];
  inputs = [Cly2, Clz2];
  vars = [states inputs];
end

function c = get_constraints_first_knot(states, inputs, slack, lbs, ubs, ic)
  syms h;

  pn1 = states(1,1);
  pe1 = states(1,2);
  pd1 = states(1,3);
  s1 = states(1,4);
  q01 = states(1,5);
  q11 = states(1,6);
  q21 = states(1,7);
  q31 = states(1,8);
  Cly1 = inputs(1,1);
  Clz1 = inputs(1,2);
  s_Cly_lb = slack(1);
  s_Clz_lb = slack(2);
  s_pd_ub = slack(3);
  s_Cly_ub = slack(4);
  s_Clz_ub = slack(5);
  lb_Cly = lbs(1);
  lb_Clz = lbs(2);
  ub_pd = ubs(1);
  ub_Cly = ubs(2);
  ub_Clz = ubs(3);
  pn_ic = ic(1);
  pe_ic = ic(2);
  pd_ic = ic(3);
  s_ic = ic(4);
  q0_ic = ic(5);
  q1_ic = ic(6);
  q2_ic = ic(7);
  q3_ic = ic(8);

  pn2 = states(2,1);
  pe2 = states(2,2);
  pd2 = states(2,3);
  s2 = states(2,4);
  q02 = states(2,5);
  q12 = states(2,6);
  q22 = states(2,7);
  q32 = states(2,8);
  Cly2 = inputs(2,1);
  Clz2 = inputs(2,2);

  dx1 = get_dx1(pn1, pe1, pd1, s1, q01, q11, q21, q31, Cly1, Clz1);
  dx2 = get_dx2(pn2, pe2, pd2, s2, q02, q12, q22, q32, Cly2, Clz2);
  c = [...
    % the trapezoidal defect constraints for the knot
    pn2 - pn1 - 0.5*h*(dx1(1) + dx2(1));...
    pe2 - pe1 - 0.5*h*(dx1(2) + dx2(2));...
    pd2 - pd1 - 0.5*h*(dx1(3) + dx2(3));...
    s2  - s1  - 0.5*h*(dx1(4) + dx2(4));...
    q02 - q01 - 0.5*h*(dx1(5) + dx2(5));...
    q12 - q11 - 0.5*h*(dx1(6) + dx2(6));...
    q22 - q21 - 0.5*h*(dx1(7) + dx2(7));...
    q32 - q31 - 0.5*h*(dx1(8) + dx2(8));...
    % the unit length constraint on the quaternion
    q01^2 + q11^2 + q21^2 + q31^2 - 1;...
    % the slack variable constraints for the knot
    % first lower bound
    lb_Cly - Cly1 + s_Cly_lb;...
    lb_Clz - Clz1 + s_Clz_lb;...
    % then upper bound
    pd1 - ub_pd + s_pd_ub;...
    Cly1 - ub_Cly + s_Cly_ub;...
    Clz1 - ub_Clz + s_Clz_ub;...
    % the initial condition constraints for the knot
    pn_ic - pn1;...
    pe_ic - pe1;...
    pd_ic - pd1;...
    s_ic - s1;...
    q0_ic - q01;...
    q1_ic - q11;...
    q2_ic - q21;...
    q3_ic - q31];
end

function c = get_constraints_middle_knot(states, inputs, slack, lbs, ubs)
  syms h;

  pn1 = states(1,1);
  pe1 = states(1,2);
  pd1 = states(1,3);
  s1 = states(1,4);
  q01 = states(1,5);
  q11 = states(1,6);
  q21 = states(1,7);
  q31 = states(1,8);
  Cly1 = inputs(1,1);
  Clz1 = inputs(1,2);
  s_Cly_lb = slack(1);
  s_Clz_lb = slack(2);
  s_pd_ub = slack(3);
  s_Cly_ub = slack(4);
  s_Clz_ub = slack(5);
  lb_Cly = lbs(1);
  lb_Clz = lbs(2);
  ub_pd = ubs(1);
  ub_Cly = ubs(2);
  ub_Clz = ubs(3);

  pn2 = states(2,1);
  pe2 = states(2,2);
  pd2 = states(2,3);
  s2 = states(2,4);
  q02 = states(2,5);
  q12 = states(2,6);
  q22 = states(2,7);
  q32 = states(2,8);
  Cly2 = inputs(2,1);
  Clz2 = inputs(2,2);

  dx1 = get_dx1(pn1, pe1, pd1, s1, q01, q11, q21, q31, Cly1, Clz1);
  dx2 = get_dx2(pn2, pe2, pd2, s2, q02, q12, q22, q32, Cly2, Clz2);
  c = [...
    % the trapezoidal defect constraints for the knot
    pn2 - pn1 - 0.5*h*(dx1(1) + dx2(1));...
    pe2 - pe1 - 0.5*h*(dx1(2) + dx2(2));...
    pd2 - pd1 - 0.5*h*(dx1(3) + dx2(3));...
    s2  - s1  - 0.5*h*(dx1(4) + dx2(4));...
    q02 - q01 - 0.5*h*(dx1(5) + dx2(5));...
    q12 - q11 - 0.5*h*(dx1(6) + dx2(6));...
    q22 - q21 - 0.5*h*(dx1(7) + dx2(7));...
    q32 - q31 - 0.5*h*(dx1(8) + dx2(8));...
    % the unit length constraint on the quaternion
    q01^2 + q11^2 + q21^2 + q31^2 - 1;...
    % the slack variable constraints for the knot
    % first lower bound
    lb_Cly - Cly1 + s_Cly_lb;...
    lb_Clz - Clz1 + s_Clz_lb;...
    % then upper bound
    pd1 - ub_pd + s_pd_ub;...
    Cly1 - ub_Cly + s_Cly_ub;...
    Clz1 - ub_Clz + s_Clz_ub];
end

function c = get_constraints_end_knot(states, inputs, slack, lbs, ubs)

  pn1 = states(1,1);
  pe1 = states(1,2);
  pd1 = states(1,3);
  s1 = states(1,4);
  q01 = states(1,5);
  q11 = states(1,6);
  q21 = states(1,7);
  q31 = states(1,8);
  Cly1 = inputs(1,1);
  Clz1 = inputs(1,2);
  s_Cly_lb = slack(1);
  s_Clz_lb = slack(2);
  s_pd_ub = slack(3);
  s_Cly_ub = slack(4);
  s_Clz_ub = slack(5);
  lb_Cly = lbs(1);
  lb_Clz = lbs(2);
  ub_pd = ubs(1);
  ub_Cly = ubs(2);
  ub_Clz = ubs(3);

  c = [...
    % the unit length constraint on the quaternion
    q01^2 + q11^2 + q21^2 + q31^2 - 1;...
    % the slack variable constraints for the knot
    % first lower bound
    lb_Cly - Cly1 + s_Cly_lb;...
    lb_Clz - Clz1 + s_Clz_lb;...
    % then upper bound
    pd1 - ub_pd + s_pd_ub;...
    Cly1 - ub_Cly + s_Cly_ub;...
    Clz1 - ub_Clz + s_Clz_ub];
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

syms pn1 pe1 pd1 s1 q01 q11 q21 q31 Cly1 Clz1...
     pn2 pe2 pd2 s2 q02 q12 q22 q32 Cly2 Clz2...
     s_Cly_lb s_Clz_lb s_pd_ub s_Cly_ub s_Clz_ub...
     lb_Cly lb_Clz ub_pd ub_Cly ub_Clz...
     pn_ic pe_ic pd_ic s_ic q0_ic q1_ic q2_ic q3_ic xd1 xd2 xd3

states = [pn1 pe1 pd1 s1 q01 q11 q21 q31;...
          pn2 pe2 pd2 s2 q02 q12 q22 q32];
inputs = [Cly1 Clz1;...
          Cly2 Clz2];
slack = [s_Cly_lb s_Clz_lb s_pd_ub s_Cly_ub s_Clz_ub];
lbs = [lb_Cly lb_Clz];
ubs = [ub_pd ub_Cly ub_Clz];
ic = [pn_ic pe_ic pd_ic s_ic q0_ic q1_ic q2_ic q3_ic];
xd = [xd1 xd2 xd3];

z = [pn1 pe1 pd1 s1 q01 q11 q21 q31 Cly1 Clz1 s_Cly_lb s_Clz_lb s_pd_ub s_Cly_ub s_Clz_ub pn2 pe2 pd2 s2 q02 q12 q22 q32 Cly2 Clz2];

cf = get_constraints_first_knot(states, inputs, slack, lbs, ubs, ic);

cm = get_constraints_middle_knot(states, inputs, slack, lbs, ubs);

ce = get_constraints_end_knot(states, inputs, slack, lbs, ubs);

##Lf = get_lagrangian_first_knot(states, inputs, slack, lbs, ubs, ic);
##Lm = get_lagrangian_middle_knot(states, inputs, slack, lbs, ubs);
##Le = get_lagrangian_end_knot(states, inputs, slack, lbs, ubs, xd);

clc
disp(cf)
for i = 1:length(cf)
  for j = 1:length(z)
    blockcf = simplify(diff(cf(i), z(j)));
    if blockcf ~= 0
      fprintf('jac_block(%d,%d) = %s;\n',i,j,char(blockcf));
    end
  end
end

fprintf('\n\n\n')
disp(cm)
for i = 1:length(cm)
  for j = 1:length(z)
    blockcm = simplify(diff(cm(i), z(j)));
    if blockcm ~= 0
      fprintf('jac_block(%d,%d) = %s;\n',i,j,char(blockcm));
    end
  end
end

fprintf('\n\n\n')
disp(ce)
z = [pn1 pe1 pd1 s1 q01 q11 q21 q31 Cly1 Clz1 s_Cly_lb s_Clz_lb s_pd_ub s_Cly_ub s_Clz_ub];
for i = 1:length(ce)
  for j = 1:length(z)
    blockce = simplify(diff(ce(i), z(j)));
    if blockce ~= 0
      fprintf('jac_block(%d,%d) = %s;\n',i,j,char(blockce));
    end
  end
end

##z = [pn1 pe1 pu1 s1 psi1 gam1 Clp1 Clg1 s_pu_lb s_Clp_lb s_Clg_lb s_Clp_ub s_Clg_ub];
##fprintf('\n\n\n')
##disp(Lf)
##for i = 1:length(z)
##  for j = 1:length(z)
##    dif1 = diff(Lf, z(i));
##    blocklf(i,j) = simplify(diff(dif1, z(j)));
##    if blocklf(i,j) ~= 0
##      if i > j
##        fprintf('hes_block(%d,%d) = hes_block(%d,%d);\n',i,j,j,i);
##      else
##        fprintf('hes_block(%d,%d) = %s;\n',i,j,char(blocklf(i,j)));
##      end
##    end
##  end
##end
##
##fprintf('\n\n\n')
##disp(Lm)
##for i = 1:length(z)
##  for j = 1:length(z)
##    dif1 = diff(Lm, z(i));
##    blocklm(i,j) = simplify(diff(dif1, z(j)));
##    if blocklm(i,j) ~= 0
##      if i > j
##        fprintf('hes_block(%d,%d) = hes_block(%d,%d);\n',i,j,j,i);
##      else
##        fprintf('hes_block(%d,%d) = %s;\n',i,j,char(blocklm(i,j)));
##      end
##    end
##  end
##end
##
##fprintf('\n\n\n')
##disp(Le)
##for i = 1:length(z)
##  for j = 1:length(z)
##    dif1 = diff(Le, z(i));
##    blockle(i,j) = simplify(diff(dif1, z(j)));
##    if blockle(i,j) ~= 0
##      if i > j
##        fprintf('hes_block(%d,%d) = hes_block(%d,%d);\n',i,j,j,i);
##      else
##        fprintf('hes_block(%d,%d) = %s;\n',i,j,char(blockle(i,j)));
##      end
##    end
##  end
##end
