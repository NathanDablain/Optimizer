clear
clc

% Should be able to supply differential equations, break out states and inputs
% as well as what variables will be bounded, and then be able to build the jacobian
% and hessian from there.

% Make a new ordering that groups knot vars together so slack isnt afterwords but
% instead packed with the states and inputs of a knot

% This model is a 3D rocket with NED position, speed, NED to velocity quaternion, and angular rates
% as states and perpendicular and vertical lift as an inputs. There is an upper bound on lift and an upper bound on
% the down position. The goal is to hit a target position at a target time

pkg load symbolic

function [dx, w] = get_dx1(pn1, pe1, pd1, s1, q01, q11, q21, q31, wy1, wz1, Cly1, Clz1)
  syms c_D(s1) rho(pd1) A g T1 m1

  q = [q01; q11; q21; q31];
  D = 0.5*rho(pd1)*A*c_D(s1)*s1*s1;

  v_NED = [s1*(-2*q21^2 - 2*q31^2 + 1); 2*s1*(q01*q31 + q11*q21); 2*s1*(-q01*q21 + q11*q31)];
  omega = [[0 0 -wy1 -wz1]; [0 0 wz1 -wy1]; [wy1 -wz1 0 0]; [wz1 wy1 0 0]];

  dx(1:3) = v_NED;
  dx(4) = ((T1 - D)/m1) + 2*g*(-q01*q21 + q11*q31);
  dx(5:8) = 0.5.*omega*q;

  w = [(-(0.5*rho(pd1)*A*Clz1*s1*s1/m1) - g*(-2*q11^2 - 2*q21^2 + 1))/s1;...
       ((0.5*rho(pd1)*A*Cly1*s1*s1/m1) + 2*g*(q01*q11 + q21*q31))/s1];
end

function [dx, w] = get_dx2(pn2, pe2, pd2, s2, q02, q12, q22, q32, wy2, wz2, Cly2, Clz2)
  syms c_D(s2) rho(pd2) A g T2 m2

  q = [q02; q12; q22; q32];
  D = 0.5*rho(pd2)*A*c_D(s2)*s2*s2;

  v_NED = [s2*(-2*q22^2 - 2*q32^2 + 1); 2*s2*(q02*q32 + q12*q22); 2*s2*(-q02*q22 + q12*q32)];
  omega = [[0 0 -wy2 -wz2]; [0 0 wz2 -wy2]; [wy2 -wz2 0 0]; [wz2 wy2 0 0]];

  dx(1:3) = v_NED;
  dx(4) = ((T2 - D)/m2) + 2*g*(-q02*q22 + q12*q32);
  dx(5:8) = 0.5.*omega*q;

  w = [(-(0.5*rho(pd2)*A*Clz2*s2*s2/m2) - g*(-2*q12^2 - 2*q22^2 + 1))/s2;...
       ((0.5*rho(pd2)*A*Cly2*s2*s2/m2) + 2*g*(q02*q12 + q22*q32))/s2];
end

function c = get_constraints_first_knot(states, inputs, slack, lbs, ubs, ic)
  syms h g;

  pn1 = states(1,1);
  pe1 = states(1,2);
  pd1 = states(1,3);
  s1 = states(1,4);
  q01 = states(1,5);
  q11 = states(1,6);
  q21 = states(1,7);
  q31 = states(1,8);
  wy1 = states(1,9);
  wz1 = states(1,10);
  Cly1 = inputs(1,1);
  Clz1 = inputs(1,2);
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
  wy2 = states(2,9);
  wz2 = states(2,10);
  Cly2 = inputs(2,1);
  Clz2 = inputs(2,2);

  [dx1, w_guess] = get_dx1(pn1, pe1, pd1, s1, q01, q11, q21, q31, wy1, wz1, Cly1, Clz1);
  dx1 = simplify(dx1);
  [dx2, w_guess2] = get_dx2(pn2, pe2, pd2, s2, q02, q12, q22, q32, wy2, wz2, Cly2, Clz2);
  dx2 = simplify(dx2);

##  wy = 0.5*(wy1 + wy2);
##  wz = 0.5*(wz1 + wz2);
  w = [wy1 wz1];

  theta = norm(w)*h;

  delta_q = [cos(theta/2);0;(w(1)/norm(w))*sin(theta/2);(w(2)/norm(w))*sin(theta/2)];

  constraint = [q02;q12;q22;q32] - [q01 -q11 -q21 -q31; q11 q01 -q31 q21; q21 q31 q01 -q11; q31 -q21 q11 q01]*delta_q;


  c = [...
    % the trapezoidal defect constraints for the knot
    pn2 - pn1 - 0.5*h*(dx1(1) + dx2(1));...
    pe2 - pe1 - 0.5*h*(dx1(2) + dx2(2));...
    pd2 - pd1 - 0.5*h*(dx1(3) + dx2(3));...
    s2  - s1  - 0.5*h*(dx1(4) + dx2(4));...
    constraint(1);...
    constraint(2);...
    constraint(3);...
    constraint(4);...
    % the angular rate defect constraints for the knot
    wy1 - 0.5*(w_guess(1) + w_guess2(1));...
    wz1 - 0.5*(w_guess(2) + w_guess2(2));...
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
  syms h g;

  pn1 = states(1,1);
  pe1 = states(1,2);
  pd1 = states(1,3);
  s1 = states(1,4);
  q01 = states(1,5);
  q11 = states(1,6);
  q21 = states(1,7);
  q31 = states(1,8);
  wy1 = states(1,9);
  wz1 = states(1,10);
  Cly1 = inputs(1,1);
  Clz1 = inputs(1,2);

  pn2 = states(2,1);
  pe2 = states(2,2);
  pd2 = states(2,3);
  s2 = states(2,4);
  q02 = states(2,5);
  q12 = states(2,6);
  q22 = states(2,7);
  q32 = states(2,8);
  wy2 = states(2,9);
  wz2 = states(2,10);
  Cly2 = inputs(2,1);
  Clz2 = inputs(2,2);

  [dx1, w_guess] = get_dx1(pn1, pe1, pd1, s1, q01, q11, q21, q31, wy1, wz1, Cly1, Clz1);
  dx1 = simplify(dx1);
  [dx2, w_guess2] = get_dx2(pn2, pe2, pd2, s2, q02, q12, q22, q32, wy2, wz2, Cly2, Clz2);
  dx2 = simplify(dx2);

##  wy = 0.5*(wy1 + wy2);
##  wz = 0.5*(wz1 + wz2);
  w = [wy1 wz1];

  theta = norm(w)*h;

  delta_q = [cos(theta/2);0;(w(1)/norm(w))*sin(theta/2);(w(2)/norm(w))*sin(theta/2)];

  constraint = [q02;q12;q22;q32] - [q01 -q11 -q21 -q31; q11 q01 -q31 q21; q21 q31 q01 -q11; q31 -q21 q11 q01]*delta_q;

  c = [...
    % the trapezoidal defect constraints for the knot
    pn2 - pn1 - 0.5*h*(dx1(1) + dx2(1));...
    pe2 - pe1 - 0.5*h*(dx1(2) + dx2(2));...
    pd2 - pd1 - 0.5*h*(dx1(3) + dx2(3));...
    s2  - s1  - 0.5*h*(dx1(4) + dx2(4));...
    constraint(1);...
    constraint(2);...
    constraint(3);...
    constraint(4);...
    % the angular rate defect constraints for the knot
    wy1 - 0.5*h*(w_guess(1) + w_guess2(1));...
    wz1 - 0.5*h*(w_guess(2) + w_guess2(2))];
end

function c = get_constraints_end_knot(states, inputs, slack, lbs, ubs)
  syms h g;

  pn1 = states(1,1);
  pe1 = states(1,2);
  pd1 = states(1,3);
  s1 = states(1,4);
  q01 = states(1,5);
  q11 = states(1,6);
  q21 = states(1,7);
  q31 = states(1,8);
  wy1 = states(1,9);
  wz1 = states(1,10);
  Cly1 = inputs(1,1);
  Clz1 = inputs(1,2);

  c = [];
end

function j = get_trapezoidal_integration(states, inputs)
  syms h g;

  pn1 = states(1,1);
  pe1 = states(1,2);
  pd1 = states(1,3);
  s1 = states(1,4);
  q01 = states(1,5);
  q11 = states(1,6);
  q21 = states(1,7);
  q31 = states(1,8);
  wy1 = states(1,9);
  wz1 = states(1,10);
  Cly1 = inputs(1,1);
  Clz1 = inputs(1,2);

  pn2 = states(2,1);
  pe2 = states(2,2);
  pd2 = states(2,3);
  s2 = states(2,4);
  q02 = states(2,5);
  q12 = states(2,6);
  q22 = states(2,7);
  q32 = states(2,8);
  wy2 = states(2,9);
  wz2 = states(2,10);
  Cly2 = inputs(2,1);
  Clz2 = inputs(2,2);

  [dx1, w_guess] = get_dx1(pn1, pe1, pd1, s1, q01, q11, q21, q31, wy1, wz1, Cly1, Clz1);
  dx1 = simplify(dx1);
  [dx2, w_guess2] = get_dx2(pn2, pe2, pd2, s2, q02, q12, q22, q32, wy2, wz2, Cly2, Clz2);
  dx2 = simplify(dx2);

  w = [wy1 wz1];

  theta = norm(w)*h;

  delta_q = [cos(theta/2);0;(w(1)/norm(w))*sin(theta/2);(w(2)/norm(w))*sin(theta/2)];

  constraint = [q02;q12;q22;q32] - [q01 -q11 -q21 -q31; q11 q01 -q31 q21; q21 q31 q01 -q11; q31 -q21 q11 q01]*delta_q;

  j = [...
    % the trapezoidal defect constraints for the knot
    pn2 - pn1 - 0.5*h*(dx1(1) + dx2(1));...
    pe2 - pe1 - 0.5*h*(dx1(2) + dx2(2));...
    pd2 - pd1 - 0.5*h*(dx1(3) + dx2(3));...
    s2  - s1  - 0.5*h*(dx1(4) + dx2(4));...
    constraint(1);...
    constraint(2);...
    constraint(3);...
    constraint(4)];
end

function lagrangian = get_lagrangian_first_knot(states, inputs, slack, lbs, ubs, ic)
  syms mu lambda1 lambda2 lambda3 lambda4 lambda5 lambda6 lambda7...
          lambda8 lambda9 lambda10 lambda11 lambda12 lambda13 lambda14...
          lambda15 lambda16 lambda17 lambda18 real;

  lambda = [lambda1 lambda2 lambda3 lambda4 lambda5 lambda6 lambda7...
            lambda8 lambda9 lambda10 lambda11 lambda12 lambda13 lambda14...
            lambda15 lambda16 lambda17 lambda18];

  cost = -mu*(log(inputs(1,1) - lbs(1)) +...
                  log(inputs(1,2) - lbs(2)) +...
                  log(ubs(2) - inputs(1,1)) +...
                  log(ubs(3) - inputs(1,2)));

  cf = get_constraints_first_knot(states, inputs, slack, lbs, ubs, ic);
  lagrangian = cost + lambda*cf;
end

function lagrangian = get_lagrangian_middle_knot(states, inputs, slack, lbs, ubs)
  syms mu mu_q lambda1 lambda2 lambda3 lambda4 lambda5 lambda6 lambda7...
          lambda8 lambda9 lambda10 real;

  lambda = [lambda1 lambda2 lambda3 lambda4 lambda5 lambda6 lambda7...
            lambda8 lambda9 lambda10];
##  cost = -mu*(log(slack(1)) + log(slack(2)) + log(slack(3)) + log(slack(4)) + log(slack(5)));
  cost = -mu*(log(inputs(1,1) - lbs(1)) +...
                  log(inputs(1,2) - lbs(2)) +...
                  log(ubs(2) - inputs(1,1)) +...
                  log(ubs(3) - inputs(1,2))) +...
          mu_q*(states(1,5)^2 + states(1,6)^2 + states(1,7)^2 + states(1,8)^2 - 1);
  cm = get_constraints_middle_knot(states, inputs, slack, lbs, ubs);
  lagrangian = cost + lambda*cm;
end

syms pn1 pe1 pd1 s1 q01 q11 q21 q31 wy1 wz1 Cly1 Clz1...
     pn2 pe2 pd2 s2 q02 q12 q22 q32 wy2 wz2 Cly2 Clz2...
     s_Cly_lb s_Clz_lb s_pd_ub s_Cly_ub s_Clz_ub...
     lb_Cly lb_Clz ub_pd ub_Cly ub_Clz...
     pn_ic pe_ic pd_ic s_ic q0_ic q1_ic q2_ic q3_ic xd1 xd2 xd3 real

states = [pn1 pe1 pd1 s1 q01 q11 q21 q31 wy1 wz1;...
          pn2 pe2 pd2 s2 q02 q12 q22 q32 wy2 wz2];
inputs = [Cly1 Clz1;...
          Cly2 Clz2];
slack = [s_Cly_lb s_Clz_lb s_pd_ub s_Cly_ub s_Clz_ub];
lbs = [lb_Cly lb_Clz];
ubs = [ub_pd ub_Cly ub_Clz];
ic = [pn_ic pe_ic pd_ic s_ic q0_ic q1_ic q2_ic q3_ic];
xd = [xd1 xd2 xd3];

##z = [pn1 pe1 pd1 s1 q01 q11 q21 q31 wy1 wz1 Cly1 Clz1 pn2 pe2 pd2 s2 q02 q12 q22 q32 wy2 wz2 Cly2 Clz2];

##cf = simplify(get_constraints_first_knot(states, inputs, slack, lbs, ubs, ic));

##cm = simplify(get_constraints_middle_knot(states, inputs, slack, lbs, ubs));
##
##ce = simplify(get_constraints_end_knot(states, inputs, slack, lbs, ubs));

##Lf = get_lagrangian_first_knot(states, inputs, slack, lbs, ubs, ic);

##Lm = simplify(get_lagrangian_middle_knot(states, inputs, slack, lbs, ubs));

J = get_trapezoidal_integration(states, inputs);

clc
##disp(cf)
##for i = 1:length(cf)
##  for j = 1:length(z)
##    block_con = simplify(diff(cf(i), z(j)));
##    if block_con ~= 0
##      fprintf('jac_block(%d,%d) = %s;\n',i,j,char(block_con));
##    end
##  end
##end
##fprintf('\n\n\n')
##
##for i = 1:length(cm)
##  for j = 1:length(z)
##    block_con = simplify(diff(cm(i), z(j)));
##    if block_con ~= 0
##      fprintf('jac_block(%d,%d) = %s;\n',i,j,char(block_con));
##    end
##  end
##end
##fprintf('\n\n\n')
##
##z = [pn1 pe1 pd1 s1 q01 q11 q21 q31 wy1 wz1 Cly1 Clz1 s_Cly_lb s_Clz_lb s_pd_ub s_Cly_ub s_Clz_ub];
##for i = 1:length(ce)
##  for j = 1:length(z)
##    block_con = simplify(diff(ce(i), z(j)));
##    if block_con ~= 0
##      fprintf('jac_block(%d,%d) = %s;\n',i,j,char(block_con));
##    end
##  end
##end
##fprintf('\n\n\n')

##z = [pn1 pe1 pd1 s1 q01 q11 q21 q31 wy1 wz1 Cly1 Clz1];
##fprintf('\n\n\n')
##for i = 1:length(z)
##  dif1 = simplify(diff(Lf, z(i)));
##  for j = 1:length(z)
##    blocklf = simplify(diff(dif1, z(j)));
##    if blocklf ~= 0
##      if i > j
##        fprintf('hes_block(%d,%d) = hes_block(%d,%d);\n',i,j,j,i);
##      else
##        fprintf('hes_block(%d,%d) = %s;\n',i,j,char(blocklf));
##      end
##    end
##  end
##end

##fprintf('\n\n\n')
##for i = 1:length(z)
##  dif1 = simplify(diff(Lm, z(i)));
##  for j = 1:length(z)
##    blocklf = simplify(diff(dif1, z(j)));
##    if blocklf ~= 0
##      if i > j
##        fprintf('hes_block(%d,%d) = hes_block(%d,%d);\n',i,j,j,i);
##      else
##        fprintf('hes_block(%d,%d) = %s;\n',i,j,char(blocklf));
##      end
##    end
##  end
##end

z = [pn2 pe2 pd2 s2 q02 q12 q22 q32];
for i = 1:length(J)
  for j = 1:length(z)
    block_con = simplify(diff(J(i), z(j)));
    if block_con ~= 0
      fprintf('jac(%d,%d) = %s;\n',i,j,char(block_con));
    end
  end
end
fprintf('\n\n\n')
