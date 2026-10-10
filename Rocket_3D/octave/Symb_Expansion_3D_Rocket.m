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

function dx = get_dx1(pn1, pe1, pd1, vn1, ve1, vd1, an1, ae1, ad1, c_L1)
  syms c_D1 rho1 A g1 T1 m1
  v_vec = [vn1; ve1; vd1];
  v_uv = v_vec ./ norm(v_vec);
  a_uv = [an1; ae1; ad1];
  g_vec = [0;0;g1];

  D = 0.5*rho1*A*c_D1*norm(v_vec)^2;
  L = 0.5*rho1*A*c_L1*norm(v_vec)^2;
  dx(1:3) = v_vec;
  dx(4:6) = ((T1 - D)*v_uv + L*a_uv)/m1 + g_vec;
end

function dx = get_dx2(pn2, pe2, pd2, vn2, ve2, vd2, an2, ae2, ad2, c_L2)
  syms c_D2 rho2 A g2 T2 m2
  v_vec = [vn2; ve2; vd2];
  v_uv = v_vec ./ norm(v_vec);
  a_uv = [an2; ae2; ad2];
  g_vec = [0;0;g2];

  D = 0.5*rho2*A*c_D2*norm(v_vec)^2;
  L = 0.5*rho2*A*c_L2*norm(v_vec)^2;
  dx(1:3) = v_vec;
  dx(4:6) = ((T2 - D)*v_uv + L*a_uv)/m2 + g_vec;
end

function c = get_constraints_first_knot(states, inputs, ic)
  syms h

  pn1 = states(1,1);
  pe1 = states(1,2);
  pd1 = states(1,3);
  vn1 = states(1,4);
  ve1 = states(1,5);
  vd1 = states(1,6);

  an1 = inputs(1,1);
  ae1 = inputs(1,2);
  ad1 = inputs(1,3);
  c_L1 = inputs(1,4);

  pn_ic = ic(1);
  pe_ic = ic(2);
  pd_ic = ic(3);
  vn_ic = ic(4);
  ve_ic = ic(5);
  vd_ic = ic(6);

  pn2 = states(2,1);
  pe2 = states(2,2);
  pd2 = states(2,3);
  vn2 = states(2,4);
  ve2 = states(2,5);
  vd2 = states(2,6);

  an2 = inputs(2,1);
  ae2 = inputs(2,2);
  ad2 = inputs(2,3);
  c_L2 = inputs(2,4);

  dx1 = get_dx1(pn1, pe1, pd1, vn1, ve1, vd1, an1, ae1, ad1, c_L1);
  dx1 = simplify(dx1);
  dx2 = get_dx2(pn2, pe2, pd2, vn2, ve2, vd2, an2, ae2, ad2, c_L2);
  dx2 = simplify(dx2);

  c = [...
    % the trapezoidal defect constraints for the knot
    pn2 - pn1 - 0.5*h*(dx1(1) + dx2(1));...
    pe2 - pe1 - 0.5*h*(dx1(2) + dx2(2));...
    pd2 - pd1 - 0.5*h*(dx1(3) + dx2(3));...
    vn2 - vn1 - 0.5*h*(dx1(4) + dx2(4));...
    ve2 - ve1 - 0.5*h*(dx1(5) + dx2(5));...
    vd2 - vd1 - 0.5*h*(dx1(6) + dx2(6));...
    % the orthogonality constraint
    dot([an1 ae1 ad1],[vn1 ve1 vd1]);...
    % the unit length constraint
    norm([an1 ae1 ad1]) - 1.0;...
    % the initial condition constraints for the knot
    pn_ic - pn1;...
    pe_ic - pe1;...
    pd_ic - pd1;...
    vn_ic - vn1;...
    ve_ic - ve1;...
    vd_ic - vd1];
end

function c = get_constraints_middle_knot(states, inputs)
  syms h

  pn1 = states(1,1);
  pe1 = states(1,2);
  pd1 = states(1,3);
  vn1 = states(1,4);
  ve1 = states(1,5);
  vd1 = states(1,6);

  an1 = inputs(1,1);
  ae1 = inputs(1,2);
  ad1 = inputs(1,3);
  c_L1 = inputs(1,4);

  pn2 = states(2,1);
  pe2 = states(2,2);
  pd2 = states(2,3);
  vn2 = states(2,4);
  ve2 = states(2,5);
  vd2 = states(2,6);

  an2 = inputs(2,1);
  ae2 = inputs(2,2);
  ad2 = inputs(2,3);
  c_L2 = inputs(2,4);

  dx1 = get_dx1(pn1, pe1, pd1, vn1, ve1, vd1, an1, ae1, ad1, c_L1);
  dx1 = simplify(dx1);
  dx2 = get_dx2(pn2, pe2, pd2, vn2, ve2, vd2, an2, ae2, ad2, c_L2);
  dx2 = simplify(dx2);

  c = [...
    % the trapezoidal defect constraints for the knot
    pn2 - pn1 - 0.5*h*(dx1(1) + dx2(1));...
    pe2 - pe1 - 0.5*h*(dx1(2) + dx2(2));...
    pd2 - pd1 - 0.5*h*(dx1(3) + dx2(3));...
    vn2 - vn1 - 0.5*h*(dx1(4) + dx2(4));...
    ve2 - ve1 - 0.5*h*(dx1(5) + dx2(5));...
    vd2 - vd1 - 0.5*h*(dx1(6) + dx2(6));...
    % the orthogonality constraint
    dot([an1 ae1 ad1],[vn1 ve1 vd1]);...
    % the unit length constraint
    norm([an1 ae1 ad1]) - 1.0];
end

function c = get_constraints_end_knot(states, inputs)
  syms h

  pn1 = states(1,1);
  pe1 = states(1,2);
  pd1 = states(1,3);
  vn1 = states(1,4);
  ve1 = states(1,5);
  vd1 = states(1,6);

  an1 = inputs(1,1);
  ae1 = inputs(1,2);
  ad1 = inputs(1,3);
  c_L1 = inputs(1,4);

  c = [...
    % the orthogonality constraint
    dot([an1 ae1 ad1],[vn1 ve1 vd1]);...
    % the unit length constraint
    norm([an1 ae1 ad1]) - 1.0];
end

function j = get_trapezoidal_integration(states, inputs)
  syms h

  pn1 = states(1,1);
  pe1 = states(1,2);
  pd1 = states(1,3);
  vn1 = states(1,4);
  ve1 = states(1,5);
  vd1 = states(1,6);

  an1 = inputs(1,1);
  ae1 = inputs(1,2);
  ad1 = inputs(1,3);
  c_L1 = inputs(1,4);

  pn2 = states(2,1);
  pe2 = states(2,2);
  pd2 = states(2,3);
  vn2 = states(2,4);
  ve2 = states(2,5);
  vd2 = states(2,6);

  an2 = inputs(2,1);
  ae2 = inputs(2,2);
  ad2 = inputs(2,3);
  c_L2 = inputs(2,4);

  dx1 = get_dx1(pn1, pe1, pd1, vn1, ve1, vd1, an1, ae1, ad1, c_L1);
  dx1 = simplify(dx1);
  dx2 = get_dx2(pn2, pe2, pd2, vn2, ve2, vd2, an2, ae2, ad2, c_L2);
  dx2 = simplify(dx2);

  j = [...
    pn2 - pn1 - 0.5*h*(dx1(1) + dx2(1));...
    pe2 - pe1 - 0.5*h*(dx1(2) + dx2(2));...
    pd2 - pd1 - 0.5*h*(dx1(3) + dx2(3));...
    vn2 - vn1 - 0.5*h*(dx1(4) + dx2(4));...
    ve2 - ve1 - 0.5*h*(dx1(5) + dx2(5));...
    vd2 - vd1 - 0.5*h*(dx1(6) + dx2(6))];
end

function lagrangian = get_lagrangian_first_knot(states, inputs, lbs, ubs, ic)
  syms mu lambda1 lambda2 lambda3 lambda4 lambda5 lambda6 lambda7...
          lambda8 lambda9 lambda10 lambda11 lambda12 lambda13 lambda14 real;

  lambda = [lambda1 lambda2 lambda3 lambda4 lambda5 lambda6 lambda7...
            lambda8 lambda9 lambda10 lambda11 lambda12 lambda13 lambda14];

  cost = -mu*(log(inputs(1,4) - lbs(1)) + log(ubs(1) - states(1,3)) + log(ubs(2) - inputs(1,4)));

  cf = get_constraints_first_knot(states, inputs, ic);
  lagrangian = cost + lambda*cf;
end

function lagrangian = get_lagrangian_middle_knot(states, inputs, lbs, ubs)
  syms mu mu_q lambda1 lambda2 lambda3 lambda4 lambda5 lambda6 lambda7 lambda8 real;

  lambda = [lambda1 lambda2 lambda3 lambda4 lambda5 lambda6 lambda7 lambda8];
  cost = -mu*(log(inputs(1,4) - lbs(1)) + log(ubs(1) - states(1,3)) + log(ubs(2) - inputs(1,4)));

  cm = get_constraints_middle_knot(states, inputs);
  lagrangian = cost + lambda*cm;
end

syms pn1 pe1 pd1 vn1 ve1 vd1 an1 ae1 ad1 c_L1...
     pn2 pe2 pd2 vn2 ve2 vd2 an2 ae2 ad2 c_L2...
     lb_c_L ub_pd ub_c_L...
     pn_ic pe_ic pd_ic vn_ic ve_ic vd_ic xd1 xd2 xd3 real

states = [pn1 pe1 pd1 vn1 ve1 vd1;...
          pn2 pe2 pd2 vn2 ve2 vd2];
inputs = [an1 ae1 ad1 c_L1;...
          an2 ae2 ad2 c_L2];
lbs = [lb_c_L];
ubs = [ub_pd ub_c_L];
ic = [pn_ic pe_ic pd_ic vn_ic ve_ic vd_ic];
xd = [xd1 xd2 xd3];

z = [pn1 pe1 pd1 vn1 ve1 vd1 an1 ae1 ad1 c_L1 pn2 pe2 pd2 vn2 ve2 vd2 an2 ae2 ad2 c_L2];

##cf = simplify(get_constraints_first_knot(states, inputs, ic));

##cm = simplify(get_constraints_middle_knot(states, inputs));
##
##ce = simplify(get_constraints_end_knot(states, inputs));

##Lf = get_lagrangian_first_knot(states, inputs, lbs, ubs, ic);

##Lm = simplify(get_lagrangian_middle_knot(states, inputs, lbs, ubs));

J = get_trapezoidal_integration(states, inputs);

clc
##disp(cf);
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
##z = [pn1 pe1 pd1 vn1 ve1 vd1 an1 ae1 ad1 c_L1];
##for i = 1:length(ce)
##  for j = 1:length(z)
##    block_con = simplify(diff(ce(i), z(j)));
##    if block_con ~= 0
##      fprintf('jac_block(%d,%d) = %s;\n',i,j,char(block_con));
##    end
##  end
##end
##fprintf('\n\n\n')

##z = [pn1 pe1 pd1 vn1 ve1 vd1 an1 ae1 ad1 c_L1];
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

z = [pn2 pe2 pd2 vn2 ve2 vd2];
for i = 1:length(J)
  for j = 1:length(z)
    block_con = simplify(diff(J(i), z(j)));
    if block_con ~= 0
      fprintf('jac(%d,%d) = %s;\n',i,j,char(block_con));
    end
  end
end
fprintf('\n\n\n')
