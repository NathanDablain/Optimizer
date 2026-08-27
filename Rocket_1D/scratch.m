clc
clear
close all

% syms r v m T h lambda1 lambda2 lambda3 rho0 r0 C_D A C gravity real
% 
% rho = rho0*exp(-r/r0);
% L = lambda1*(-0.5*h*v) +...
%   lambda2*(-0.5*h*((T - 0.5*C_D*rho*A*v^2)/m - gravity)) +...
%   lambda3*(-0.5*h*(-T/C));
%
% row4 = diff(L, T);
% entry1 = simplify(diff(row4, r))
% entry2 = simplify(diff(row4, v))
% entry3 = simplify(diff(row4, m))
% entry4 = simplify(diff(row4, T))


syms r1 v1 m1 T1 r2 v2 m2 T2 h rho0 r0 C_D A C gravity real

rho1 = rho0*exp(-r1/r0);
rho2 = rho0*exp(-r2/r0);

c = v2 - v1 - 0.5*h*(((T1 - 0.5*C_D*A*rho1*v1^2)/m1 - gravity) + ((T2 - 0.5*C_D*A*rho2*v2^2)/m2 - gravity));
entry1 = simplify(diff(c, r1))
entry2 = simplify(diff(c, v1))
entry3 = simplify(diff(c, m1))
entry4 = simplify(diff(c, T1))
entry5 = simplify(diff(c, r2))
entry6 = simplify(diff(c, v2))
entry7 = simplify(diff(c, m2))
entry8 = simplify(diff(c, T2))

