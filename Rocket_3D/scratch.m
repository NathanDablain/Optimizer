clear
clc
pkg load symbolic

##syms rn re rd s psi gam Clp Clg Thrust m rho c_D A g rho0 r0 h lambda1 lambda2 lambda3 lambda4 lambda5 lambda6 real
##
##rho = rho0*exp(rd/r0);
##Drag = 0.5*rho*c_D*A*s*s;
##L = -0.5*h*(lambda1*(s * cos(gam) * cos(psi)) +...
##    lambda2*(s * cos(gam) * sin(psi)) +...
##    lambda3*(-s * sin(gam)) +...
##    lambda4*((Thrust - Drag)/m - g*sin(gam)) +...
##    lambda5*((rho * A * s * Clp) / (2 * cos(gam))) +...
##    lambda6*((0.5 * rho * A * (s^2) * Clg - g * cos(gam)) / s));
##
##row8 = diff(L, Clg);
##entry1 = simplify(diff(row8, rn))
##entry2 = simplify(diff(row8, re))
##entry3 = simplify(diff(row8, rd))
##entry4 = simplify(diff(row8, s))
##entry5 = simplify(diff(row8, psi))
##entry6 = simplify(diff(row8, gam))
##entry7 = simplify(diff(row8, Clp))
##entry8 = simplify(diff(row8, Clg))

##
##  dx(1) = s * cos(gam) * cos(psi);
##  dx(2) = s * cos(gam) * sin(psi);
##  dx(3) = -s * sin(gam);
##  dx(4) = (Thrust - Drag)/m - gravity*sin(gam);
##  dx(5) = (rho * A * s * input(1)) / (2 * cos(gam));
##  dx(6) = (0.5 * rho * A * (s^2) * input(2) - gravity * cos(gam)) / s;

syms rn1 rn2 re1 re2 rd1 rd2 s1 s2 psi1 psi2 gam1 gam2 Clp1 Clp2 Clg1 Clg2 T1 T2 m1 m2 rho1 rho2 c_D A g rho0 r0 h real

rho1 = rho0*exp(rd1/r0);
rho2 = rho0*exp(rd2/r0);
D1 = 0.5 * rho1 * A * c_D * s1 * s1;
D2 = 0.5 * rho2 * A * c_D * s2 * s2;

c = gam2 - gam1 - 0.5*h*((0.5*rho1*A*s1*s1*Clg1 - g*cos(gam1))/s1 + (0.5*rho2*A*s2*s2*Clg2 - g*cos(gam2))/s2);

entry1 = simplify(diff(c, rn1))
entry2 = simplify(diff(c, re1))
entry3 = simplify(diff(c, rd1))
entry4 = simplify(diff(c, s1))
entry5 = simplify(diff(c, psi1))
entry6 = simplify(diff(c, gam1))
entry7 = simplify(diff(c, Clp1))
entry8 = simplify(diff(c, Clg1))
entry9 = simplify(diff(c, rn2))
entry10 = simplify(diff(c, re2))
entry11 = simplify(diff(c, rd2))
entry12 = simplify(diff(c, s2))
entry13 = simplify(diff(c, psi2))
entry14 = simplify(diff(c, gam2))
entry15 = simplify(diff(c, Clp2))
entry16 = simplify(diff(c, Clg2))
