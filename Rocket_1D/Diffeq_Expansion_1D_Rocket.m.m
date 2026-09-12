clear
clc
close all

pkg load symbolic

syms Re gravity0 h c1 c2 c3 c4 c5 c6 c7 a s real;

g_part_h = -2*Re*gravity0/((Re + h)^3);

simplify(diff(g_part_h, h))

rho_part_h = (c2*c3*((c1 - c2*h + c4)/c5)^c6*...
                             (1 - c6))/(c7*(c1-c2*h+c4)^2);

simplify(diff(rho_part_h, h))

rho_part_h = (-c2*c4*exp(c3 - c4*h))/(c5*(c1+c6));

simplify(diff(rho_part_h, h))

rho_part_h = (c2*c3*((c1 + c2*h + c4)/c5)^c6*...
                             (c6 - 1))/(c7*(c1+c2*h+c4)^2);

simplify(diff(rho_part_h, h))

c_D_part_s = (-c2*c3*((s/a)^-c3))/s

simplify(diff(c_D_part_s, s))

c_D_part_s = (pi*c2*sin(2*pi*s/(a*c3)))/(a*c3)

simplify(diff(c_D_part_s, s))
