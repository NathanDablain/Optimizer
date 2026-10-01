clear
clc

pkg load symbolic

syms mu q0 q1 q2 q3 ub lb u real

J = (q0^2 + q1^2 + q2^2 + q3^2 - 1.0)^2

simplify(diff(J,q0))

