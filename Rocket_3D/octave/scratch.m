clear
clc

pkg load symbolic

syms Cly Clz dCly dClz alpha ub mu real

Cly_new = Cly + alpha*dCly;
Clz_new = Clz + alpha*dClz;

x = Cly_new^2 + Clz_new^2 == ub^2

S = simplify(solve(x, alpha))

J = -mu*log(ub^2 - Cly^2 - Clz^2);
simplify(diff(diff(J,Cly),Clz))
simplify(diff(diff(J,Clz),Cly))
simplify(diff(diff(J,Cly),Cly))
simplify(diff(diff(J,Clz),Clz))

