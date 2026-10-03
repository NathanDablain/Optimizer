clear
clc

pkg load symbolic

syms q01 q11 q21 q31 q02 q12 q22 q32 wy1 wz1 wy2 wz2 h real

##wy = 0.5*(wy1 + wy2);
##wz = 0.5*(wz1 + wz2);
##w = [wy wz];
##
##theta = norm(w)*h;
##
##delta_q = [cos(theta/2);0;(w(1)/norm(w))*sin(theta/2);(w(2)/norm(w))*sin(theta/2)];
##q2 = [q02; q12; q22; q32];
##constraint = q2 -   [q01 -q11 -q21 -q31;
##                     q11  q01 -q31  q21;
##                     q21  q31  q01 -q11;
##                     q31 -q21  q11  q01]*delta_q;
##z = [q01 q11 q21 q31 wy1 wz1 q02 q12 q22 q32 wy2 wz2];
##for i = 1:length(constraint)
##  for j = 1:length(z)
##    block_con = simplify(diff(constraint(i), z(j)));
##    fprintf('jac_block(%d,%d) = %s;\n',i,j,char(block_con));
##  end
##end



