% Specific function to problem being solved, describes the continuous time differential equations
function dX = diffeq(state, input)
  dX(1) = state(2);
  dX(2) = input(1);
end
