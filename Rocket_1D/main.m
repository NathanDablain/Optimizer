clear
close all
clc
addpath('../');

% Lets say I want to move a rocket to 50km altitude in 120seconds

N_knots = 30;
% height, velocity, mass
N_states = 3;
% thrust
N_inputs = 1;

tf = 120;
t = linspace(0, tf, N_knots);
x1_f = 50e3;
IC = [0 0 425e3];

method = 'Trapezoidal';
if method == 'Trapezoidal'
  dim = (N_states+N_inputs)*N_knots;
end

lb_states = [zeros(N_knots,1), zeros(N_knots,1), 25e3.*ones(N_knots,1)];
ub_states = [inf(N_knots,1), inf(N_knots,1), IC(3).*ones(N_knots,1)];
lb_inputs = zeros(N_knots,1);
ub_inputs = 8e6.*ones(N_knots,1);

ic_states = [zeros(N_knots,1), zeros(N_knots,1), IC(3).*ones(N_knots,1)];
ic_inputs = ub_inputs;

[states_guess, inputs_guess] = build_traj(t, ic_states(1,:), ic_inputs);


X0 = pack_trap(states_guess, inputs_guess);
lb = pack_trap(lb_states, lb_inputs);
ub = pack_trap(ub_states, ub_inputs);
for i = 1:dim
  if X0(i) < lb(i)
    disp(['State ' num2str(i) ' Violates lower bound of ' num2str(lb(i))]);
  elseif X0(i) > ub(i)
    disp(['State ' num2str(i) ' Violates upper bound of ' num2str(ub(i))]);
  endif
endfor

[ineq, eq] = constraints_trap(X0, t, N_states, N_inputs, IC);

tic
[x_opt2, J2] = fmincon(@(x)cost_step1(x, N_knots), X0, [], [], [], [], lb, ub, @(x)constraints_trap(x, t, N_states, N_inputs, IC));
toc

[states_2, inputs_2] = unpack_trap(x_opt2, N_knots, N_states, N_inputs);

figure(1)
subplot(3,1,1)
plot(t, states_guess(:,1))
hold on
plot(t, states_2(:,1))
legend('First pass', 'Second pass')

subplot(3,1,2)
plot(t, states_guess(:,2))
hold on
plot(t, states_2(:,2))

subplot(3,1,3)
plot(t, inputs_guess(:,1))
hold on
plot(t, inputs_2(:,1))
