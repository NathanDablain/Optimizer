clear
close all
clc
addpath('../');

% Lets say I want to move a block 2 meters in 2 seconds and have zero speed at 2 seconds

N_knots = 10;
N_states = 2;
N_inputs = 1;
tf = 10;
t = linspace(0, tf, N_knots);
x1_f = 2.0;
x2_f = 0.0;
IC = [0 0];

method = 'Trapezoidal';
if method == 'Trapezoidal'
  dim = (N_states+N_inputs)*N_knots;
end

X0 = zeros(dim, 1);
lb = -inf(dim,1);
ub = inf(dim,1);

if method == 'Trapezoidal'
  tic
  [x_opt1, J] = fmincon(@(x)cost_step1(x, N_knots), X0, [], [], [], [], lb, ub, @(x)constraints_trap(x, t, N_states, N_inputs, IC));
  toc
end
[states_1, inputs_1] = unpack_trap(x_opt1, N_knots, N_states, N_inputs);

error_x1 = x1_f - x_opt1(N_knots);
error_x2 = x2_f - x_opt1(2*N_knots);

x1_error_threshold = 0.1;
x2_error_threshold = 0.1;

if abs(error_x1) < x1_error_threshold && abs(error_x2) < x2_error_threshold
  lb(N_knots) = x1_f - x1_error_threshold;
  ub(N_knots) = x1_f + x1_error_threshold;
  lb(2*N_knots) = x2_f - x2_error_threshold;
  ub(2*N_knots) = x2_f + x2_error_threshold;

  tic
  [x_opt2, J2] = fmincon(@(x)cost_step2(x, N_knots), x_opt1, [], [], [], [], lb, ub, @(x)constraints_trap(x, t, N_states, N_inputs, IC));
  toc
end

[states_2, inputs_2] = unpack_trap(x_opt2, N_knots, N_states, N_inputs);

figure(1)
subplot(3,1,1)
plot(t, states_1(:,1))
hold on
plot(t, states_2(:,1))
legend('First pass', 'Second pass')

subplot(3,1,2)
plot(t, states_1(:,2))
hold on
plot(t, states_2(:,2))

subplot(3,1,3)
plot(t, inputs_1(:,1))
hold on
plot(t, inputs_2(:,1))
