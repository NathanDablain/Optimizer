clear
close all
clc
addpath('../');

% Lets say I want to move a rocket to 50km altitude in 120seconds

N_knots = 50;
% height, velocity, mass
N_states = 3;
% thrust
N_inputs = 1;

tf = 120;
t = linspace(0, tf, N_knots);
xd = 50e3;
IC = [0 0 425e3];

dim = (N_states+N_inputs)*N_knots;

lb_states = [-inf(N_knots,1), -inf(N_knots,1), 25e3.*ones(N_knots,1)];
ub_states = [inf(N_knots,1), inf(N_knots,1), inf(N_knots,1)];
lb_inputs = zeros(N_knots,1);
ub_inputs = 8e6.*ones(N_knots,1);

ic_states = [zeros(N_knots,1), zeros(N_knots,1), IC(3).*ones(N_knots,1)];
ic_inputs = 0.75*ub_inputs;

[states_guess, inputs_guess] = build_traj(t, ic_states(1,:), ic_inputs);

% We dont pack slack
X0 = pack_trap(states_guess, inputs_guess, []);
lb = pack_trap(lb_states, lb_inputs, []);
ub = pack_trap(ub_states, ub_inputs, []);

ub_index = [];
lb_index = [];
slack_vars = [];

for i = 1:length(lb)
  if lb(i) ~= -inf
    lb_index = [lb_index i];
    slack_vars = [slack_vars X0(i)-lb(i)];
  end
end

for i = 1:length(ub)
  if ub(i) ~= inf
    ub_index = [ub_index i];
    slack_vars = [slack_vars ub(i)-X0(i)];
  end
end
X0 = pack_trap(states_guess, inputs_guess, slack_vars);

cons = struct('N_knots', N_knots, 'N_states', N_states, 'N_inputs', N_inputs, 'N_sv', length(slack_vars),...
              't', t, 'IC', IC, 'lb', lb, 'ub', ub, 'lb_index', lb_index, 'ub_index', ub_index, 'xd', xd);

[ineq, eq, Jac] = constraints_trap(X0, cons);

tic
[x_opt, J] = kkt(@(x)cost(x, cons), X0, cons, @(x)constraints_trap(x, cons), true, true);
toc

[states, inputs] = unpack_trap(x_opt, N_knots, N_states, N_inputs);

figure(1)
subplot(3,1,1)
plot(t, states_guess(:,1))
hold on
plot(t, states(:,1))
legend('First pass', 'Second pass')

subplot(3,1,2)
plot(t, states_guess(:,2))
hold on
plot(t, states(:,2))

subplot(3,1,3)
plot(t, inputs_guess(:,1))
hold on
plot(t, inputs(:,1))