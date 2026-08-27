clear
close all
clc
addpath('../');

% Lets say I want to move a block 2 meters in 2 seconds and have zero speed at 2 seconds

N_knots = 20;
N_states = 2;
N_inputs = 1;
tf = 10;
t = linspace(0, tf, N_knots);
IC = [0 0];
xd = [2 0];

dim = (N_states+N_inputs)*N_knots;
state_lb = [-inf(N_knots,1) -inf(N_knots,1)];
input_lb = [-2000*ones(N_knots,1)];
state_ub = [inf(N_knots,1) inf(N_knots,1)];
input_ub = [2000*ones(N_knots,1)];

ic_states = [zeros(N_knots,1), zeros(N_knots,1)];
ic_inputs = zeros(N_knots,1);

X0 = pack_trap(ic_states, ic_inputs, []);
ub = pack_trap(state_ub, input_ub, []);
lb = pack_trap(state_lb, input_lb, []);

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
X0 = pack_trap(ic_states, ic_inputs, slack_vars);

cons = struct('N_knots', N_knots, 'N_states', N_states, 'N_inputs', N_inputs, 'N_sv', length(slack_vars),...
              't', t, 'IC', IC, 'lb', lb, 'ub', ub, 'lb_index', lb_index, 'ub_index', ub_index, 'xd', xd);
[ineq, eq, Jac] = constraints_trap(X0, cons);

tic
[x_opt, J] = kkt(@(x)cost(x, cons), X0, cons, @(x)constraints_trap(x, cons), true, true);
toc
[states, inputs] = unpack_trap(x_opt, N_knots, N_states, N_inputs);

figure(1)
subplot(3,1,1)
plot(t, states(:,1), 'LineWidth', 1.5)
ylabel('Position')

subplot(3,1,2)
plot(t, states(:,2), 'LineWidth', 1.5)
ylabel('Velocity')

subplot(3,1,3)
plot(t, inputs(:,1), 'LineWidth', 1.5)
ylabel('Acceleration')
