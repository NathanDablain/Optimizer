clc
clear
close all

N_states = 6;
N_inputs = 2;
N_knots = 100;

tf = 20;
t = linspace(0, tf, N_knots);
xd = [100 200 -100];
IC = [0 0 0 50 0 0];

dim = (N_states+N_inputs)*N_knots;

lb_states = [-inf(N_knots,1), -inf(N_knots,1), -inf(N_knots,1),...
             zeros(N_knots,1), -pi*ones(N_knots,1) -pi/2*ones(N_knots,1)];
ub_states = [inf(N_knots,1), inf(N_knots,1), inf(N_knots,1),...
             inf(N_knots,1), pi*ones(N_knots,1), pi/2*ones(N_knots,1)];
lb_inputs = [-0.5*ones(N_knots,1), -0.5*ones(N_knots,1)];
ub_inputs = [0.5*ones(N_knots,1), 0.5*ones(N_knots,1)];

ic_states = [IC(1)*ones(N_knots,1), IC(2)*ones(N_knots,1), IC(3)*ones(N_knots,1),...
             IC(4)*ones(N_knots,1), IC(5)*ones(N_knots,1), IC(6)*ones(N_knots,1)];
ic_inputs = [0.1*ones(N_knots,1) zeros(N_knots,1)];

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

knot_size = N_states + N_inputs;
N_colloc_constraints = N_states*(N_knots-1);
N_ic_constraints = N_states;
N_bound_constraints = length(slack_vars);
N_constraints = N_colloc_constraints + N_ic_constraints + N_bound_constraints;
N_decision_variables = length(X0);

state_scaling = [1e3 1e3 1e3 5e2 pi pi];
slack_scaling = [5e2 pi pi 0.5 0.5 pi pi 0.5 0.5];
scaling = zeros(N_constraints,1);
for i = 1:N_knots
  start_index = N_states*(i-1) + 1;
  end_index = start_index + length(state_scaling) - 1;
  scaling(start_index:end_index) = state_scaling';
end
for i = 1:length(slack_scaling)
  start_index = end_index+1;
  end_index = start_index + N_knots - 1;
  scaling(start_index:end_index) = slack_scaling(i)*ones(N_knots,1);
end

cons = struct('N_knots', N_knots, 'N_states', N_states, 'N_inputs', N_inputs, 'N_sv', length(slack_vars),...
              't', t, 'IC', IC, 'lb', lb, 'ub', ub, 'lb_index', lb_index, 'ub_index', ub_index, 'xd', xd,...
              'N_con', N_constraints, 'N_dec', N_decision_variables, 'H_block', knot_size^2,...
              'J_block', N_states*(2*knot_size), 'N_ic', N_states, 'scaling', scaling);

[~, eq] = constraints_trap(X0, cons);

##tic
##[x_opt, J] = kkt(X0, cons);
##toc
##
##[states, inputs] = unpack_trap(x_opt, N_knots, N_states, N_inputs);

figure()
subplot(3,1,1)
plot(t, states_guess(:,1))
subplot(3,1,2)
plot(t, states_guess(:,2))
subplot(3,1,3)
plot(t, -states_guess(:,3))

figure()
subplot(3,1,1)
plot(t, states_guess(:,4))
subplot(3,1,2)
plot(t, states_guess(:,5))
subplot(3,1,3)
plot(t, states_guess(:,6))
