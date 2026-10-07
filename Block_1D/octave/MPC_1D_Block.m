% We want to show how this would work in real time in an MPC context
% So will need to run an outer sim loop, say 1000Hz
% Run an inner Flight Controller loop, say 100Hz
% And run an inner inner SQP loop, say 10 Hz
% System gets some amount of initial iterations to converge and then gets to only run one iteration per inner inner loop

clear
close all
clc

sim_rate = 1000;
fc_rate = 100;
sqp_rate = 10;
tf = 2.0;
t = linspace(0, tf, tf*sim_rate);
x_true = zeros(2,length(t));
u_true = zeros(1,length(t));

t_sqp = linspace(0, tf, tf*sqp_rate);
pf = zeros(1,tf*sqp_rate);
df = zeros(1,tf*sqp_rate);
position_mean = 0.0001;
position_sigma = 0.00;
velocity_mean = 0.0;
velocity_sigma = 0.0;
acceleration_mean = 0.0;
acceleration_sigma = 0.0;

O = struct(...
  'N_states', 2,...
  'N_inputs', 1,...
  'N_knots', 10,...
  'N_lb', 1,...
  'N_ub', 1,...
  'ic', [0 0.0 0.0],...
  'xd', [2 0.4],...
  'lb', -2.0,...
  'ub', 2.0,...
  'tf', tf,...
  't', [],...
  'knot_size', [],... % derived
  'N_decision_variables', [],... % derived
  'First_knot_constraints', 5,...
  'Middle_knot_constraints', 2,...
  'End_knot_constraints', 0,...
  'N_constraints', [],... % derived
  'z', [],...
  'lambda', [],...
  'pf', 0,...
  'df', 0,...
  'iterations', [],...
  'states', [],...
  'inputs', []...
  );

% Set the derived quantities
O.knot_size = O.N_states +...
              O.N_inputs;

O.N_decision_variables = O.knot_size * O.N_knots;

O.N_constraints = O.First_knot_constraints +...
                  O.Middle_knot_constraints*(O.N_knots - 2) +...
                  O.End_knot_constraints;

O.z = zeros(O.N_decision_variables, 1);
O.lambda = zeros(O.N_constraints, 1);
O.t = linspace(0, O.tf, O.N_knots);

O.z(1:length(O.ic)) = O.ic';
O.states = zeros(2,O.N_knots);
O.inputs = zeros(1,O.N_knots);

% Initial solve
O.iterations = 5;
O = Solve_1D_Block(O);
u_command = O.inputs(1);
Ideal_states = O.states;
Ideal_inputs = O.inputs;
Ideal_time = O.t;
pf_counter = 1;
pf(pf_counter) = O.pf;
df(pf_counter) = O.df;
pf_counter = pf_counter+1;

% Simulation loop (1000Hz)
for i = 2:length(t)

  % SQP loop (50Hz)
  if mod(i,sim_rate/sqp_rate) == 0 && i ~= length(t)
    O.iterations = 1;
    O.t = linspace(t(i-1), O.tf, O.N_knots);
    O.ic = [x_true(:,i-1); min(max(u_true(1,i-1), -1.999), 1.999)];
    O.z(1:length(O.ic)) = O.ic';
    O = Solve_1D_Block(O);
    pf(pf_counter) = O.pf;
    df(pf_counter) = O.df;
    pf_counter = pf_counter+1;
  end

  % Flight Controller loop 100Hz
  if mod(i,sim_rate/fc_rate) == 0
    u_command = interp1(O.t, O.inputs, t(i)+(0.5/fc_rate));
  end

  % Noise
  position_noise = position_mean + (position_sigma*randn());
  velocity_noise = velocity_mean + (velocity_sigma*randn());
  acceleration_noise = acceleration_mean + (acceleration_sigma*randn());

  % System dynamics
  dt = t(i) - t(i-1);
  x_true(1,i) = x_true(1,i-1) + dt*x_true(2,i-1) + position_noise;
  x_true(2,i) = x_true(2,i-1) + dt*u_true(1,i-1) + velocity_noise;
  u_true(1,i) = u_command + acceleration_noise;

end


% Plot output
figure()
subplot(3,1,1)
plot(t, x_true(1,:))
hold on
plot(Ideal_time, Ideal_states(1,:))

subplot(3,1,2)
plot(t,x_true(2,:))
hold on
plot(Ideal_time, Ideal_states(2,:))

subplot(3,1,3)
plot(t,u_true(1,:))
hold on
plot(Ideal_time, Ideal_inputs(1,:))

figure()
plot(t_sqp, pf)
hold on
plot(t_sqp, df);
legend('Primal Feasability', 'Dual Feasability')
