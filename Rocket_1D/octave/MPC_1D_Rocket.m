clear
close all
clc

function dx = get_dx(z_knot)
  C = 2000.0;
  A = 10.52;
  gravity0 = 9.8065;
  Re = 6371e3;

  g = gravity0*((Re/(Re + z_knot(1)))^2);

  if z_knot(1) < 11000
    c1 = 15.04;
    c2 = 0.00649;
    c3 = 101.29;
    c4 = 273.1;
    c5 = 288.08;
    c6 = 5.256;
    c7 = 0.2869;
    Temperature = c1 - c2*z_knot(1);
    Pressure =c3*(((Temperature+c4)/c5)^c6);

  elseif z_knot(1) >= 11000 && z_knot(1) < 25000
    c1 = -56.46;
    c2 = 22.65;
    c3 = 1.73;
    c4 = 0.000157;
    c5 = 0.2869;
    c6 = 273.1;
    Temperature = c1;
    Pressure = c2*exp(c3 - c4*z_knot(1));

  else
    c1 = -131.21;
    c2 = 0.00299;
    c3 = 2.488;
    c4 = 273.1;
    c5 = 216.6;
    c6 = -11.388;
    c7 = 0.2869;
    Temperature = c1 + c2*z_knot(1);
    Pressure = c3*(((Temperature+c4)/c5)^c6);

  end
  rho = Pressure/(0.2869*(Temperature+273.1));

  speed_of_sound = sqrt(1.4*287*(Temperature+273.1));
  Mach = z_knot(2) / speed_of_sound;
  if Mach < 0.8
    c_D = 0.22;

  elseif Mach >= 0.8 && Mach < 1.2
    c1 = 0.22;
    c2 = 0.48;
    c3 = 0.8;
    c_D = c1 + c2*sin(pi*(Mach - c3)/c3)^2;

  elseif Mach >= 1.2
    c1 = 0.25;
    c2 = 0.54;
    c3 = 1.2;
    c_D = c1 + c2/(Mach^c3);

  end

  Drag = 0.5*c_D*rho*A*z_knot(2)*z_knot(2);

  dx(1) = z_knot(2);
  dx(2) = (z_knot(4) - Drag)/z_knot(3) - g;
  dx(3) = -z_knot(4) / C;
end

sim_rate = 200;
fc_rate = 100;
sqp_rate = 1;
tf = 120.0;
t = linspace(0, tf, tf*sim_rate);
x_true = zeros(3,length(t));
x_true(3,1) = 425e3;
u_true = zeros(1,length(t));

t_sqp = linspace(0, tf, tf*sqp_rate);
pf = zeros(1,tf*sqp_rate);
df = zeros(1,tf*sqp_rate);
position_mean = 0.0;
position_sigma = 0.00;
velocity_mean = 0.0;
velocity_sigma = 0.0;
mass_mean = 0.0;
mass_sigma = 0.0;
thrust_mean = 0.0;
thrust_sigma = 0.0;

O = struct(...
  'N_states', 3,...
  'N_inputs', 1,...
  'N_knots', 50,...
  'N_lb', 2,...
  'N_ub', 1,...
  'ic', [0.0 0.0 425.0e3 6.0e6],...
  'xd', 50e3,...
  'lb', [25.0e3 0.0],...
  'ub', 8e6,...
  'tf', 120.0,...
  't', [],...
  'knot_size', [],... % derived
  'N_decision_variables', [],... % derived
  'First_knot_constraints', 7,...
  'Middle_knot_constraints', 3,...
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
O.states = zeros(3,O.N_knots);
O.inputs = zeros(1,O.N_knots);

% Initial solve
O.iterations = 10;
O = Simulate_Rocket(O);
O = Solve_1D_Rocket(O);
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
    O.ic = [x_true(:,i-1); u_true(1,i-1)];
    O.z(1:length(O.ic)) = O.ic';
    O = Solve_1D_Rocket(O);
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
  mass_noise = mass_mean + (mass_sigma*randn());
  thrust_noise = thrust_mean + (thrust_sigma*randn());

  % System dynamics
  dt = t(i) - t(i-1);
  input = [x_true(1,i-1) + position_noise;...
           x_true(2,i-1) + velocity_noise;...
           x_true(3,i-1) + mass_noise;...
           u_true(1,i-1) + thrust_noise];

  dx = get_dx(input);
  x_true(:,i) = x_true(:,i-1) + dt*dx';
  u_true(1,i) = u_command + thrust_noise;

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
