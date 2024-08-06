% MATLAB: Load variables and run FFSP and particle filter
%%
% Load data from Python
data = load('data_for_matlab_1.mat');
%% 
% Extract parameters and initial conditions from the loaded data
initial_conditions = data.initial_conditions;
parameter_values = data.parameter_values;
time_of_observations = data.time_of_observations;
time_of_observations_1=time_of_observations(1:end-1);
state_sequence = data.state_sequence;
state_sequence_1= state_sequence(1:end-1);
observation_func = data.observation_func;
time_list = data.time_list;
hidden_trajectories = data.hidden_trajectories;

%%
% Initialize parameters from Python data
c1 = double(parameter_values.c1);
c2 = double(parameter_values.c2);
c3 = double(parameter_values.c3);
c4 = double(parameter_values.c4);
c5 = double(parameter_values.c5);
k1 = double(parameter_values.k1);
k2 = double(parameter_values.k2);
k3 = double(parameter_values.k3);
k4 = double(parameter_values.k4);

Omega = 10;
N = 2;  % Number of hidden species

% Initial conditions
z10 = double(initial_conditions.z1);
z20 = double(initial_conditions.z2);
x10 = double(initial_conditions.x1);

%%
% Reaction rate parameters vector C
C = [c1; c2; c3; c4; c5; k1; k2; k3; k4];

MAK = 0; % Assuming mass action kinetics parameter to be zero

% Define propensities
W_star{1} = @(z1, z2, x1) C(1);
W_star{2} = @(z1, z2, x1) C(2);
W_star{3} = @(z1, z2, x1) C(3) .* z1 .* z2; % Ensure element-wise multiplication
W_star{4} = @(z1, z2, x1) C(4) + (C(8) * (z1 / C(6))) / (C(9) + (z1 / C(6)) + (z2 / C(7))); % Ensure element-wise multiplication
W_star{5} = @(z1, z2, x1) C(5) .* x1; % Ensure element-wise multiplication

% Debugging step: Check the dimensions of the propensity functions
test_z1 = 1;
test_z2 = 1;
test_x1 = 1;
disp(W_star{1}(test_z1, test_z2, test_x1)); % Should be scalar
disp(W_star{2}(test_z1, test_z2, test_x1)); % Should be scalar
disp(W_star{3}(test_z1, test_z2, test_x1)); % Should be scalar
disp(W_star{4}(test_z1, test_z2, test_x1)); % Should be scalar
disp(W_star{5}(test_z1, test_z2, test_x1)); % Should be scalar

% Additional debugging: Verify inputs and outputs of W_star functions during execution
debug_W_star = @(f, z1, z2, x1) disp(['Inputs: z1=', num2str(z1), ', z2=', num2str(z2), ', x1=', num2str(x1), ' | Output: ', num2str(f(z1, z2, x1))]);

% Run some test cases with debugging
for i = 1:5
    disp(['Debugging W_star{', num2str(i), '}']);
    debug_W_star(W_star{i}, test_z1, test_z2, test_x1);
end

% Stoichiometric matrices
S = [1, 0, -1, 0, 0; 0, 1, -1, 0, 0; 0, 0, 0, 1, -1];
S_bis = [0, 0, 1, 0, 0; 0, 0, 1, 0, 0; 0, 0, 0, 0, 1];

% Hidden state space
N_state = 40;
Hidden_Species = Hidden_State(0:N_state, 0:N_state); % Adjust to your function definition
Shape_State_Space_Hidden = size(Hidden_Species);
State_Hidden = Shape_State_Space_Hidden(1);
[~, index] = intersect(Hidden_Species, [z10 * Omega; z20 * Omega]', 'rows');
p0_1 = zeros(State_Hidden, 1);
p0_1(index) = 1;

%%
% Run FFSP
delta = diff(state_sequence_1');
[T, F, jump_times, E_FSP, Var_FSP, SD_FSP, Err_jump, E_tot, SD_tot, rho] = ...
    FFSP_2(time_of_observations_1, state_sequence_1', p0_1, C, Hidden_Species, delta, S, S_bis, N, W_star, MAK);

%%
% Particle filter run
tf = time_of_observations(end);
N_tot = 1000;
resample = 1;

[V_tot_1, w_tot, V_jump, w_jump, match_1, match_2, resampling, pt, E_PF, Var_PF, SD_PF] = ...
    particle_filter_1(time_of_observations_1, state_sequence_1', p0_1, C, tf, N_tot, Hidden_Species, S, S_bis, resample, N, W_star, MAK);

%%
% Prepare data for plotting
E_PF = [E_PF; E_PF(end,:)];
SD_PF = [SD_PF; SD_PF(end,:)];
time_of_observations_1 = [time_of_observations_1, tf];
E_FSP = [E_FSP; E_FSP(end,:)];
SD_FSP = [SD_FSP; SD_FSP(end,:)];
%%
% Plotting results for each hidden species
for i = 1:N
    figure;
    hold on;
    stairs(time_list, hidden_trajectories(:, i), 'LineWidth', 2, 'DisplayName', 'Hidden SSA');
    errorbar(time_of_observations_1, E_PF(:, i), SD_PF(:, i), 'r--', 'LineWidth', 2, 'DisplayName', 'PF Estimation');
    errorbar(time_of_observations_1, E_FSP(:, i), SD_FSP(:, i), 'b--', 'LineWidth', 2, 'DisplayName', 'FFSP Estimation');
    xlabel('Time (s)');
    ylabel(['Species ' num2str(i)]);
    title(['Estimation for Species ' num2str(i)]);
    legend show;
    hold off;
end

%%
% Save FFSP and particle filter results
save('FFSP_estimations_maybe_1.mat', 'E_FSP', 'SD_FSP');
save('particle_estimations_maybe_1.mat', 'E_PF', 'SD_PF');
