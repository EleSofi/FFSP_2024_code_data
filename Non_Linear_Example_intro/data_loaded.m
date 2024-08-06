% MATLAB Script to Load Data and Plot Observation and Filter Results

% Load the data saved from Python
data = load('data_for_matlab_with_filters_0013.mat');

% Extract the data
initial_conditions = data.initial_conditions;
parameter_values = data.parameter_values;
time_of_observations = data.time_of_observations;
state_sequence = data.state_sequence;
observation_func = data.observation_func;
time_list = data.time_list; % This has 63 points
hidden_trajectories = data.hidden_trajectories;
common_time_points = data.common_time_points; % This has 1001 points

kalman_trajectory = data.hat_z_t_trajectory;
kalman_R_trajectory = data.R_t_trajectory;
ekf_trajectory = data.hat_z_t_trajectory_1;
ekf_R_trajectory = data.R_t_trajectory_1;
E_FSP_interp = data.E_FSP_interp;
SD_FSP_interp = data.SD_FSP_interp;
E_PF_interp = data.E_PF_interp;
SD_PF_interp = data.SD_PF_interp;

hidden_species = {'Z_1', 'Z_2'};
observed_species = {'X_1'};

% Set point density
point_density = 0.1;
