import sympy as sp
from matplotlib.ticker import MaxNLocator, FuncFormatter
from scipy.interpolate import interp1d

# Import custom modules using relative imports
from .CRN import *
from .CRN_Simulation.CRN_ContinuousTimeFiltering.CRNForContinuousTimeFiltering import *
# Enable LaTeX rendering globally
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


def generate_and_plot_ssa(CF_not_scaled, parameter_values, initial_conditions, Omega, freq=1):
        scaled_initial_conditions = {key: value * Omega for key, value in initial_conditions.items()}
        trajectory = CF_not_scaled.SSA_rich_data_structure(scaled_initial_conditions, parameter_values, 0, 1)

        # Extract observable trajectories
        Y_trj = CF_not_scaled.extract_Observable_trajectory(trajectory.time_sequence, trajectory.state_sequence)
        time_list_1 = Y_trj['time_list'][::freq]
        state_list_1 = Y_trj['state_list'][::freq]
        time_list_1 = np.array(time_list_1)
        state_list_1 = np.array(state_list_1)

        # Save trajectories
        time_list = np.array(trajectory.time_sequence)
        species_names = CF_not_scaled.species_names
        species_data = {species: np.array([state[idx] for state in trajectory.state_sequence]) for idx, species in
                        enumerate(species_names)}

        observable_species = CF_not_scaled.observable_species
        observed_data = {species: state_list_1[:, idx] for idx, species in enumerate(observable_species)}

        # Plot hidden species
        hidden_species = [species for species in species_names if species not in observable_species]
        all_species = hidden_species + list(observable_species)
        n_species = len(all_species)
        ncols = 2
        nrows = (n_species + ncols - 1) // ncols
        fig, axes = plt.subplots(nrows, ncols, figsize=(10, nrows * 5), constrained_layout=True)
        axes = np.atleast_1d(axes).flatten()
        colors = plt.get_cmap('tab10').colors

        for idx, species in enumerate(hidden_species):
            ax = axes[idx]
            ax.step(time_list, species_data[species], label=f'{species} SSA', linewidth=2,
                    color=colors[idx % len(colors)], where='post')
            ax.legend(fontsize=12)
            ax.xaxis.set_major_locator(MaxNLocator(nbins=10, prune='both'))
            ax.yaxis.set_major_locator(MaxNLocator(integer=True, prune='both'))
            ax.tick_params(axis='both', which='major', labelsize=12)
            ax.set_xlabel("Time (s)", fontsize=12)
            ax.set_ylabel(species, fontsize=12)
            ax.set_aspect('auto')

        for idx, species in enumerate(observable_species, start=len(hidden_species)):
            ax = axes[idx]
            ax.step(time_list_1, observed_data[species], label=f'{species} observed', linewidth=2,
                    color=colors[idx % len(colors)], where='post')
            ax.legend(fontsize=12)
            ax.xaxis.set_major_locator(MaxNLocator(nbins=10, prune='both'))
            ax.yaxis.set_major_locator(MaxNLocator(integer=True, prune='both'))
            ax.tick_params(axis='both', which='major', labelsize=12)
            ax.set_xlabel("Time (s)", fontsize=12)
            ax.set_ylabel(species, fontsize=12)
            ax.set_aspect('auto')

        for ax in axes[len(all_species):]:
            fig.delaxes(ax)

        plt.show()
        return species_data, observed_data, time_list, time_list_1


# Compute and plot statistics from SSA trajectories
def compute_and_plot_statistics_SSA(CF_not_scaled, all_trajectories, T):
    common_time_points = np.linspace(0, T, 1001)
    num_species = len(all_trajectories[0].state_sequence[0])

    interpolated_states = [[] for _ in range(num_species)]
    interpolated_states_square = [[] for _ in range(num_species)]

    for trj in all_trajectories:
        for i in range(num_species):
            interp_func = interp1d(trj.time_sequence, np.array([state[i] for state in trj.state_sequence]),
                                   kind='linear', fill_value="extrapolate")
            interpolated_states[i].append(interp_func(common_time_points))
            interp_func_square = interp1d(trj.time_sequence, np.array([state[i] for state in trj.state_sequence]) ** 2,
                                          kind='linear', fill_value="extrapolate")
            interpolated_states_square[i].append(interp_func_square(common_time_points))

    statistics = {}
    for i, species in enumerate(CF_not_scaled.species_names):
        interpolated_states_np = np.array(interpolated_states[i])
        mean_states = np.mean(interpolated_states_np, axis=0)
        interpolated_states_square_np = np.array(interpolated_states_square[i])
        mean_square_states = np.mean(interpolated_states_square_np, axis=0)
        variance_states = mean_square_states - mean_states ** 2
        statistics[species] = {'mean': mean_states, 'variance': variance_states}

    # Plot statistics
    n_species = len(statistics)
    ncols = 2
    nrows = (n_species + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(15, nrows * 7.5), constrained_layout=True)
    axes = np.atleast_1d(axes).flatten()
    colors = plt.get_cmap('tab10').colors

    for idx, (ax, species) in enumerate(zip(axes[:n_species], statistics.keys())):
        mean_states = statistics[species]['mean']
        variance_states = statistics[species]['variance']
        ax.plot(common_time_points, mean_states, label=f'{species} mean', linewidth=2, color=colors[idx % len(colors)])
        ax.fill_between(common_time_points, mean_states - np.sqrt(variance_states),
                        mean_states + np.sqrt(variance_states), color=colors[idx % len(colors)], alpha=0.3,
                        label=f'{species} variance')
        ax.legend(fontsize=12)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=10, prune='both'))
        ax.yaxis.set_major_locator(MaxNLocator(integer=True, prune='both'))
        ax.tick_params(axis='both', which='major', labelsize=12)
        ax.set_xlabel("Time (s)", fontsize=12)
        ax.set_ylabel(species, fontsize=12)
        ax.set_aspect('auto')

    for ax in axes[len(statistics):]:
        fig.delaxes(ax)

    plt.tight_layout()
    plt.show()
    return statistics


def generate_diffusion_trajectories_and_plot(num_trajectories, T, Omega, parameter_values, initial_conditions,
                                             stoichiometry_matrix, propensities):
    S = np.array(stoichiometry_matrix)
    time_points = np.linspace(0, T, 1001)
    dt = time_points[1] - time_points[0]
    Z0 = np.array(list(initial_conditions.values())).reshape(-1, 1)
    all_trajectories = np.zeros((num_trajectories, len(time_points), Z0.shape[0]))
    D = DiagonalMatrix(propensities, parameter_values, initial_conditions.keys())

    for traj in range(num_trajectories):
        Z_t = Z0
        Z_over_time = np.zeros((len(time_points), Z0.shape[0]))
        Z_over_time[0, :] = Z0.flatten()

        for i in range(1, len(time_points)):
            species_values = {name: Z_t[j, 0] for j, name in enumerate(initial_conditions.keys())}
            propensities_vector_t = D.evaluate_propensities(species_values)
            D_t = D.evaluate(Omega, species_values)
            dW_t = np.random.normal(0, np.sqrt(dt), size=(D_t.shape[0], 1))
            Z_t = Z_t + S @ propensities_vector_t * dt + S @ D_t @ dW_t
            Z_over_time[i, :] = Z_t.flatten()

        all_trajectories[traj, :, :] = Z_over_time

    # Compute statistics
    mean_trajectory = np.mean(all_trajectories, axis=0)
    variance_trajectory = np.var(all_trajectories, axis=0)

    # Prepare statistics for output
    statistics = {}
    for i, species in enumerate(initial_conditions.keys()):
        statistics[species] = {'mean': mean_trajectory[:, i], 'variance': variance_trajectory[:, i]}

    # Plot statistics
    num_species = Z0.shape[0]
    ncols = 2
    nrows = (num_species + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(15, nrows * 7.5), constrained_layout=True)
    axes = np.atleast_1d(axes).flatten()
    colors = plt.get_cmap('tab10').colors

    for idx in range(num_species):
        ax = axes[idx]
        mean_states = mean_trajectory[:, idx]
        variance_states = variance_trajectory[:, idx]
        ax.plot(time_points, mean_states, label=f'Species {idx + 1} mean', linewidth=2, color=colors[idx % len(colors)])
        ax.fill_between(time_points, mean_states - np.sqrt(variance_states), mean_states + np.sqrt(variance_states),
                        color=colors[idx % len(colors)], alpha=0.3, label=f'Species {idx + 1} variance')
        ax.legend(fontsize=12)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=10, prune='both'))
        ax.yaxis.set_major_locator(MaxNLocator(integer=True, prune='both'))
        ax.tick_params(axis='both', which='major', labelsize=12)
        ax.set_xlabel("Time (s)", fontsize=12)
        ax.set_ylabel(f'Species {idx + 1}', fontsize=12)
        ax.set_aspect('auto')

    for ax in axes[len(mean_trajectory[0]):]:
        fig.delaxes(ax)

    plt.tight_layout()
    plt.show()

    return all_trajectories, statistics


def compare_statistics(ssa_statistics, diffusion_statistics, species_names, T, Omega):
    common_time_points = np.linspace(0, T, 1001)
    n_species = len(species_names)
    ncols = 2
    nrows = (n_species + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(28, nrows * 18), constrained_layout=True)
    axes = np.atleast_1d(axes).flatten()
    colors = plt.get_cmap('tab10').colors

    def custom_formatter(x, pos):
        if x == 0:
            return '0'
        return f'{x:g}'

    for idx, species in enumerate(species_names):
        ax = axes[idx]
        ssa_mean = ssa_statistics[species]['mean']
        ssa_variance = ssa_statistics[species]['variance']
        diffusion_mean = diffusion_statistics[species]['mean'] * Omega
        diffusion_variance = diffusion_statistics[species]['variance'] * (Omega ** 2)
        ax.plot(common_time_points, ssa_mean, label=r'SSA ${}$ mean'.format(species), linewidth=2,
                color=colors[2 * idx % len(colors)])
        ax.fill_between(common_time_points, ssa_mean - np.sqrt(ssa_variance), ssa_mean + np.sqrt(ssa_variance),
                        color=colors[2 * idx % len(colors)], alpha=0.3, label=r'SSA ${}$ variance'.format(species))

        ax.plot(common_time_points, diffusion_mean, label=r'Diffusion ${}$ mean'.format(species), linewidth=2, linestyle='--',
                color=colors[(2 * idx + 1) % len(colors)])
        ax.fill_between(common_time_points, diffusion_mean - np.sqrt(diffusion_variance),
                        diffusion_mean + np.sqrt(diffusion_variance), color=colors[(2 * idx + 1) % len(colors)],
                        alpha=0.3, label=r'Diffusion ${}$ variance'.format(species))

        ax.legend(fontsize=45)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=10, prune='both'))
        ax.xaxis.set_major_formatter(FuncFormatter(custom_formatter))
        ax.yaxis.set_major_locator(MaxNLocator(integer=True, prune='both'))
        ax.yaxis.set_major_formatter(FuncFormatter(custom_formatter))
        ax.tick_params(axis='both', which='major', labelsize=50)
        ax.set_aspect('auto')

    for ax in axes[n_species:]:
        fig.delaxes(ax)

    plt.show()
    return fig


class DiagonalMatrix:
    def __init__(self, propensities, parameter_values, species_names):
        self.propensities = propensities
        self.parameter_values = parameter_values
        self.species_names = species_names

    def evaluate(self, omega, species_values):
        diag_matrix = np.zeros((len(self.propensities), len(self.propensities)))
        for i, propensity in enumerate(self.propensities):
            sig = inspect.signature(propensity)
            prop_args = [
                self.parameter_values[param.name] if param.name in self.parameter_values else species_values[param.name]
                for param in sig.parameters.values()
            ]
            diag_matrix[i, i] = np.sqrt(np.abs(propensity(*prop_args) / omega))
        return diag_matrix

    def evaluate_propensities(self, species_values):
        propensities_vector = np.zeros((len(self.propensities), 1))
        for i, propensity in enumerate(self.propensities):
            sig = inspect.signature(propensity)
            prop_args = [
                self.parameter_values[param.name] if param.name in self.parameter_values else species_values[param.name]
                for param in sig.parameters.values()
            ]
            propensities_vector[i, 0] = propensity(*prop_args)
        return propensities_vector

    def get_symbolic_jacobian(self):
        species_symbols = sp.symbols(self.species_names)
        parameter_symbols = {name: sp.symbols(name) for name in self.parameter_values.keys()}
        propensities_sympy = []
        for propensity in self.propensities:
            sig = inspect.signature(propensity)
            prop_args = []
            for param in sig.parameters.values():
                if param.name in parameter_symbols:
                    prop_args.append(parameter_symbols[param.name])
                elif param.name in self.species_names:
                    prop_args.append(species_symbols[self.species_names.index(param.name)])
                else:
                    raise ValueError(f"Unknown parameter or species: {param.name}")
            propensities_sympy.append(propensity(*prop_args))
        jacobian_matrix_sym = sp.Matrix(propensities_sympy).jacobian(species_symbols)
        return jacobian_matrix_sym

    def evaluate_jacobian(self, species_values):
        jacobian_matrix_sym = self.get_symbolic_jacobian()
        parameter_symbols = {name: sp.symbols(name) for name in self.parameter_values.keys()}
        species_symbols = sp.symbols(self.species_names)
        jacobian_matrix_func = sp.lambdify((list(parameter_symbols.values()) + list(species_symbols)), jacobian_matrix_sym, 'numpy')
        param_values = np.array(list(self.parameter_values.values()))
        combined_values = np.concatenate((param_values, [species_values[name] for name in self.species_names]))
        jacobian_matrix_eval = jacobian_matrix_func(*combined_values)
        return jacobian_matrix_eval

    def evaluate_gradient(self, species_values):
        species_symbols = sp.symbols(self.species_names)
        parameter_symbols = {name: sp.symbols(name) for name in self.parameter_values.keys()}
        gradients_sympy = []
        for propensity in self.propensities:
            sig = inspect.signature(propensity)
            prop_args = []
            for param in sig.parameters.values():
                if param.name in parameter_symbols:
                    prop_args.append(parameter_symbols[param.name])
                elif param.name in self.species_names:
                    prop_args.append(species_symbols[self.species_names.index(param.name)])
                else:
                    raise ValueError(f"Unknown parameter or species: {param.name}")
            prop_expr = propensity(*prop_args)
            gradient_sym = sp.Matrix([prop_expr]).jacobian(species_symbols)
            gradients_sympy.append(gradient_sym)
        gradients_matrix_func = sp.lambdify((list(parameter_symbols.values()) + list(species_symbols)), gradients_sympy, 'numpy')
        param_values = np.array(list(self.parameter_values.values()))
        combined_values = np.concatenate((param_values, [species_values[name] for name in self.species_names]))
        gradients_matrix_eval = gradients_matrix_func(*combined_values)
        return np.array(gradients_matrix_eval)

class ObservationModel:
    def __init__(self, observation_func, parameter_values, species_names):
        self.observation_func = observation_func
        self.parameter_values = parameter_values
        self.species_names = species_names

    def evaluate(self, species_values):
        all_values = {**self.parameter_values, **species_values}
        obs_args = [all_values[param.name] for param in inspect.signature(self.observation_func).parameters.values()]
        obs_value = self.observation_func(*obs_args)
        return np.array([[obs_value]])

    def get_symbolic_gradient(self):
        species_symbols = sp.symbols(self.species_names)
        parameter_symbols = {name: sp.symbols(name) for name in self.parameter_values.keys()}
        sig = inspect.signature(self.observation_func)
        prop_args = []
        for param in sig.parameters.values():
            if param.name in self.species_names:
                prop_args.append(species_symbols[self.species_names.index(param.name)])
            else:
                prop_args.append(parameter_symbols[param.name])
        obs_expr = self.observation_func(*prop_args)
        gradient_matrix_sym = sp.Matrix([obs_expr]).jacobian(species_symbols)
        return gradient_matrix_sym

    def evaluate_gradient(self, species_values):
        gradient_matrix_sym = self.get_symbolic_gradient()
        species_symbols = sp.symbols(self.species_names)
        gradient_matrix_func = sp.lambdify(species_symbols, gradient_matrix_sym, 'numpy')
        gradient_matrix_eval = gradient_matrix_func(*[species_values[name] for name in self.species_names])
        return np.array(gradient_matrix_eval)

def kalman_filter(initial_conditions, observations, parameter_values, propensities, stoichiometric_matrix, observation_func, Omega, rho, N):
    # Define common time points
    common_time_points = np.linspace(0, observations[0]['time'][-1], N)

    # Initial conditions
    hat_z_0 = np.array(initial_conditions['hat_z_0']).reshape(-1, 1)
    R_0 = initial_conditions['R_0']
    interpolated_observations = []

    # Interpolate observations and stack them
    for obs in observations:
        interp_obs = interp1d(obs['time'], obs['state'], axis=0, kind='linear', fill_value='extrapolate')
        interpolated_observations.append(interp_obs(common_time_points))
    interpolated_observations = np.stack(interpolated_observations, axis=0)




    # Initial trajectories
    hat_z_t = hat_z_0
    R_t = R_0
    hat_z_t_trajectory = [hat_z_t]
    R_t_trajectory = [R_t]
    S = np.array(stoichiometric_matrix)
    species_names = list(initial_conditions['species_values'].keys())
    species_values = initial_conditions['species_values']
    bar_z_t = np.array(list(species_values.values())).reshape(-1, 1)
    bar_z_t_trajectory = [bar_z_t]
    D_matrix = DiagonalMatrix(propensities, parameter_values, species_names)
    h_model = ObservationModel(observation_func, parameter_values, species_names)

    # Loop through time points
    for i in range(len(common_time_points) - 1):
        dt = common_time_points[i + 1] - common_time_points[i]

        # Update species values for current time step
        for j, name in enumerate(species_names):
            species_values[name] = bar_z_t[j, 0]

        # Calculate propensities and their Jacobian for current \bar{z}_t

        propensities_vector_bar_t = D_matrix.evaluate_propensities(species_values)


        # Evaluate Jacobian and other necessary terms for current \bar{z}_t
        J_t = D_matrix.evaluate_jacobian(species_values)
        D_bar_t = D_matrix.evaluate(Omega, species_values)
        h_bar = h_model.evaluate(species_values) / rho
        h_prime_bar = h_model.evaluate_gradient(species_values) / rho



        dbar_z_t = S @ propensities_vector_bar_t * dt
        bar_z_t = bar_z_t + dbar_z_t
        bar_z_t_trajectory.append(bar_z_t)

        # Innovation term using current observations and predictions
        innovation = (interpolated_observations[:, i].reshape(-1, 1) / (Omega * rho)) * dt - ((h_prime_bar @ (hat_z_t - bar_z_t) + h_bar) * dt)
        # print("innovation is ", innovation)
        # print("first term in the innovation is", (interpolated_observations[:, i].reshape(-1, 1) / (Omega * rho)) * dt)
        # print("second term in the innovation is", ((h_prime_bar @ (hat_z_t - bar_z_t) + h_bar) * dt))
        # Update \hat{z}_t using current information
        d_hat_z_t = (S @ J_t @ hat_z_t + S @ propensities_vector_bar_t - S @ J_t @ bar_z_t) * dt + R_t @ h_prime_bar.T @ innovation
        hat_z_t = hat_z_t + d_hat_z_t
        hat_z_t_trajectory.append(hat_z_t)

        # Update \mathbb{R}_t using current information
        dR_t = (S @ J_t @ R_t + R_t @ J_t.T @ S.T + S @ D_bar_t @ D_bar_t.T @ S.T - R_t @ h_prime_bar.T @ h_prime_bar @ R_t) * dt
        R_t = R_t + dR_t
        R_t_trajectory.append(R_t)


    # Scale trajectories by Omega and Omega^2
    hat_z_t_trajectory = np.hstack(hat_z_t_trajectory)
    R_t_trajectory = np.array(R_t_trajectory)
    hat_z_t_trajectory *= Omega
    R_t_trajectory *= Omega**2

    return common_time_points, hat_z_t_trajectory, R_t_trajectory, bar_z_t_trajectory


def plot_filters_results(common_time_points, hat_z_t_trajectory, R_t_trajectory, time_list, true_trajectories, species_names):
    num_species = hat_z_t_trajectory.shape[0]
    fig, axes = plt.subplots(1, num_species, figsize=(6 * num_species, 6), constrained_layout=True)
    for i in range(num_species):
        ax = axes[i] if num_species > 1 else axes
        ax.errorbar(common_time_points, hat_z_t_trajectory[i, :], yerr=np.sqrt(R_t_trajectory[:, i, i]), label=f'Estimated {species_names[i]}', fmt='-')
        ax.step(time_list, true_trajectories[i], label=f'True {species_names[i]} (SSA)', linestyle='--', linewidth=2, where='post')
        ax.legend(fontsize=15)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=10, prune='both'))
        ax.yaxis.set_major_locator(MaxNLocator(integer=True, prune='both'))
        ax.tick_params(axis='both', which='major', labelsize=15)
        ax.set_xlabel("Time (s)", fontsize=15)
        ax.set_ylabel(species_names[i], fontsize=15)
    plt.show()


def extended_kalman_filter(initial_conditions, observations, parameter_values, propensities, stoichiometric_matrix, observation_func, Omega, rho, N):
    # Define common time points
    common_time_points = np.linspace(0, observations[0]['time'][-1], N)

    # Initial conditions
    hat_z_0 = np.array(initial_conditions['hat_z_0']).reshape(-1, 1)
    R_0 = initial_conditions['R_0']
    interpolated_observations = []

    # Interpolate observations and stack them
    for obs in observations:
        interp_obs = interp1d(obs['time'], obs['state'], axis=0, kind='linear', fill_value='extrapolate')
        interpolated_observations.append(interp_obs(common_time_points))
    interpolated_observations = np.stack(interpolated_observations, axis=0)




    # Initial trajectories
    hat_z_t = hat_z_0
    R_t = R_0
    hat_z_t_trajectory = [hat_z_t]
    R_t_trajectory = [R_t]
    S = np.array(stoichiometric_matrix)
    species_names = list(initial_conditions['species_values'].keys())
    species_values = initial_conditions['species_values']
    D_matrix = DiagonalMatrix(propensities, parameter_values, species_names)
    h_model = ObservationModel(observation_func, parameter_values, species_names)

    # Loop through time points
    for i in range(len(common_time_points) - 1):
        dt = common_time_points[i + 1] - common_time_points[i]

        # Update species values for current time step
        for j, name in enumerate(species_names):
            species_values[name] = hat_z_t[j, 0]

        # Calculate propensities and their Jacobian for current \hat{z}_t

        propensities_vector_bar_t = D_matrix.evaluate_propensities(species_values)


        # Evaluate Jacobian and other necessary terms for current \hat{z}_t
        J_t = D_matrix.evaluate_jacobian(species_values)
        D_bar_t = D_matrix.evaluate(Omega, species_values)
        h_hat = h_model.evaluate(species_values) / rho
        h_prime_hat = h_model.evaluate_gradient(species_values) / rho


        # Innovation term using current observations and predictions
        innovation = ((interpolated_observations[:, i].reshape(-1, 1) / (Omega * rho))-h_hat) * dt

        # Update \hat{z}_t using current information
        d_hat_z_t = S @ propensities_vector_bar_t * dt + R_t @ h_prime_hat.T @ innovation
        hat_z_t = hat_z_t + d_hat_z_t
        hat_z_t_trajectory.append(hat_z_t)

        # Update \mathbb{R}_t using current information
        dR_t = (S @ J_t @ R_t + R_t @ J_t.T @ S.T + S @ D_bar_t @ D_bar_t.T @ S.T - R_t @ h_prime_hat.T @ h_prime_hat @ R_t) * dt
        R_t = R_t + dR_t
        R_t_trajectory.append(R_t)


    # Scale trajectories by Omega and Omega^2
    hat_z_t_trajectory = np.hstack(hat_z_t_trajectory)
    R_t_trajectory = np.array(R_t_trajectory)
    hat_z_t_trajectory *= Omega
    R_t_trajectory *= Omega**2

    return common_time_points, hat_z_t_trajectory, R_t_trajectory

# Custom formatter for the x-axis to remove trailing zeros after the decimal point

def custom_formatter(x, pos):
    if x == int(x):
        return f'{int(x)}'
    else:
        return f'{x:.1f}'


def plot_observation_trajectory(time_points, observed_trajectories, observable_species):
    fig, ax = plt.subplots(figsize=(10, 10), dpi=300, constrained_layout=True)

    for i, species in enumerate(observable_species):
        ax.step(time_points, observed_trajectories[:, i], label=r'${}$ (SSA)'.format(species), linestyle='--', linewidth=3, where='post')

    ax.legend(fontsize=45)
    # ax.set_title(r'Observation Trajectory', fontsize=24)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=10, prune='both'))
    ax.xaxis.set_major_formatter(FuncFormatter(custom_formatter))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True, prune='both'))
    ax.tick_params(axis='both', which='major', labelsize=50)
    #ax.set_xlabel(r'Time (s)', fontsize=30)
    # ax.set_ylabel(r'${}$'.format(observable_species[0]), fontsize=24)  # Update as needed
    ax.set_aspect(aspect='auto', adjustable='box')

    plt.show()
    return fig


def plot_filters_results_total(common_time_points, kalman_trajectory, kalman_R_trajectory, ekf_trajectory, ekf_R_trajectory, particle_trajectory, particle_sd, ffsp_trajectory, ffsp_sd, time_list, true_trajectories, observed_trajectories, time_list_1, hidden_species, point_density=1):
    num_hidden_species = len(hidden_species)
    fig, axes = plt.subplots(num_hidden_species, 4, figsize=(40, 10 * num_hidden_species), constrained_layout=True)

    filter_colors = {
        'Kalman': 'blue',
        'Extended Kalman': 'green',
        'Bootstrap Particle Filter': 'red',
        'FFSP': 'purple'
    }

    # Adjust the number of points to plot based on the point_density parameter
    num_points = int(len(common_time_points) * point_density)
    indices = np.linspace(0, len(common_time_points) - 1, num_points).astype(int)
    common_time_points_dense = common_time_points[indices]

    def interpolate_to_dense(time_points, data):
        interp_func = interp1d(time_points, data, axis=0, fill_value="extrapolate")
        return interp_func(common_time_points_dense)

    # Interpolating all trajectories and standard deviations
    kalman_trajectory_dense = interpolate_to_dense(common_time_points, kalman_trajectory.T).T
    ekf_trajectory_dense = interpolate_to_dense(common_time_points, ekf_trajectory.T).T
    particle_trajectory_dense = interpolate_to_dense(common_time_points, particle_trajectory)
    ffsp_trajectory_dense = interpolate_to_dense(common_time_points, ffsp_trajectory)
    particle_sd_dense = interpolate_to_dense(common_time_points, particle_sd)
    ffsp_sd_dense = interpolate_to_dense(common_time_points, ffsp_sd)

    # Custom formatter for the x-axis
    def custom_formatter(x, pos):
        if x == int(x):
            return f'{int(x)}'
        else:
            return f'{x:.1f}'

    for i in range(num_hidden_species):
        ssa_color = 'orange' if i == 0 else 'brown'

        # Kalman Filter
        ax_kalman = axes[i, 0] if num_hidden_species > 1 else axes[0]
        ax_kalman.errorbar(common_time_points_dense, kalman_trajectory_dense[i, :], yerr=np.sqrt(np.maximum(kalman_R_trajectory[indices, i, i], 0)), label=r'Estimated ${}$'.format(hidden_species[i]), fmt='-', color=filter_colors['Kalman'])
        ax_kalman.step(time_list, true_trajectories[i, :], label=r'True ${}$ (SSA)'.format(hidden_species[i]), linestyle='--', linewidth=2, where='post', color=ssa_color)
        ax_kalman.legend(fontsize=30)
        ax_kalman.set_title(r'Kalman Filter $' + f'{hidden_species[i]}' + r'$', fontsize=34)
        ax_kalman.xaxis.set_major_locator(MaxNLocator(nbins=10, prune='both'))
        ax_kalman.xaxis.set_major_formatter(FuncFormatter(custom_formatter))
        ax_kalman.yaxis.set_major_locator(MaxNLocator(integer=True, prune='both'))
        ax_kalman.tick_params(axis='both', which='major', labelsize=30)
        ax_kalman.set_xlabel(r'Time (s)', fontsize=34)
        ax_kalman.set_ylabel(r'${}$'.format(hidden_species[i]), fontsize=34)

        # Extended Kalman Filter
        ax_ekf = axes[i, 1] if num_hidden_species > 1 else axes[1]
        ax_ekf.errorbar(common_time_points_dense, ekf_trajectory_dense[i, :], yerr=np.sqrt(np.maximum(ekf_R_trajectory[indices, i, i], 0)), label=r'Estimated ${}$'.format(hidden_species[i]), fmt='-', color=filter_colors['Extended Kalman'])
        ax_ekf.step(time_list, true_trajectories[i, :], label=r'True ${}$ (SSA)'.format(hidden_species[i]), linestyle='--', linewidth=2, where='post', color=ssa_color)
        ax_ekf.legend(fontsize=30)
        ax_ekf.set_title(r'Extended Kalman Filter $' + f'{hidden_species[i]}' + r'$', fontsize=34)
        ax_ekf.xaxis.set_major_locator(MaxNLocator(nbins=10, prune='both'))
        ax_ekf.xaxis.set_major_formatter(FuncFormatter(custom_formatter))
        ax_ekf.yaxis.set_major_locator(MaxNLocator(integer=True, prune='both'))
        ax_ekf.tick_params(axis='both', which='major', labelsize=30)
        ax_ekf.set_xlabel(r'Time (s)', fontsize=34)
        ax_ekf.set_ylabel(r'${}$'.format(hidden_species[i]), fontsize=34)

        # Bootstrap Particle Filter
        ax_particle = axes[i, 2] if num_hidden_species > 1 else axes[2]
        ax_particle.errorbar(common_time_points_dense, particle_trajectory_dense[:, i], yerr=particle_sd_dense[:, i], label=r'Estimated ${}$'.format(hidden_species[i]), fmt='-', color=filter_colors['Bootstrap Particle Filter'])
        ax_particle.step(time_list, true_trajectories[i, :], label=r'True ${}$ (SSA)'.format(hidden_species[i]), linestyle='--', linewidth=2, where='post', color=ssa_color)
        ax_particle.legend(fontsize=30)
        ax_particle.set_title(r'BPF $' + f'{hidden_species[i]}' + r'$', fontsize=34)
        ax_particle.xaxis.set_major_locator(MaxNLocator(nbins=10, prune='both'))
        ax_particle.xaxis.set_major_formatter(FuncFormatter(custom_formatter))
        ax_particle.yaxis.set_major_locator(MaxNLocator(integer=True, prune='both'))
        ax_particle.tick_params(axis='both', which='major', labelsize=30)
        ax_particle.set_xlabel(r'Time (s)', fontsize=34)
        ax_particle.set_ylabel(r'${}$'.format(hidden_species[i]), fontsize=34)

        # FFSP
        ax_ffsp = axes[i, 3] if num_hidden_species > 1 else axes[3]
        ax_ffsp.errorbar(common_time_points_dense, ffsp_trajectory_dense[:, i], yerr=ffsp_sd_dense[:, i], label=r'Estimated ${}$'.format(hidden_species[i]), fmt='-', color=filter_colors['FFSP'])
        ax_ffsp.step(time_list, true_trajectories[i, :], label=r'True ${}$ (SSA)'.format(hidden_species[i]), linestyle='--', linewidth=2, where='post', color=ssa_color)
        ax_ffsp.legend(fontsize=30)
        ax_ffsp.set_title(r'FFSP $' + f'{hidden_species[i]}' + r'$', fontsize=34)
        ax_ffsp.xaxis.set_major_locator(MaxNLocator(nbins=10, prune='both'))
        ax_ffsp.xaxis.set_major_formatter(FuncFormatter(custom_formatter))
        ax_ffsp.yaxis.set_major_locator(MaxNLocator(integer=True, prune='both'))
        ax_ffsp.tick_params(axis='both', which='major', labelsize=30)
        ax_ffsp.set_xlabel(r'Time (s)', fontsize=34)
        ax_ffsp.set_ylabel(r'${}$'.format(hidden_species[i]), fontsize=34)

    plt.show()
    return fig


def plot_filters_results_total_no_labels(common_time_points, kalman_trajectory, kalman_R_trajectory, ekf_trajectory, ekf_R_trajectory, particle_trajectory, particle_sd, ffsp_trajectory, ffsp_sd, time_list, true_trajectories, observed_trajectories, time_list_1, hidden_species, point_density=1):
    num_hidden_species = len(hidden_species)
    fig, axes = plt.subplots(num_hidden_species, 4, figsize=(40, 10 * num_hidden_species), constrained_layout=True)

    filter_colors = {
        'Kalman': 'blue',
        'Extended Kalman': 'green',
        'Bootstrap Particle Filter': 'red',
        'FFSP': 'purple'
    }

    # Adjust the number of points to plot based on the point_density parameter
    num_points = int(len(common_time_points) * point_density)
    indices = np.linspace(0, len(common_time_points) - 1, num_points).astype(int)
    common_time_points_dense = common_time_points[indices]

    def interpolate_to_dense(time_points, data):
        interp_func = interp1d(time_points, data, axis=0, fill_value="extrapolate")
        return interp_func(common_time_points_dense)

    # Interpolating all trajectories and standard deviations
    kalman_trajectory_dense = interpolate_to_dense(common_time_points, kalman_trajectory.T).T
    ekf_trajectory_dense = interpolate_to_dense(common_time_points, ekf_trajectory.T).T
    particle_trajectory_dense = interpolate_to_dense(common_time_points, particle_trajectory)
    ffsp_trajectory_dense = interpolate_to_dense(common_time_points, ffsp_trajectory)
    particle_sd_dense = interpolate_to_dense(common_time_points, particle_sd)
    ffsp_sd_dense = interpolate_to_dense(common_time_points, ffsp_sd)

    # Custom formatter for the x-axis
    def custom_formatter(x, pos):
        if x == int(x):
            return f'{int(x)}'
        else:
            return f'{x:.1f}'

    for i in range(num_hidden_species):
        ssa_color = 'orange' if i == 0 else 'brown'

        # Kalman Filter
        ax_kalman = axes[i, 0] if num_hidden_species > 1 else axes[0]
        ax_kalman.errorbar(common_time_points_dense, kalman_trajectory_dense[i, :], yerr=np.sqrt(np.maximum(kalman_R_trajectory[indices, i, i], 0)), label=r'Estimated ${}$'.format(hidden_species[i]), fmt='-', color=filter_colors['Kalman'])
        ax_kalman.step(time_list, true_trajectories[i, :], label=r'True ${}$ (SSA)'.format(hidden_species[i]), linestyle='--', linewidth=3.5, where='post', color=ssa_color)
        ax_kalman.legend(fontsize=45)  # Increased font size
        ax_kalman.xaxis.set_major_locator(MaxNLocator(nbins=10, prune='both'))
        ax_kalman.xaxis.set_major_formatter(FuncFormatter(custom_formatter))
        ax_kalman.yaxis.set_major_locator(MaxNLocator(integer=True, prune='both'))
        ax_kalman.tick_params(axis='both', which='major', labelsize=50)  # Increased tick font size

        # Extended Kalman Filter
        ax_ekf = axes[i, 1] if num_hidden_species > 1 else axes[1]
        ax_ekf.errorbar(common_time_points_dense, ekf_trajectory_dense[i, :], yerr=np.sqrt(np.maximum(ekf_R_trajectory[indices, i, i], 0)), label=r'Estimated ${}$'.format(hidden_species[i]), fmt='-', color=filter_colors['Extended Kalman'])
        ax_ekf.step(time_list, true_trajectories[i, :], label=r'True ${}$ (SSA)'.format(hidden_species[i]), linestyle='--', linewidth=3.5, where='post', color=ssa_color)
        ax_ekf.legend(fontsize=45)  # Increased font size
        ax_ekf.xaxis.set_major_locator(MaxNLocator(nbins=10, prune='both'))
        ax_ekf.xaxis.set_major_formatter(FuncFormatter(custom_formatter))
        ax_ekf.yaxis.set_major_locator(MaxNLocator(integer=True, prune='both'))
        ax_ekf.tick_params(axis='both', which='major', labelsize=50)  # Increased tick font size

        # Bootstrap Particle Filter
        ax_particle = axes[i, 2] if num_hidden_species > 1 else axes[2]
        ax_particle.errorbar(common_time_points_dense, particle_trajectory_dense[:, i], yerr=particle_sd_dense[:, i], label=r'Estimated ${}$'.format(hidden_species[i]), fmt='-', color=filter_colors['Bootstrap Particle Filter'])
        ax_particle.step(time_list, true_trajectories[i, :], label=r'True ${}$ (SSA)'.format(hidden_species[i]), linestyle='--', linewidth=3.5, where='post', color=ssa_color)
        ax_particle.legend(fontsize=45)  # Increased font size
        ax_particle.xaxis.set_major_locator(MaxNLocator(nbins=10, prune='both'))
        ax_particle.xaxis.set_major_formatter(FuncFormatter(custom_formatter))
        ax_particle.yaxis.set_major_locator(MaxNLocator(integer=True, prune='both'))
        ax_particle.tick_params(axis='both', which='major', labelsize=50)  # Increased tick font size

        # FFSP
        ax_ffsp = axes[i, 3] if num_hidden_species > 1 else axes[3]
        ax_ffsp.errorbar(common_time_points_dense, ffsp_trajectory_dense[:, i], yerr=ffsp_sd_dense[:, i], label=r'Estimated ${}$'.format(hidden_species[i]), fmt='-', color=filter_colors['FFSP'])
        ax_ffsp.step(time_list, true_trajectories[i, :], label=r'True ${}$ (SSA)'.format(hidden_species[i]), linestyle='--', linewidth=3.5, where='post', color=ssa_color)
        ax_ffsp.legend(fontsize=45)  # Increased font size
        ax_ffsp.xaxis.set_major_locator(MaxNLocator(nbins=10, prune='both'))
        ax_ffsp.xaxis.set_major_formatter(FuncFormatter(custom_formatter))
        ax_ffsp.yaxis.set_major_locator(MaxNLocator(integer=True, prune='both'))
        ax_ffsp.tick_params(axis='both', which='major', labelsize=50)  # Increased tick font size

    plt.show()
    return fig


