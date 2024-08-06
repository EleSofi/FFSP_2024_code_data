function plot_filters_results(common_time_points, kalman_trajectory, kalman_R_trajectory, ekf_trajectory, ekf_R_trajectory, particle_trajectory, particle_sd, ffsp_trajectory, ffsp_sd, time_list, true_trajectories, hidden_species, point_density)
    num_hidden_species = length(hidden_species);
    figure('Position', [100, 100, 3200, 900 * num_hidden_species], 'PaperPositionMode', 'auto', 'Color', 'w');  % Adjusted figure size and set background to white

    kalman_color = "#4DBEEE";
    ekf_color = "#D95319";
    particle_filter_color = "#EDB120";
    ffsp_color = 'r';
    ssa_color = 'blue';  % Custom color for SSA trajectories

    % Adjust the number of points to plot based on the point_density parameter
    num_points = floor(length(common_time_points) * point_density);
    indices = round(linspace(1, length(common_time_points), num_points));

    % Subplot positions
    left_margin = 0.1;
    bottom_margin = 0.1;
    plot_width = 0.18;
    plot_height = 0.67 / num_hidden_species;  % Adjust for total number of rows
    horizontal_spacing = 0.02;
    vertical_spacing = 0.1;

    for i = 1:num_hidden_species
        for j = 1:4
            subplot_position = [left_margin + (j-1)*(plot_width + horizontal_spacing), ...
                                bottom_margin + (num_hidden_species-i)*(plot_height + vertical_spacing), ...
                                plot_width, plot_height];

            ax = subplot('Position', subplot_position);
            hold on;

            switch j
                case 1
                    errorbar(common_time_points(indices), kalman_trajectory(i, indices), sqrt(max(kalman_R_trajectory(indices, i, i), 0)), ...
                        'DisplayName', sprintf('Estimated %s', hidden_species{i}), 'Color', kalman_color, 'LineWidth', 1.5);  % Thicker error bars
                case 2
                    errorbar(common_time_points(indices), ekf_trajectory(i, indices), sqrt(max(ekf_R_trajectory(indices, i, i), 0)), ...
                        'DisplayName', sprintf('Estimated %s', hidden_species{i}), 'Color', ekf_color, 'LineWidth', 1.5);  % Thicker error bars
                case 3
                    errorbar(common_time_points(indices), particle_trajectory(indices, i), particle_sd(indices, i), ...
                        'DisplayName', sprintf('Estimated %s', hidden_species{i}), 'Color', particle_filter_color, 'LineWidth', 1.5);  % Thicker error bars
                case 4
                    errorbar(common_time_points(indices), ffsp_trajectory(indices, i), ffsp_sd(indices, i), ...
                        'DisplayName', sprintf('Estimated %s', hidden_species{i}), 'Color', ffsp_color, 'LineWidth', 1.5);  % Thicker error bars
            end

            stairs(time_list, true_trajectories(i, :), 'DisplayName', sprintf('True %s (SSA)', hidden_species{i}), ...
                'LineWidth', 3, 'LineStyle', '-', 'Color', ssa_color);  % Thicker line for true trajectories

            hold off;
            legend('show', 'FontSize', 15);  % Legend font size 15
            set(gca, 'FontSize', 20);  % Axis tick font size 20
            xlim([common_time_points(1), common_time_points(end)]);
            set(ax, 'PlotBoxAspectRatio', [1 1 1]);
        end
    end
end
