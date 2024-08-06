% Plot observation trajectory and filter results
plot_filters_results_with_observation_1(common_time_points, kalman_trajectory, kalman_R_trajectory, ...
    ekf_trajectory, ekf_R_trajectory, E_PF_interp, SD_PF_interp, E_FSP_interp, SD_FSP_interp, ...
    time_list, hidden_trajectories, time_of_observations, state_sequence, hidden_species, observable_species, point_density);
% Add a button to save the figure
    uicontrol('Style', 'pushbutton', 'String', 'Save Figure', ...
              'Position', [10 10 100 40], 'Callback', @saveFigure);

    function saveFigure(~, ~)
        % Set printing options for A3 landscape format
        set(gcf, 'PaperOrientation', 'landscape');
        set(gcf, 'PaperType', 'a3');
        set(gcf, 'PaperPositionMode', 'auto');
        % Save the figure as a PDF
        print(gcf, 'filter_results.pdf', '-dpdf', '-r300');  % Save as PDF with 300 dpi resolution
    end

