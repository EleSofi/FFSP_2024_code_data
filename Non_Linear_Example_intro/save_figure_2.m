% Plot observation trajectory and filter results
plot_filters_results(common_time_points, kalman_trajectory, kalman_R_trajectory, ...
    ekf_trajectory, ekf_R_trajectory, E_PF_interp, SD_PF_interp, E_FSP_interp, SD_FSP_interp, ...
    time_list, hidden_trajectories, hidden_species, point_density);

disp('Adjust the legends manually and press any key to save the figure.');
pause;

% Call the save function
saveFigure();

function saveFigure()
    % Get the current figure handle
    hFig = gcf;
    
    % Set printing options for A3 landscape format
    set(hFig, 'PaperOrientation', 'landscape');
    set(hFig, 'PaperType', 'a3');
    set(hFig, 'PaperUnits', 'centimeters');
    
    % Adjust PaperPosition for better margins
    leftMargin = 3.47; % Left margin
    bottomMargin = 0.80; % Bottom margin
    width = 50.80; % Width of the figure
    height = 28.08; % Height of the figure
    set(hFig, 'PaperPositionMode', 'manual');
    set(hFig, 'PaperPosition', [leftMargin, bottomMargin, width, height]);

    % Save the figure using exportgraphics
    exportgraphics(hFig, 'filter_results.pdf', 'ContentType', 'vector', 'Resolution', 300);
end

