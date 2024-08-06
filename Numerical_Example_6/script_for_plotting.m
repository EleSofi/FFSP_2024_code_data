


%%
fluo = 73;
figure_positions = [1, 5, 11]; % The specific positions to iterate over

for m = figure_positions
    X_1 = relevantCellTrajectory{1, m}.spotIntensity'./fluo;
    t_ob = data_to_save{m}.t_ob;
    Y = data_to_save{m}.Y;
    E_FSP = data_to_save{m}.E_FSP;
    SD_FSP = data_to_save{m}.SD_FSP;

    f = figure('Units', 'normalized', 'OuterPosition', [0 0 1 1], 'Name', sprintf('Filtering_%d_cell', m)); % Maximize the figure window

    freq = 1;

    subplot(2,1,1)
    stairs(t_ob, Y(1,:), 'm', 'LineWidth', 2)
    xlim([0 t_ob(end)])
    xticks(0:50:max(t_ob)) % Set x-ticks every 50 units
    xlabel("t (min)")
    ylabel('Y(t)')
    title('Observation Process')
    set(gca, 'FontSize', 20)
    set(gca, 'Position', [0.1, 0.58, 0.85, 0.35]) % Adjust position and size of the subplot
    axis square

    subplot(2,1,2)
    errorbar(t_ob([1 [1:freq:end]]), E_FSP([1 [1:freq:end]], 3)', SD_FSP([1 [1:freq:end]], 3), 'r--', 'LineWidth', 3)
    hold on
    stairs(timeAxis, X_1, 'b', 'LineWidth', 2)
    hold off
    title('Hidden Process (mRNA)')
    xlabel("t (min)")
    ylabel('Molecular Counts')
    xlim([0 t_ob(end)])
    xticks(0:50:max(t_ob)) % Set x-ticks every 50 units
    legend('FFSP','Exact Trajectory', 'Location', 'northeast') % Initial legend position

    set(gca, 'FontSize', 20)
    set(gca, 'Position', [0.1, 0.12, 0.85, 0.35]) % Adjust position and size of the subplot
    axis square

    % Create a save button
    uicontrol('Style', 'pushbutton', 'String', 'Save Figure', ...
              'Units', 'normalized', 'Position', [0.45, 0.02, 0.1, 0.05], ...
              'Callback', {@saveFigure, m});
end

function saveFigure(~, ~, m)
    f = gcf;
    % Save the adjusted figure as a .tif file with tight cropping
    exportgraphics(f, sprintf('Filtering_%d_cell_tight.tif', m), 'Resolution', 300, 'ContentType', 'image', 'BackgroundColor', 'white');
    disp(['Figure ' num2str(m) ' saved successfully.']);
end