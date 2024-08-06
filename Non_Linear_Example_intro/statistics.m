%% Load data
ssa_stats = load('ssa_stats.mat');
diffusion_stats = load('diffusion_stats.mat');

%%  Species names
species_names = fieldnames(ssa_stats);

%% Time points
T = 1; %final time
Omega=10; %volume parameter 
common_time_points = linspace(0, T, 1001);

%%
n_species = length(species_names);

% Create the figure for visualization
figure('Position', [100, 100, 1400, 1980]);

colors = lines(10);

% Manually define positions to ensure alignment and equal size
positions = {
    [0.2, 0.55, 0.4, 0.4],  % Position for first subplot
    [0.45, 0.55, 0.4, 0.4], % Position for second subplot, closer to the first subplot
    [0.325, 0.05, 0.4, 0.4] % Position for third subplot centered in the second row
};

for idx = 1:n_species
    species = species_names{idx};
    
    % Create the axes with the specified position
    axes('Position', positions{idx});
    
    ssa_mean = ssa_stats.(species).mean;
    ssa_variance = ssa_stats.(species).variance;
    diffusion_mean = diffusion_stats.(species).mean * Omega;
    diffusion_variance = diffusion_stats.(species).variance * (Omega ^ 2);
    
    hold on;
    
    % Plot SSA mean and variance
    plot(common_time_points, ssa_mean, 'LineWidth', 2, 'Color', colors(mod(2*idx-2, size(colors, 1))+1, :), 'DisplayName', sprintf('SSA %s mean', species));
    fill([common_time_points, fliplr(common_time_points)], [ssa_mean - sqrt(ssa_variance), fliplr(ssa_mean + sqrt(ssa_variance))], ...
        colors(mod(2*idx-2, size(colors, 1))+1, :), 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', sprintf('SSA %s variance', species));
    
    % Plot Diffusion mean and variance
    plot(common_time_points, diffusion_mean, 'LineWidth', 2, 'LineStyle', '--', 'Color', colors(mod(2*idx-1, size(colors, 1))+1, :), 'DisplayName', sprintf('Diffusion %s mean', species));
    fill([common_time_points, fliplr(common_time_points)], [diffusion_mean - sqrt(diffusion_variance), fliplr(diffusion_mean + sqrt(diffusion_variance))], ...
        colors(mod(2*idx-1, size(colors, 1))+1, :), 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', sprintf('Diffusion %s variance', species));
    
    legend('show', 'FontSize', 12);
    set(gca, 'FontSize', 16); % Increase tick font size
    % Remove xlabel line
    ylabel(''); % Optionally add y-axis label if needed
    axis square; % Make subplots square
    hold off;
end

% Instructions for manual positioning of legends
disp('Position the legends as desired and press any key to save the figure.');

% Wait for a key press
pause;

% Capture the current figure as an image
frame = getframe(gcf);
im = frame2im(frame);

% Create a new figure sized for A3 paper
a3_width = 11.7; % A3 width in inches
a3_height = 16.5; % A3 height in inches
new_fig = figure('Units', 'inches', 'Position', [0, 0, a3_width, a3_height], 'PaperUnits', 'inches', 'PaperSize', [a3_width, a3_height]);

% Display the captured image in the new figure
imshow(im, 'Border', 'tight');

% Save the new figure as a PDF in A3 portrait format
saveas(new_fig, 'figure_a3_portrait.pdf');
