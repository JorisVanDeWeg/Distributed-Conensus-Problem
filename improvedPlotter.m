function improvedPlotter(Difference, alg_name, title_var)
    numIterations = size(Difference,1);
    num_alg = size(alg_name,2);
    colors = lines(num_alg); % Get unique colors for each algorithm

    % Define unique line styles for different algorithms
    lineStyles = {'-', '--', ':', '-.'};
    
    % Extract main algorithm names to assign line styles
    mainAlgNames = cell(1, num_alg);
    for i = 1:num_alg
        parts = strsplit(alg_name{i}, ' ');
        mainAlgNames{i} = parts{1};
    end

    uniqueAlgNames = unique(mainAlgNames);
    numUniqueAlgs = numel(uniqueAlgNames);
    % Map each main algorithm to a unique line style
    styleMap = containers.Map(uniqueAlgNames, lineStyles(1:numUniqueAlgs));

    % Plotting all differences on the same graph
    figure;
    hold on; % Hold on to add multiple plots to the same figure
    % Plot each algorithm's performance
    for i = 1:num_alg
        mainAlg = mainAlgNames{i};
        semilogy(1:numIterations, Difference(:, i), 'Color', colors(i, :), ...
                 'LineStyle', styleMap(mainAlg), 'DisplayName', alg_name(i), 'LineWidth', 1.35);
    end
    
    % Finding maximum limit for plotting based on the max of results
    jor = max(max(Difference));
    bur = ceil(log10(jor));

    % Adding plot details
    title(['Performance per Iteration for ', title_var]);
    xlabel('Iteration');
    ylabel('MSE');
    legend show;
    set(gca, 'YScale', 'log');
    ylim([(10^-12) (10^bur)]);
    hold off;
end