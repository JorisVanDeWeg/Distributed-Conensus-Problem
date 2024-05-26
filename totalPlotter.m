function totalPlotter(Difference, alg_name, title_var)
    numIterations = size(Difference,1);
    num_alg = size(alg_name,2);
    colors = lines(num_alg); % Get unique colors for each algorithm

    % Plotting all differences on the same graph
    figure;
    hold on; % Hold on to add multiple plots to the same figure
    % Plot each algorithm's performance
    for i = 1:num_alg
        semilogy(1:numIterations, Difference(:, i), 'Color', colors(i, :), 'DisplayName', alg_name(i), 'LineWidth', 1.35);
    end
    
    jor = max(max(Difference));
    bur = ceil(log10(jor));
  
    % Adding plot details
    title('Performance per Iteration for ' + title_var);
    xlabel('Iteration');
    ylabel('MSE');
    legend show; 
    set(gca, 'YScale', 'log')
    ylim([(10^-12) (10^bur)])
    hold off;
end