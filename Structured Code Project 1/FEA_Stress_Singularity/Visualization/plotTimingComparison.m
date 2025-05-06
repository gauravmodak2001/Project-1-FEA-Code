function plotTimingComparison(elementCounts4Node, elementCounts8Node, cpuTimes4Node, cpuTimes8Node)
    % Create figure for timing comparison
    figure;
    loglog(elementCounts4Node, cpuTimes4Node, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
    hold on;
    loglog(elementCounts8Node, cpuTimes8Node, 'rs-', 'LineWidth', 2, 'MarkerSize', 8);
    grid on;
    
    % Add labels and legend
    xlabel('Number of Elements', 'FontSize', 12);
    ylabel('CPU Time (seconds)', 'FontSize', 12);
    title('CPU Time vs Number of Elements', 'FontSize', 14);
    legend('4-node Linear Elements', '8-node Quadratic Elements', 'Location', 'northwest');
    
    % Fit power law to data
    coeffs_4node = polyfit(log10(elementCounts4Node), log10(cpuTimes4Node), 1);
    coeffs_8node = polyfit(log10(elementCounts8Node), log10(cpuTimes8Node), 1);
    
    % Add text showing the power law relationships
    text(elementCounts4Node(2), cpuTimes4Node(2)*1.5, ...
        sprintf('Time \\propto n^{%.2f}', coeffs_4node(1)), ...
        'FontSize', 12, 'Color', 'blue');
    text(elementCounts8Node(2), cpuTimes8Node(2)*1.5, ...
        sprintf('Time \\propto n^{%.2f}', coeffs_8node(1)), ...
        'FontSize', 12, 'Color', 'red');
    
    % Save the figure
    saveas(gcf, 'timing_comparison.png');
    
    % Also create a table with the timing data
    fprintf('=== Timing Comparison ===\n');
    fprintf('Mesh Size | # Elements (4-node) | CPU Time (4-node) | # Elements (8-node) | CPU Time (8-node)\n');
    fprintf('---------+--------------------+------------------+--------------------+------------------\n');
    for i = 1:length(elementCounts4Node)
        fprintf('  %.4f  |        %4d        |      %.4f s     |        %4d        |      %.4f s\n', ...
            meshSizes(i), elementCounts4Node(i), cpuTimes4Node(i), ...
            elementCounts8Node(i), cpuTimes8Node(i));
    end
    
    % Calculate and display the computational efficiency ratio
    efficiency = zeros(length(elementCounts4Node), 1);
    for i = 1:length(elementCounts4Node)
        % Time per DOF ratio between 4-node and 8-node elements
        dof_4node = elementCounts4Node(i) * 4 * 2; % 4 nodes per element, 2 DOFs per node
        dof_8node = elementCounts8Node(i) * 8 * 2; % 8 nodes per element, 2 DOFs per node
        
        time_per_dof_4node = cpuTimes4Node(i) / dof_4node;
        time_per_dof_8node = cpuTimes8Node(i) / dof_8node;
        
        efficiency(i) = time_per_dof_4node / time_per_dof_8node;
    end
    
    fprintf('\n=== Computational Efficiency ===\n');
    fprintf('Element Size | Time per DOF ratio (4-node/8-node)\n');
    fprintf('------------+----------------------------------\n');
    for i = 1:length(elementCounts4Node)
        fprintf('    %.4f    |              %.4f\n', meshSizes(i), efficiency(i));
    end
end