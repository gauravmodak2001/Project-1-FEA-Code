function plotMaxStressVsMeshSize(results4Node, results8Node, meshSizes)
    % Create figure for comparing maximum stress vs mesh size
    figure;
    
    % Get maximum stress values for 4-node elements
    maxStress4Node = zeros(length(meshSizes), 1);
    for i = 1:length(meshSizes)
        % Find maximum stress in all elements
        maxStress = -Inf;
        for j = 1:length(results4Node{i}.Displacementandstraincalc)
            elementStress = results4Node{i}.Displacementandstraincalc(j).Stress(:,:,2); % sigma_yy
            maxStress = max(maxStress, max(elementStress(:)));
        end
        maxStress4Node(i) = maxStress;
    end
    
    % Get maximum stress values for 8-node elements
    maxStress8Node = zeros(length(meshSizes), 1);
    for i = 1:length(meshSizes)
        % Find maximum stress in all elements
        maxStress = -Inf;
        for j = 1:length(results8Node{i}.Displacementandstraincalc)
            elementStress = results8Node{i}.Displacementandstraincalc(j).Stress(:,:,2); % sigma_yy
            maxStress = max(maxStress, max(elementStress(:)));
        end
        maxStress8Node(i) = maxStress;
    end
    
    % Plot on loglog scale
    loglog(meshSizes, maxStress4Node, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
    hold on;
    loglog(meshSizes, maxStress8Node, 'rs-', 'LineWidth', 2, 'MarkerSize', 8);
    
    % Add labels and legend
    xlabel('Element Size');
    ylabel('Maximum \sigma_{yy}');
    title('Maximum \sigma_{yy} vs Element Size');
    legend('4-node Linear Elements', '8-node Quadratic Elements', 'Location', 'best');
    grid on;
    
    % Save the figure
    saveas(gcf, 'max_stress_vs_mesh_size_comparison.png');
end