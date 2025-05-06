% Stress Singularity Analysis - Main Script
% Analyzes a 2D problem with 4-node and 8-node elements
% Problem: 2×6 rectangle with force at top-right, constrained at bottom and right

clear all; close all; clc;

% Add all folders to MATLAB path
addpath('Model Creation');
addpath('Element Functions');
addpath('Visualization');
addpath('Helper Functions');
mkdir('Results'); % Create results folder if it doesn't exist


% Define material properties and problem geometry
E = 70000;        % Young's modulus
v = 0.33;         % Poisson's ratio  
thickness = 0.3;  % Thickness
problemWidth = 2; % Width of rectangle (x-direction)
problemHeight = 6; % Height of rectangle (y-direction)
forceValue = 20;  % Concentrated force value

% Array of mesh sizes (each one half the previous)
meshSizes = [1, 0.5, 0.25, 0.125];

% Arrays to store results
results4Node = cell(length(meshSizes), 1);
results8Node = cell(length(meshSizes), 1);
cpuTimes4Node = zeros(length(meshSizes), 1);
cpuTimes8Node = zeros(length(meshSizes), 1);
elementCounts4Node = zeros(length(meshSizes), 1);
elementCounts8Node = zeros(length(meshSizes), 1);

% Part 1: Analysis with 4-node elements
fprintf('\n===== ANALYSIS WITH 4-NODE LINEAR ELEMENTS =====\n');

for i = 1:length(meshSizes)
    fprintf('\nProcessing 4-node mesh with size = %f\n', meshSizes(i));
    
    % Start timing
    tic;
    
    % Create model with 4-node elements
    model = create4NodeModel(meshSizes(i), problemWidth, problemHeight, E, v, thickness);
    
    % Visualize mesh
    figure;
    plotMesh(model);
    title(sprintf('4-Node Mesh with element size = %f', meshSizes(i)));
    saveas(gcf, sprintf('4node_mesh_size_%f.png', meshSizes(i)));
    
    % Apply concentrated force at top-right corner
    model = applyConcentratedForce(model, forceValue);
    
    % Create elements structure
    elements_structure = createElementsStructure(model);
    
    % Create global stiffness matrix and force vector
    Global_elements_structure = Global_Stifness_Matrix(@convert_to_index, elements_structure, model);
    
    % Apply boundary conditions
    [K, F] = assemblyBC(model, Global_elements_structure);
    
    % Solve the system (K*Q = F)
    Q = K \ F;
    Global_elements_structure.Q = Q;
    
    % Calculate displacements and stresses
    Displacementandstraincalc_results = Displacementandstraincalc(Global_elements_structure, model);
    
    % Plot sigma_yy field directly
    figure;
    hold on;
    title(sprintf('\\sigma_{yy} field with 4-node elements, size = %f', meshSizes(i)));
    xlabel('X-coordinate');
    ylabel('Y-coordinate');
    colormap(jet);
    
    % Find min/max stress values
    minStress = Inf;
    maxStress = -Inf;
    
    % Loop through elements to extract and plot stress data
    for elemIdx = 1:model.numElements
        % Get element coordinates
        x = Displacementandstraincalc_results(elemIdx).x;
        y = Displacementandstraincalc_results(elemIdx).y;
        
        % Get stress data - use component 2 for sigma_yy
        stressComponent_values = Displacementandstraincalc_results(elemIdx).Stress(:,:,2);
        
        % Filter out any NaN or inf values
        validValues = stressComponent_values(isfinite(stressComponent_values));
        if ~isempty(validValues)
            minStress = min(minStress, min(validValues));
            maxStress = max(maxStress, max(validValues));
        end
        
        % Create filled contour for this element using pcolor
        pcolor(x, y, stressComponent_values);
        shading interp;
    end
    
    % Add colorbar
    cb = colorbar;
    
    % Set colorbar range
    if isfinite(minStress) && isfinite(maxStress) && (maxStress > minStress)
        clim([minStress, maxStress]);
    else
        % Default range if we have invalid or constant values
        clim([-1, 1]);
    end
    
    % Add element outlines for clarity
    for elemIdx = 1:model.numElements
        elementNodes = model.Connectivity(elemIdx, :);
        nodeCoords = model.GlobalNodes(elementNodes, :);
        
        % Close the element by repeating the first node
        nodeCoords = [nodeCoords; nodeCoords(1, :)];
        
        % Plot element outline
        plot(nodeCoords(:, 1), nodeCoords(:, 2), 'k-', 'LineWidth', 0.5);
    end
    
    % Add text with max/min values
    if isfinite(minStress) && isfinite(maxStress)
        text(0.02*problemWidth, 0.02*problemHeight, sprintf('Min: %.2f', minStress), ...
             'FontSize', 10, 'BackgroundColor', 'white');
        text(0.02*problemWidth, 0.06*problemHeight, sprintf('Max: %.2f', maxStress), ...
             'FontSize', 10, 'BackgroundColor', 'white');
    end
    
    % Set axis limits
    axis equal;
    axis([0, problemWidth, 0, problemHeight]);
    
    % Save the stress plot
    saveas(gcf, sprintf('4node_sigma_yy_size_%f.png', meshSizes(i)));
    
    % Store timing
    cpuTimes4Node(i) = toc;
    
    % Store results for later analysis
    results4Node{i}.model = model;
    results4Node{i}.elements_structure = elements_structure;
    results4Node{i}.Global_elements_structure = Global_elements_structure;
    results4Node{i}.K = K;
    results4Node{i}.F = F;
    results4Node{i}.Q = Q;
    results4Node{i}.Displacementandstraincalc = Displacementandstraincalc_results;
    results4Node{i}.meshSize = meshSizes(i);
    results4Node{i}.minStress = minStress;
    results4Node{i}.maxStress = maxStress;
    
    % Track element count
    elementCounts4Node(i) = model.numElements;
    
    % Report timing
    fprintf('Mesh size = %.4f, Number of elements = %d, CPU time = %.4f seconds\n', ...
            meshSizes(i), model.numElements, cpuTimes4Node(i));
end

% Part 2: Analysis with 8-node elements
fprintf('\n===== ANALYSIS WITH 8-NODE QUADRATIC ELEMENTS =====\n');

for i = 1:length(meshSizes)
    fprintf('\nProcessing 8-node mesh with size = %f\n', meshSizes(i));
    
    % Start timing
    tic;
    
    % Create model with 8-node elements
    model = create8NodeModel(meshSizes(i), problemWidth, problemHeight, E, v, thickness);
    
    % Visualize mesh
    figure;
    plotMesh(model);
    title(sprintf('8-Node Mesh with element size = %f', meshSizes(i)));
    saveas(gcf, sprintf('8node_mesh_size_%f.png', meshSizes(i)));
    
    % Apply concentrated force at top-right corner
    model = applyConcentratedForce(model, forceValue);
    
    % Create elements structure
    elements_structure = createElementsStructure(model);
    
    % Create global stiffness matrix and force vector
    Global_elements_structure = Global_Stifness_Matrix(@convert_to_index, elements_structure, model);
    
    % Apply boundary conditions
    [K, F] = assemblyBC(model, Global_elements_structure);
    
    % Solve the system (K*Q = F)
    Q = K \ F;
    Global_elements_structure.Q = Q;
    
    % Calculate displacements and stresses
    Displacementandstraincalc_results = Displacementandstraincalc(Global_elements_structure, model);
    


    % Plot sigma_yy field for 8-node elements
figure;
hold on;
title(sprintf('\\sigma_{yy} field with 8-node elements, size = %f', meshSizes(i)));
xlabel('X-coordinate');
ylabel('Y-coordinate');
colormap(jet);

% Find min/max stress values
minStress = Inf;
maxStress = -Inf;

% Create a regular grid for interpolation
gridSize = 50;
[xq, yq] = meshgrid(linspace(0, problemWidth, gridSize), linspace(0, problemHeight, gridSize));
vq = nan(size(xq));

% Collect stress data from all elements
for elemIdx = 1:model.numElements
    % Get stress data and coordinates
    elemStress = Displacementandstraincalc_results(elemIdx).Stress(:,:,2);
    elemX = Displacementandstraincalc_results(elemIdx).x;
    elemY = Displacementandstraincalc_results(elemIdx).y;
    
    % Update min/max values
    validValues = elemStress(isfinite(elemStress));
    if ~isempty(validValues)
        minStress = min(minStress, min(validValues));
        maxStress = max(maxStress, max(validValues));
    end
    
    % For each point in the grid, check if it's in this element
    elemNodes = model.Connectivity(elemIdx, :);
    elemCoords = model.GlobalNodes(elemNodes(1:4), :); % Just corner nodes for boundary check
    
    % Find points inside this element
    in = inpolygon(xq, yq, elemCoords(:,1), elemCoords(:,2));
    
    % For points inside, estimate stress by finding closest point in element data
    [yi, xi] = find(in);
    for j = 1:length(yi)
        % Current grid point
        px = xq(yi(j), xi(j));
        py = yq(yi(j), xi(j));
        
        % Find closest point in element data
        [~, idx] = min(sqrt((elemX(:) - px).^2 + (elemY(:) - py).^2));
        [r, c] = ind2sub(size(elemX), idx);
        
        % Assign stress value
        if isfinite(elemStress(r, c))
            vq(yi(j), xi(j)) = elemStress(r, c);
        end
    end
end

% Create contour plot
contourf(xq, yq, vq, 20, 'LineStyle', 'none');

% Add colorbar with proper range
cb = colorbar;
if isfinite(minStress) && isfinite(maxStress) && minStress < maxStress
    caxis([minStress, maxStress]);
else
    % Default range if we have issues
    caxis([-1, 1]);
end

% Add element outlines for clarity
for elemIdx = 1:model.numElements
    elementNodes = model.Connectivity(elemIdx, :);
    
    % Plot corner nodes and edges with thinner lines
    for j = 1:4
        n1 = elementNodes(j);
        n2 = elementNodes(mod(j, 4) + 1);
        n5 = elementNodes(j + 4);
        
        % Plot corner to mid-side node to corner
        plot([model.GlobalNodes(n1,1), model.GlobalNodes(n5,1), model.GlobalNodes(n2,1)], ...
             [model.GlobalNodes(n1,2), model.GlobalNodes(n5,2), model.GlobalNodes(n2,2)], ...
             'k-', 'LineWidth', 0.5);
    end
end

% Add text with max/min values
if isfinite(minStress) && isfinite(maxStress)
    text(0.02*problemWidth, 0.02*problemHeight, sprintf('Min: %.2f', minStress), ...
         'FontSize', 10, 'BackgroundColor', 'white');
    text(0.02*problemWidth, 0.06*problemHeight, sprintf('Max: %.2f', maxStress), ...
         'FontSize', 10, 'BackgroundColor', 'white');
end

% Set axis properties
axis equal;
axis([0, problemWidth, 0, problemHeight]);

% Save figure
saveas(gcf, sprintf('8node_sigma_yy_size_%f.png', meshSizes(i)));



    
    % Store timing
    cpuTimes8Node(i) = toc;
    
    % Store results for later analysis
    results8Node{i}.model = model;
    results8Node{i}.elements_structure = elements_structure;
    results8Node{i}.Global_elements_structure = Global_elements_structure;
    results8Node{i}.K = K;
    results8Node{i}.F = F;
    results8Node{i}.Q = Q;
    results8Node{i}.Displacementandstraincalc = Displacementandstraincalc_results;
    results8Node{i}.meshSize = meshSizes(i);
    results8Node{i}.minStress = minStress;
    results8Node{i}.maxStress = maxStress;
    
    % Track element count
    elementCounts8Node(i) = model.numElements;
    
    % Report timing
    fprintf('Mesh size = %.4f, Number of elements = %d, CPU time = %.4f seconds\n', ...
            meshSizes(i), model.numElements, cpuTimes8Node(i));
end

% Part 3: Plot sigma_yy along x=2 for all mesh sizes and both element types
fprintf('\nPlotting stress profiles along x=2...\n');
figure;
hold on;
title('\sigma_{yy} along x=2 for different mesh sizes');
xlabel('y-coordinate');
ylabel('\sigma_{yy}');
grid on;

% Line styles and colors for different mesh sizes
lineStyles4Node = {'-o', '-s', '-^', '-d'};
lineStyles8Node = {'--o', '--s', '--^', '--d'};
lineColors = {'b', 'r', 'g', 'm'};
legendEntries = cell(length(meshSizes)*2, 1);

% Plot for 4-node elements
for i = 1:length(meshSizes)
    % Extract stress data along x=2
    [yCoords, stressValues] = extractStressProfile(results4Node{i}, 2);
    
    % Plot the profile
    plot(yCoords, stressValues, lineStyles4Node{i}, 'Color', lineColors{i}, 'LineWidth', 1.5);
    legendEntries{i} = sprintf('4-node, size = %.3f', meshSizes(i));
end

% Plot for 8-node elements
for i = 1:length(meshSizes)
    % Extract stress data along x=2
    [yCoords, stressValues] = extractStressProfile(results8Node{i}, 2);
    
    % Plot the profile
    plot(yCoords, stressValues, lineStyles8Node{i}, 'Color', lineColors{i}, 'LineWidth', 1.5);
    legendEntries{i+length(meshSizes)} = sprintf('8-node, size = %.3f', meshSizes(i));
end

% Add legend
legend(legendEntries, 'Location', 'best');
saveas(gcf, 'sigma_yy_vs_y_at_x2_comparison.png');

% Part 4: Plot maximum stress vs mesh size
fprintf('\nPlotting maximum stress vs mesh size...\n');
figure;
maxStress4Node = zeros(length(meshSizes), 1);
maxStress8Node = zeros(length(meshSizes), 1);

for i = 1:length(meshSizes)
    maxStress4Node(i) = results4Node{i}.maxStress;
    maxStress8Node(i) = results8Node{i}.maxStress;
end

loglog(meshSizes, maxStress4Node, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
loglog(meshSizes, maxStress8Node, 'rs-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Element Size');
ylabel('Maximum \sigma_{yy}');
title('Maximum \sigma_{yy} vs Element Size');
legend('4-node Linear Elements', '8-node Quadratic Elements', 'Location', 'best');
grid on;
saveas(gcf, 'max_stress_vs_mesh_size_comparison.png');

% Part 5: Plot timing information
fprintf('\nPlotting timing comparison...\n');
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

% Print timing comparison table
fprintf('\n=== Timing Comparison ===\n');
fprintf('Mesh Size | # Elements (4-node) | CPU Time (4-node) | # Elements (8-node) | CPU Time (8-node)\n');
fprintf('---------+--------------------+------------------+--------------------+------------------\n');
for i = 1:length(meshSizes)
    fprintf('  %.4f  |        %4d        |      %.4f s     |        %4d        |      %.4f s\n', ...
        meshSizes(i), elementCounts4Node(i), cpuTimes4Node(i), ...
        elementCounts8Node(i), cpuTimes8Node(i));
end

% Save all results
save('singularity_analysis_results.mat', 'results4Node', 'results8Node', 'meshSizes', ...
     'cpuTimes4Node', 'cpuTimes8Node', 'elementCounts4Node', 'elementCounts8Node');

fprintf('\nAnalysis complete. Results saved to singularity_analysis_results.mat\n');

% Helper function to extract stress profile along x=constant
function [yCoords, stressValues] = extractStressProfile(results, xValue)
    % Extract all nodes along x=xValue
    model = results.model;
    nodes = model.GlobalNodes;
    
    % Find nodes along x=xValue
    x2Nodes = find(abs(nodes(:,1) - xValue) < 1e-6);
    
    % Get y-coordinates and sort them
    yCoords = nodes(x2Nodes, 2);
    [yCoords, sortIdx] = sort(yCoords);
    x2Nodes = x2Nodes(sortIdx);
    
    % Initialize array for stress values
    stressValues = zeros(size(yCoords));
    
    % Get elements
    elements = model.Connectivity;
    numElements = model.numElements;
    
    % For each node, find and average stress values
    for j = 1:length(x2Nodes)
        node = x2Nodes(j);
        yVal = yCoords(j);
        
        % Find which elements contain this node
        containingElements = [];
        for elemIdx = 1:numElements
            if any(elements(elemIdx,:) == node)
                containingElements = [containingElements, elemIdx];
            end
        end
        
        % Average the stresses from elements that contain this node
        if ~isempty(containingElements)
            stressSum = 0;
            validCount = 0;
            
            for elemIdx = containingElements
                % Get the element's stress field for sigma_yy
                stressField = results.Displacementandstraincalc(elemIdx).Stress(:,:,2);
                
                % Get coordinates of the element
                elemX = results.Displacementandstraincalc(elemIdx).x;
                elemY = results.Displacementandstraincalc(elemIdx).y;
                
                % Find the closest point in the element's stress field to this node
                [~, minDistIdx] = min(sqrt((elemX(:) - xValue).^2 + (elemY(:) - yVal).^2));
                [row, col] = ind2sub(size(elemX), minDistIdx);
                
                % Add the stress value if it's valid
                stressValue = stressField(row, col);
                if isfinite(stressValue)
                    stressSum = stressSum + stressValue;
                    validCount = validCount + 1;
                end
            end
            
            % Calculate average stress (if any valid values found)
            if validCount > 0
                stressValues(j) = stressSum / validCount;
            else
                stressValues(j) = NaN;
            end
        else
            stressValues(j) = NaN;
        end
    end
    
    % Remove any NaN values
    validIdx = ~isnan(stressValues);
    yCoords = yCoords(validIdx);
    stressValues = stressValues(validIdx);
end