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
    
    % Plot sigma_yy field
    figure;
    plotStressField(Displacementandstraincalc_results, model, 2); % 2 for sigma_yy
    title(sprintf('\\sigma_{yy} field with 4-node elements, size = %f', meshSizes(i)));
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
    
    % Plot sigma_yy field
    figure;
    plotStressField(Displacementandstraincalc_results, model, 2); % 2 for sigma_yy
    title(sprintf('\\sigma_{yy} field with 8-node elements, size = %f', meshSizes(i)));
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
    
    % Track element count
    elementCounts8Node(i) = model.numElements;
    
    % Report timing
    fprintf('Mesh size = %.4f, Number of elements = %d, CPU time = %.4f seconds\n', ...
            meshSizes(i), model.numElements, cpuTimes8Node(i));
end

% Part 3: Plot sigma_yy along x=2 for all mesh sizes and both element types
plotSigmaYYvsY(results4Node, results8Node, meshSizes);

% Part 4: Plot maximum stress vs mesh size
plotMaxStressVsMeshSize(results4Node, results8Node, meshSizes);

% Part 5: Plot timing information
plotTimingComparison(elementCounts4Node, elementCounts8Node, cpuTimes4Node, cpuTimes8Node);

% Save all results
save('singularity_analysis_results.mat', 'results4Node', 'results8Node', 'meshSizes', ...
     'cpuTimes4Node', 'cpuTimes8Node', 'elementCounts4Node', 'elementCounts8Node');

fprintf('\nAnalysis complete. Results saved to singularity_analysis_results.mat\n');