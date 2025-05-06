function plotSigmaYYvsY(results4Node, results8Node, meshSizes)
    % Create a new figure for plotting sigma_yy vs y at x=2
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
    
    % Loop through each mesh size for 4-node elements
    for i = 1:length(meshSizes)
        % Get the mesh size for legend
        meshSize = meshSizes(i);
        legendEntries{i} = sprintf('4-node, size = %.3f', meshSize);
        
        % Extract data for this mesh
        Displacementandstraincalc_results = results4Node{i}.Displacementandstraincalc;
        nodeCoords = results4Node{i}.model.GlobalNodes;
        
        % Find nodes along x=2 (with small tolerance for numerical precision)
        xLine = 2.0;
        x2Nodes = find(abs(nodeCoords(:,1) - xLine) < 1e-6);
        
        % Get y-coordinates and sort them
        y2Coords = nodeCoords(x2Nodes, 2);
        [y2Coords, sortIdx] = sort(y2Coords);
        x2Nodes = x2Nodes(sortIdx);
        
        % Initialize arrays for storing stress values
        sigma_yy_values = zeros(size(y2Coords));
        
        % Get elements (elements that contain nodes on x=2)
        elements = results4Node{i}.model.Connectivity;
        numElements = results4Node{i}.model.numElements;
        
        % Loop through nodes at x=2 to extract stress values
        for j = 1:length(x2Nodes)
            node = x2Nodes(j);
            yVal = y2Coords(j);
            
            % Find which element contains this node
            containingElements = [];
            for elemIdx = 1:numElements
                if any(elements(elemIdx,:) == node)
                    containingElements = [containingElements, elemIdx];
                end
            end
            
            % Average the stresses from elements that contain this node
            if ~isempty(containingElements)
                stressSum = 0;
                for elemIdx = containingElements
                    % Get the element's stress field for sigma_yy
                    stressField = Displacementandstraincalc_results(elemIdx).Stress(:,:,2);
                    
                    % Get coordinates of the element
                    elemX = Displacementandstraincalc_results(elemIdx).x;
                    elemY = Displacementandstraincalc_results(elemIdx).y;
                    
                    % Find the closest point in the element's stress field to this node
                    [~, minDistIdx] = min(sqrt((elemX(:) - xLine).^2 + (elemY(:) - yVal).^2));
                    [row, col] = ind2sub(size(elemX), minDistIdx);
                    
                    % Add the stress value at this point
                    stressSum = stressSum + stressField(row, col);
                end
                
                % Calculate average stress
                sigma_yy_values(j) = stressSum / length(containingElements);
            end
        end
        
        % Plot sigma_yy vs y for this mesh size
        plot(y2Coords, sigma_yy_values, lineStyles4Node{i}, 'Color', lineColors{i}, 'LineWidth', 1.5, 'MarkerSize', 6);
    end
    
    % Loop through each mesh size for 8-node elements
    for i = 1:length(meshSizes)
        % Get the mesh size for legend
        meshSize = meshSizes(i);
        legendEntries{i+length(meshSizes)} = sprintf('8-node, size = %.3f', meshSize);
        
        % Extract data for this mesh
        Displacementandstraincalc_results = results8Node{i}.Displacementandstraincalc;
        nodeCoords = results8Node{i}.model.GlobalNodes;
        
        % Find nodes along x=2 (with small tolerance for numerical precision)
        xLine = 2.0;
        x2Nodes = find(abs(nodeCoords(:,1) - xLine) < 1e-6);
        
        % Get y-coordinates and sort them
        y2Coords = nodeCoords(x2Nodes, 2);
        [y2Coords, sortIdx] = sort(y2Coords);
        x2Nodes = x2Nodes(sortIdx);
        
        % Initialize arrays for storing stress values
        sigma_yy_values = zeros(size(y2Coords));
        
        % Get elements (elements that contain nodes on x=2)
        elements = results8Node{i}.model.Connectivity;
        numElements = results8Node{i}.model.numElements;
        
        % Loop through nodes at x=2 to extract stress values
        for j = 1:length(x2Nodes)
            node = x2Nodes(j);
            yVal = y2Coords(j);
            
            % Find which element contains this node
            containingElements = [];
            for elemIdx = 1:numElements
                if any(elements(elemIdx,:) == node)
                    containingElements = [containingElements, elemIdx];
                end
            end
            
            % Average the stresses from elements that contain this node
            if ~isempty(containingElements)
                stressSum = 0;
                for elemIdx = containingElements
                    % Get the element's stress field for sigma_yy
                    stressField = Displacementandstraincalc_results(elemIdx).Stress(:,:,2);
                    
                    % Get coordinates of the element
                    elemX = Displacementandstraincalc_results(elemIdx).x;
                    elemY = Displacementandstraincalc_results(elemIdx).y;
                    
                    % Find the closest point in the element's stress field to this node
                    [~, minDistIdx] = min(sqrt((elemX(:) - xLine).^2 + (elemY(:) - yVal).^2));
                    [row, col] = ind2sub(size(elemX), minDistIdx);
                    
                    % Add the stress value at this point
                    stressSum = stressSum + stressField(row, col);
                end
                
                % Calculate average stress
                sigma_yy_values(j) = stressSum / length(containingElements);
            end
        end
        
        % Plot sigma_yy vs y for this mesh size
        plot(y2Coords, sigma_yy_values, lineStyles8Node{i}, 'Color', lineColors{i}, 'LineWidth', 1.5, 'MarkerSize', 6);
    end
    
    % Add legend
    legend(legendEntries, 'Location', 'best');
    
    % Save the figure
    saveas(gcf, 'sigma_yy_vs_y_at_x2_comparison.png');
end