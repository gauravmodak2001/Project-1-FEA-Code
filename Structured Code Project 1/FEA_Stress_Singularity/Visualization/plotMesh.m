function plotMesh(model)
    % Function to plot the mesh (nodes and elements)
    % Works for both 4-node and 8-node elements
    
    % Get node coordinates
    nodes = model.GlobalNodes;
    
    % Get element connectivity
    connectivity = model.Connectivity;
    numNodesPerElem = size(connectivity, 2); % 4 or 8
    
    % Create figure
    hold on;
    grid on;
    
    % Define colors for elements (cycling through a set of colors)
    colorSet = ['r', 'g', 'b', 'c', 'm', 'y'];
    
    % Plot each element
    for i = 1:model.numElements
        % Get element nodes
        elementNodes = connectivity(i, :);
        nodeCoords = nodes(elementNodes, :);
        
        if numNodesPerElem == 4
            % For 4-node elements, draw straight lines between corners
            corners = [nodeCoords(1:4,:); nodeCoords(1,:)]; % Close the element
            
            % Plot element with specific color
            colorIdx = mod(i-1, length(colorSet)) + 1;
            plot(corners(:,1), corners(:,2), colorSet(colorIdx), 'LineWidth', 1);
            
        else % 8-node elements
            % For 8-node elements, we need to draw the shape with midside nodes
            % Plot corner nodes and edges
            for j = 1:4
                n1 = j;
                n2 = mod(j, 4) + 1;
                n5 = j + 4;
                
                % Plot corner to mid-side node
                plot([nodeCoords(n1, 1), nodeCoords(n5, 1)], ...
                     [nodeCoords(n1, 2), nodeCoords(n5, 2)], ...
                     colorSet(mod(i-1, length(colorSet)) + 1), 'LineWidth', 1);
                
                % Plot mid-side node to corner node
                plot([nodeCoords(n5, 1), nodeCoords(n2, 1)], ...
                     [nodeCoords(n5, 2), nodeCoords(n2, 2)], ...
                     colorSet(mod(i-1, length(colorSet)) + 1), 'LineWidth', 1);
            end
        end
        
        % Calculate element center for element number
        centerX = mean(nodeCoords(1:min(4,numNodesPerElem), 1));
        centerY = mean(nodeCoords(1:min(4,numNodesPerElem), 2));
        text(centerX, centerY, sprintf('%d', i), 'FontSize', 8, 'Color', 'k', ...
             'HorizontalAlignment', 'center');
    end
    
    % Plot and label all nodes
    scatter(nodes(:, 1), nodes(:, 2), 30, 'k', 'filled');
    
    % Label nodes - but only if not too many
    if size(nodes, 1) < 100
        for j = 1:size(nodes, 1)
            text(nodes(j, 1), nodes(j, 2), sprintf(' %d', j), 'FontSize', 8, 'Color', 'b');
        end
    end
    
    % Add labels
    xlabel('X-coordinate');
    ylabel('Y-coordinate');
    
    % Set axis limits with some padding
    maxWidth = max(nodes(:, 1));
    maxHeight = max(nodes(:, 2));
    axis([-0.1, maxWidth*1.1, -0.1, maxHeight*1.1]);
    
    % Make axis equal for proper visualization
    axis equal;
    hold off;
end