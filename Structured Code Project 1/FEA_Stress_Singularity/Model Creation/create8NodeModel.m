function model = create8NodeModel(meshSize, width, height, E, v, thickness)
    % Create a model with 8-node quadratic elements at specified mesh size
    
    % Initialize model data structure
    model = struct();
    model.E = E;
    model.v = v;
    model.thikness = thickness;  % Note: keep original spelling for compatibility
    model.nDOFPNode = 2;    % Degrees of freedom per node (x,y)
    
    % Define plane stress condition
    model.D = D_matrix_Stress(model);
    
    % Calculate number of elements in each direction
    numElemsX = ceil(width/meshSize);
    numElemsY = ceil(height/meshSize);
    
    % Adjust mesh size to fit domain exactly
    actualMeshSizeX = width/numElemsX;
    actualMeshSizeY = height/numElemsY;
    
    % For 8-node elements, we first create the 4-node base mesh
    % and then add mid-side nodes
    
    % Calculate number of corner nodes
    numNodesX_corner = numElemsX + 1;
    numNodesY_corner = numElemsY + 1;
    totalCornerNodes = numNodesX_corner * numNodesY_corner;
    
    % Generate corner node coordinates
    cornerNodeCoords = zeros(totalCornerNodes, 2);
    nodeCounter = 1;
    
    for j = 1:numNodesY_corner
        for k = 1:numNodesX_corner
            x = (k-1) * actualMeshSizeX;
            y = (j-1) * actualMeshSizeY;
            cornerNodeCoords(nodeCounter, :) = [x, y];
            nodeCounter = nodeCounter + 1;
        end
    end
    
    % Define 4-node element connectivity first
    elemConnectivity4 = zeros(numElemsX * numElemsY, 4);
    elemCounter = 1;
    
    for j = 1:numElemsY
        for k = 1:numElemsX
            n1 = (j-1)*numNodesX_corner + k;           % Bottom left
            n2 = (j-1)*numNodesX_corner + (k+1);       % Bottom right
            n3 = j*numNodesX_corner + (k+1);           % Top right
            n4 = j*numNodesX_corner + k;               % Top left
            
            elemConnectivity4(elemCounter, :) = [n1, n2, n3, n4];
            elemCounter = elemCounter + 1;
        end
    end
    
    % Now add mid-side nodes for 8-node quadratic elements
    
    % First, calculate total number of nodes in 8-node mesh
    totalMidSideNodes = numElemsX * (numNodesY_corner) + numElemsY * (numNodesX_corner);
    totalNodes = totalCornerNodes + totalMidSideNodes;
    
    % Initialize complete node array
    nodeCoords = zeros(totalNodes, 2);
    
    % Copy corner nodes to complete node array
    nodeCoords(1:totalCornerNodes, :) = cornerNodeCoords;
    
    % Define 8-node element connectivity
    elemConnectivity = zeros(numElemsX * numElemsY, 8);
    
    % Fill in mid-side nodes
    currentNodeIdx = totalCornerNodes + 1;
    
    % Create horizontal edges mid-side nodes first
    horizEdgeNodes = zeros(numElemsX, numNodesY_corner);
    
    for j = 1:numNodesY_corner
        for k = 1:numElemsX
            % Create mid-side node between corner nodes
            n1 = (j-1)*numNodesX_corner + k;           % Left node
            n2 = (j-1)*numNodesX_corner + (k+1);       % Right node
            
            % Calculate mid-side node coordinates
            x = (cornerNodeCoords(n1, 1) + cornerNodeCoords(n2, 1)) / 2;
            y = (cornerNodeCoords(n1, 2) + cornerNodeCoords(n2, 2)) / 2;
            
            nodeCoords(currentNodeIdx, :) = [x, y];
            horizEdgeNodes(k, j) = currentNodeIdx;
            currentNodeIdx = currentNodeIdx + 1;
        end
    end
    
    % Create vertical edges mid-side nodes
    vertEdgeNodes = zeros(numNodesX_corner, numElemsY);
    
    for j = 1:numElemsY
        for k = 1:numNodesX_corner
            % Create mid-side node between corner nodes
            n1 = (j-1)*numNodesX_corner + k;           % Bottom node
            n2 = j*numNodesX_corner + k;               % Top node
            
            % Calculate mid-side node coordinates
            x = (cornerNodeCoords(n1, 1) + cornerNodeCoords(n2, 1)) / 2;
            y = (cornerNodeCoords(n1, 2) + cornerNodeCoords(n2, 2)) / 2;
            
            nodeCoords(currentNodeIdx, :) = [x, y];
            vertEdgeNodes(k, j) = currentNodeIdx;
            currentNodeIdx = currentNodeIdx + 1;
        end
    end
    
    % Now build 8-node element connectivity
    for j = 1:numElemsY
        for k = 1:numElemsX
            elemIdx = (j-1)*numElemsX + k;
            
            % Get corner nodes from 4-node connectivity
            n1 = elemConnectivity4(elemIdx, 1);  % Bottom left
            n2 = elemConnectivity4(elemIdx, 2);  % Bottom right
            n3 = elemConnectivity4(elemIdx, 3);  % Top right
            n4 = elemConnectivity4(elemIdx, 4);  % Top left
            
            % Get mid-side nodes
            n5 = horizEdgeNodes(k, j);          % Bottom edge
            n6 = vertEdgeNodes(k+1, j);         % Right edge
            n7 = horizEdgeNodes(k, j+1);        % Top edge
            n8 = vertEdgeNodes(k, j);           % Left edge
            
            % Serendipity 8-node ordering: counterclockwise from bottom-left
            elemConnectivity(elemIdx, :) = [n1, n2, n3, n4, n5, n6, n7, n8];
        end
    end
    
    % Store mesh in data structure
    model.GlobalNodes = nodeCoords;
    model.Connectivity = elemConnectivity;
    model.numElements = size(elemConnectivity, 1);
    model.nTotalDOF = size(nodeCoords, 1) * model.nDOFPNode;
    model.type_of_mesh = "Undistorted";
    model.updateNodalLocations = nodeCoords;
    
    % Define boundary conditions
    % Find nodes on the bottom edge (y = 0)
    bottomNodes = find(abs(nodeCoords(:, 2)) < 1e-6);
    
    % Find nodes on the right edge (x = width)
    rightNodes = find(abs(nodeCoords(:, 1) - width) < 1e-6);
    
    % Initialize BC arrays
    totalBCs = length(bottomNodes) + length(rightNodes);
    model.BCnodes = zeros(totalBCs, 1);
    model.BCDof = zeros(totalBCs, 1);
    model.BC = zeros(totalBCs, 1);
    
    % Set boundary conditions:
    % - Bottom edge: restrict y-direction (vertical)
    % - Right edge: restrict x-direction (horizontal)
    bcCounter = 1;
    
    for i = 1:length(bottomNodes)
        model.BCnodes(bcCounter) = bottomNodes(i);
        model.BCDof(bcCounter) = 2;  % y-direction constraint
        model.BC(bcCounter) = 0;     % zero displacement
        bcCounter = bcCounter + 1;
    end
    
    for i = 1:length(rightNodes)
        model.BCnodes(bcCounter) = rightNodes(i);
        model.BCDof(bcCounter) = 1;  % x-direction constraint
        model.BC(bcCounter) = 0;     % zero displacement
        bcCounter = bcCounter + 1;
    end
    
    % Store mesh size in model for reference
    model.meshSize = meshSize;
    
    % Set no body forces initially
    model.f = [0, 0, 0, 0];
end