function elements_structure = createElementsStructure(model)
    % Create elements structure for FEA analysis
    % Works for both 4-node and 8-node elements
    
    connectivityMatrix = model.Connectivity;
    allNodalLocations = model.GlobalNodes;
    D = model.D;
    nDOFPNode = model.nDOFPNode;
    
    numElements = size(connectivityMatrix, 1);
    numNodes = size(connectivityMatrix, 2); % 4 or 8 nodes per element
    
    % Initialize output structure
    elements_structure(numElements).shapeFunctions = [];
    elements_structure(numElements).DiffShapefunctions = [];
    elements_structure(numElements).Stiffnessmatrix = [];
    elements_structure(numElements).SurfaceForce = [];
    elements_structure(numElements).x = [];
    elements_structure(numElements).y = [];
    
    % Setup Gaussian quadrature points and weights
    if numNodes == 4
        Gauss_in_Xi = 2;  % Number of Gauss points in ξ direction for 4-node elements
        Gauss_in_Eta = 2; % Number of Gauss points in η direction for 4-node elements
        [Integration_points, Weights] = Gauss_quadrant(2); % 2-point rule for 4-node elements
    else % numNodes == 8
        Gauss_in_Xi = 3;  % Number of Gauss points in ξ direction for 8-node elements
        Gauss_in_Eta = 3; % Number of Gauss points in η direction for 8-node elements
        [Integration_points, Weights] = Gauss_quadrant(3); % 3-point rule for 8-node elements
    end
    
    % Loop through each element
    for i = 1:numElements
        % Get element data
        nodeIndices = connectivityMatrix(i, :);
        X = allNodalLocations(nodeIndices, :);
        X = reshape(X.', [], 1);
        
        try
            % Calculate element stiffness matrix
            K_temp = Calculate_K_matrix(numNodes, X, D, Gauss_in_Xi, Gauss_in_Eta, Integration_points, Weights, model);
            
            % Initialize surface force vector
            SurfaceForce = zeros(nDOFPNode * numNodes, 1);
            
            % Calculate shape functions and coordinates at integration points
            n = 10; % Grid size for visualization
            x = zeros(n+1, n+1);
            y = zeros(n+1, n+1);
            
            % Store shape functions for last point (will be overwritten in loop)
            shapeFunctions = [];
            DiffShapefunctions = [];
            
            % Loop through visualization grid
            for j = 1:n+1
                Xi_j = -1 + 2*(j-1)/n;
                for k = 1:n+1
                    Eta_k = -1 + 2*(k-1)/n;
                    
                    % Calculate shape functions
                    shapeFunctions = Shape_function([Xi_j, Eta_k], numNodes);
                    DiffShapefunctions = Diff_Shape_functions([Xi_j, Eta_k], numNodes);
                    
                    % Calculate physical coordinates
                    temp_C = calculate_coordinates(@Shape_function, [Xi_j, Eta_k], X, numNodes);
                    x(j, k) = temp_C(1, 1);
                    y(j, k) = temp_C(2, 1);
                end
            end
            
            % Store element properties in output structure
            elements_structure(i).SurfaceForce = SurfaceForce;
            elements_structure(i).shapeFunctions = shapeFunctions;
            elements_structure(i).DiffShapefunctions = DiffShapefunctions;
            elements_structure(i).Stiffnessmatrix = K_temp;
            elements_structure(i).x = x;
            elements_structure(i).y = y;
            
        catch err
            % Error handler for problematic elements
            warning('Problem creating element %d: %s', i, err.message);
            
            % Create default values
            elements_structure(i).SurfaceForce = zeros(nDOFPNode * numNodes, 1);
            elements_structure(i).shapeFunctions = zeros(2, nDOFPNode * numNodes);
            elements_structure(i).DiffShapefunctions = zeros(4, nDOFPNode * numNodes);
            elements_structure(i).Stiffnessmatrix = eye(nDOFPNode * numNodes);
            elements_structure(i).x = ones(n+1, n+1);
            elements_structure(i).y = ones(n+1, n+1);
        end
    end
end