function [Difference, sensorData_PDMM] = PDMM(adjMatrix, sensorPositions, sensorData, numIterations, Extras, dropoutInd, additionPos, gamma_p)

    numSensors = size(sensorData, 1);
    averageVector = mean(sensorData) * ones(numSensors, 1);
    Difference = zeros(numIterations, 1);
    sensorData_PDMM = sensorData;
    Orig_adjMatrix = adjMatrix;
    
    % We assume it stays the same
    sensorRange = sqrt(2*log10(numSensors)/numSensors)*100; 

    % Define step sizes for primal and dual updates in the slides this is c
    % gamma_p = 1
    gamma_d = 1/gamma_p; 
    
    % Initialization of the primal and dual variables
    lambda_hat = zeros(numSensors, numSensors); % Dual variable for each edge
    w = zeros(numSensors, 1);

    % Perform asynchronous updates
    for k = 1:numIterations
        % Check for Transmission Error for Chosen Node
        % Extras(1,1) = Transmission error can occur
        % Extras(1,2) = Probability of transmission error
        if Extras(1,1) && rand(1) < Extras(1,2)
            % Compute the MSE of the difference between the vectors
            Difference(k) = mean((sensorData_PDMM - averageVector).^2);
            continue;
        end

        % Short Term Dropout per every Iteration 
        % Extras(2,1) = Short term dropout can occur (true/false)
        % Extras(2,2) = Probability of short term dropout
        % Extras(2,3) = Short term dropout occurs until iteration <.>
        if Extras(2,1) && k < Extras(2,3)          
            adjMatrix = Orig_adjMatrix; % dropout takes only 1 iter to recover
            Dropout = (rand(size(Orig_adjMatrix,1),1) < Extras(2,2));
            % lambda_hat(Dropout,:) = [];
            % lambda_hat(:,Dropout) = [];
            % w(dropoutInd) = [];

            % Find indices where Dropout is true
            idx = find(Dropout == 1);
            % Set rows and columns to zero
            adjMatrix(idx, :) = 0;  
            adjMatrix(:, idx) = 0;
        end
        if Extras(2,1) && k == Extras(2,3)
            adjMatrix = Orig_adjMatrix; % dropout takes only 1 iter to recover
        end

        % Long Term Dropout
        % Extras(3,1) = Long term dropout can occur (true/false)
        % Extras(3,2) = Dropout will happen at iteration <.>
        % Extras(3,3) = Number of nodes that will dropout
        if Extras(3,1) && k == Extras(3,2)
            % Set rows and columns to zero
            adjMatrix(dropoutInd, :) = [];  
            adjMatrix(:, dropoutInd) = [];
            Orig_adjMatrix = adjMatrix;

            numSensors = numSensors - length(dropoutInd);
            sensorData(dropoutInd,:) = [];
            sensorData_PDMM(dropoutInd,:) = [];
            
            averageVector = mean(sensorData) * ones(numSensors, 1);
            lambda_hat(dropoutInd,:) = [];
            lambda_hat(:,dropoutInd) = [];
            w(dropoutInd) = [];
        end

        % Consistent Addition per every Iteration 
        % Extras(4,1) = Consistent addition can occur (true/false)
        % Extras(4,2) = Probability of consistent addition
        % Extras(4,3) = Consistent addition occurs until iteration <.>
        if Extras(4,1) && rand(1) < Extras(4,2) && k <= Extras(4,3)

            % Randomly place sensors in the area
            sensorPositions = [sensorPositions;
                                100 * rand(1,2)];
            numSensors = numSensors + 1;
            % Check for connectivity using graph theory
            adjMatrix = zeros(numSensors, numSensors);
            for i = 1:numSensors
                for j = i+1:numSensors
                    if norm(sensorPositions(i,:) - sensorPositions(j,:)) <= sensorRange
                        adjMatrix(i,j) = 1;
                        adjMatrix(j,i) = 1;
                    end
                end
            end
            % Generating Sensor Data
            new_sensorData = genSensorData("Random Uniform", 1, sensorPositions(end,:));
            sensorData = [sensorData; new_sensorData];
            sensorData_PDMM = [sensorData_PDMM; new_sensorData];
            averageVector = mean(sensorData) * ones(numSensors, 1);
            
            lambda_hat = [lambda_hat; zeros(1, size(lambda_hat,2))];
            lambda_hat = [lambda_hat, zeros(size(lambda_hat, 1), 1)];
            w = [w; 0];
        end

        % One time Addition 
        % Extras(5,1) = One time addition can occur (true/false)
        % Extras(5,2) = One time addition will happen at iteration <.>
        % Extras(5,3) = Number of nodes that will be added
        if Extras(5,1) && k == Extras(5,2)
            % Randomly place sensors in the area
            sensorPositions = [sensorPositions; 
                                additionPos];
            numSensors = numSensors + Extras(5,3);
            % Check for connectivity using graph theory
            adjMatrix = zeros(numSensors, numSensors);
            for i = 1:numSensors
                for j = i+1:numSensors
                    if norm(sensorPositions(i,:) - sensorPositions(j,:)) <= sensorRange
                        adjMatrix(i,j) = 1;
                        adjMatrix(j,i) = 1;
                    end
                end
            end
           
            % Generating Sensor Data
            new_sensorData = genSensorData("Random Uniform", Extras(5,3), additionPos);
            sensorData = [sensorData; new_sensorData];
            sensorData_PDMM = [sensorData_PDMM; new_sensorData];
            averageVector = mean(sensorData) * ones(numSensors, 1);

            lambda_hat = [lambda_hat; zeros(Extras(5,3), size(lambda_hat,2))];
            lambda_hat = [lambda_hat, zeros(size(lambda_hat, 1), Extras(5,3))];
            w = [w; zeros(Extras(5,3),1)];
        end
        
        % Select the node i to be activated according to paper round robin loop
        %i = mod(k, numSensors) + 1;
        i = randi(size(adjMatrix,1));
    
        % Get the neighbors of node i Assuming Connected graph so every node
        % has atleast 1 neighbor
        neighbors_i = find(adjMatrix(i, :) == 1);    
        Ni = length(neighbors_i); % Amount of neighbours of node i
        
        commonElements = intersect(i, neighbors_i);
        if ~isempty(commonElements)
            Difference(k) = mean((sensorData_ADMM - averageVector).^2);
            continue;
        end

        if ~isempty(neighbors_i)
            sum_x_hat_j = sum(sensorData_PDMM(neighbors_i));
            sum_lambda_hat_ji = 0;
        
            for j = neighbors_i
                if j > i
                    sum_lambda_hat_ji = sum_lambda_hat_ji + lambda_hat(j, i);
                else 
                    sum_lambda_hat_ji = sum_lambda_hat_ji - lambda_hat(j, i);
                end
            end
        
            % Compute the auxiliary variable w for node i 
            w(i) = (sum_x_hat_j + gamma_d * sum_lambda_hat_ji + gamma_d * sensorData(i)) /( Ni + gamma_d);
            
            % Update x_hat for node i using equation (64)
            sensorData_PDMM(i) = (sensorData(i) + gamma_p * sum_x_hat_j + sum_lambda_hat_ji) / (1 + Ni * gamma_p);
            %fprintf('Difference %.2f\n',   w(i) - x_hat(i));
        
            % Update lambda_hat for edges connected to node i using equation (65)
            for j = neighbors_i
                if j > i
                    lambda_hat(i, j) = lambda_hat(j, i) - (1 / gamma_d) * (w(i) - sensorData_PDMM(j)); 
                else 
                    lambda_hat(i, j) = lambda_hat(j, i) - (1 / gamma_d) * (-w(i) + sensorData_PDMM(j));
                end
            end
            %lambda_hat(neighbors_i, i) = lambda_hat(i, neighbors_i);
        end
        % Track the L2 norm difference for convergence
        % Difference(k) = norm(sensorData_PDMM - averageVector,2)^2;

        Difference(k) = mean((sensorData_PDMM - averageVector).^2);
        if Difference(k) < 10^-12
            break;
        end
    end
    fprintf("%d \n",numSensors)
end