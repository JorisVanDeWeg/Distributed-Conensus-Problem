function [Difference, sensorData_ADMM] = ADMM(adjMatrix, sensorPositions, sensorData, numIterations, Extras, dropoutInd, additionPos, rho)

    Difference = zeros(numIterations, 1);
    numSensors = size(sensorData, 1);
    averageVector = mean(sensorData) * ones(numSensors, 1);
    sensorData_ADMM = sensorData;
    Orig_adjMatrix = adjMatrix;

    % We assume it stays the same
    sensorRange = sqrt(2*log10(numSensors)/numSensors)*100; 

    z = zeros(numSensors, numSensors);
    nu = zeros(numSensors, numSensors);
    
    % Parameters
    % rho = 1.7;
    
    for k = 1:numIterations
        % Check for Transmission Error for Chosen Node
        % Extras(1,1) = Transmission error can occur
        % Extras(1,2) = Probability of transmission error
        if Extras(1,1) && rand(1) < Extras(1,2)
            % Compute the MSE of the difference between the vectors
            Difference(k) = mean((sensorData_ADMM - averageVector).^2);
            continue;
        end
        
        % Short Term Dropout per every Iteration 
        % Extras(2,1) = Short term dropout can occur (true/false)
        % Extras(2,2) = Probability of short term dropout
        % Extras(2,3) = Short term dropout occurs until iteration <.>
        if Extras(2,1) && k < Extras(2,3)
            adjMatrix = Orig_adjMatrix; % dropout takes only 1 iter to recover
            Dropout = (rand(size(Orig_adjMatrix,1),1) < Extras(2,2));
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
            sensorData_ADMM(dropoutInd,:) = [];
            % averageVector(dropoutInd) = [];
            averageVector = mean(sensorData) * ones(numSensors, 1);

            z(dropoutInd, :) = [];
            z(:, dropoutInd) = [];
            nu(dropoutInd, :) = [];
            nu(:, dropoutInd) = [];
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
            sensorData_ADMM = [sensorData_ADMM; new_sensorData];
            averageVector = mean(sensorData) * ones(numSensors, 1);
            
            z = [z; zeros(1, size(z,2))];
            z = [z, zeros(size(z, 1), 1)];
            nu = [nu; zeros(1, size(nu,2))];
            nu = [nu, zeros(size(nu, 1), 1)];
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
            sensorData_ADMM = [sensorData_ADMM; new_sensorData];
            averageVector = mean(sensorData) * ones(numSensors, 1);

            z = [z; zeros(Extras(5,3), size(z,2))];
            z = [z, zeros(size(z, 1), Extras(5,3))];
            nu = [nu; zeros(Extras(5,3), size(nu,2))];
            nu = [nu, zeros(size(nu, 1), Extras(5,3))];
        end

        % Randomly select a sensor i
        i = randi(numSensors);
    
        % Find its neighbors of node i
        neighbors = find(adjMatrix(i, :) == 1);

        commonElements = intersect(i, neighbors);
        if ~isempty(commonElements)
            Difference(k) = mean((sensorData_ADMM - averageVector).^2);
            continue;
        end

        if ~isempty(neighbors)
            % Randomly select one neighbor, j, uniformly
            j = neighbors(randi(length(neighbors)));
            
            % update x_i
            sensorData_ADMM(i) = (sensorData(i) + sum((rho*z(i,neighbors)) - nu(i,neighbors)))/(1+rho*length(neighbors));
            % sensorData_ADMM(i) = (sensorData(i) + rho*z(i,j) - nu(i,j))/(1+rho*length(neighbors));
            % Find its neighbors of node j
            neighbors_j = find(adjMatrix(j, :) == 1);
            
            % update x_j
            sensorData_ADMM(j) = (sensorData(j) + sum((rho*z(j,neighbors_j)) - nu(j,neighbors_j)))/(1+rho*length(neighbors_j));
            % sensorData_ADMM(j) = (sensorData(j) + rho*z(j,i) - nu(j,i))/(1+rho*length(neighbors_j));
            
            % update z_ij
            z(i,j) = 0.5*(sensorData_ADMM(i) + sensorData_ADMM(j));
            z(j,i) = 0.5*(sensorData_ADMM(i) + sensorData_ADMM(j));

            % update dual variables
            nu(i,j) = nu(i,j) + rho * (sensorData_ADMM(i) - z(i,j));
            nu(j,i) = nu(j,i) + rho * (sensorData_ADMM(j) - z(j,i));
        end
    
        Difference(k) = mean((sensorData_ADMM - averageVector).^2);
        
        if Difference(k) < 10^-12
            break;
        end
    end
    fprintf("%d \n",numSensors)
end