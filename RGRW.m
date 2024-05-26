function [Difference, sensorData_RGRW] = RGRW(adjMatrix, sensorPositions, sensorData, numIterations, Extras, dropoutInd, additionPos, p_restart)

    Difference = zeros(numIterations, 1);
    numSensors = size(sensorData, 1);
    averageVector = mean(sensorData) * ones(numSensors, 1);
    sensorData_RGRW = sensorData;
    Orig_adjMatrix = adjMatrix;

    % We assume it stays the same
    sensorRange = sqrt(2*log10(numSensors)/numSensors)*100; 
    % p_restart = 0.50; % Tune
    
    for k = 1:numIterations
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
            sensorData_RGRW(dropoutInd,:) = [];
            
            averageVector = mean(sensorData_RGRW) * ones(numSensors, 1);
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
            sensorData_RGRW = [sensorData_RGRW; new_sensorData];
            averageVector = mean(sensorData) * ones(numSensors, 1);
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
            sensorData_RGRW = [sensorData_RGRW; new_sensorData];
            averageVector = mean(sensorData) * ones(numSensors, 1);
        end
        
        % Randomly select a starting sensor
        activeNode = randi(numSensors);
        pathSensors = activeNode;
    
        visited = false(1, numSensors); % Track visited nodes
        visited(activeNode) = true;


        % Build the path of connected sensors
        while true
            % Check for Transmission Error for Chosen Node
            if Extras(1,1) && rand(1) < Extras(1,2)
                break;
            end

            neighbors = find(adjMatrix(pathSensors(end), :) == 1 & ~visited);
            % Exclude the sensors that are already in the path
            neighbors = setdiff(neighbors, pathSensors(end));
            
            % If no neighbors or time to restart; end the path 
            if isempty(neighbors) || rand() < p_restart
                break;
            end
    
            % Randomly select the next sensor from the remaining neighbors
            nextSensor = neighbors(randi(length(neighbors)));
            pathSensors = [pathSensors, nextSensor];
            visited(nextSensor) = true;
        end
        
        % Average the sensor data along the path
        if length(pathSensors) > 1 % Ensure there is more than one sensor in the path
            pathAverage = mean(sensorData_RGRW(pathSensors));
            
            % Update the sensor data for all sensors in the path
            sensorData_RGRW(pathSensors) = pathAverage;
        end
    
        % Compute the L2 norm of the difference between the vectors
        % Difference(k) = (norm(sensorData_RGRW - averageVector))^2;
        Difference(k) = mean((sensorData_RGRW - averageVector).^2);

        if Difference(k) < 10^-12
            break;
        end
    end
end