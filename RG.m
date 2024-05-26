function [Difference, sensorData_RG] = RG(adjMatrix, sensorPositions, sensorData, numIterations, Extras, dropoutInd, additionPos)
    
    numSensors = size(sensorData, 1);
    averageVector = mean(sensorData) * ones(numSensors, 1);
    Difference = zeros(numIterations, 1);
    sensorData_RG = sensorData;
    
    % We assume it stays the same
    sensorRange = sqrt(2*log10(numSensors)/numSensors)*100; 

    Orig_adjMatrix = adjMatrix;
    for k = 1:numIterations
        % Check for Transmission Error for Chosen Node
        % Extras(1,1) = Transmission error can occur
        % Extras(1,2) = Probability of transmission error
        if Extras(1,1) && rand(1) < Extras(1,2)
            % Compute the MSE of the difference between the vectors
           Difference(k) = mean((sensorData_RG - averageVector).^2);
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
            
            sensorPositions(dropoutInd,:) = [];
            numSensors = numSensors - Extras(3,3);

            sensorData(dropoutInd) = [];
            sensorData_RG(dropoutInd) = [];
            averageVector = mean(sensorData_RG) * ones(numSensors, 1);
            %averageVector(dropoutInd) = [];
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
            sensorData_RG = [sensorData_RG; new_sensorData];
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
            sensorData_RG = [sensorData_RG; new_sensorData];
            averageVector = mean(sensorData) * ones(numSensors, 1);
        end
        %sensor = mod(k, numSensors) + 1;
        sensor = randi(numSensors); % Randomly select a sensor
        neighbors = find(adjMatrix(sensor, :) == 1); % Find its neighbors

        if ~isempty(neighbors)
            % Randomly select one neighbor
            neighbor = neighbors(randi(length(neighbors)));
            
            % Average their sensor data
            averageValue = mean([sensorData_RG(sensor), sensorData_RG(neighbor)]);
            
            % Update both sensors' data
            sensorData_RG(sensor) = averageValue;
            sensorData_RG(neighbor) = averageValue;
        end

        % Compute the MSE of the difference between the vectors
        Difference(k) = mean((sensorData_RG - averageVector).^2);
        if Difference(k) < 10^-12
            break;
        end
    end
end
