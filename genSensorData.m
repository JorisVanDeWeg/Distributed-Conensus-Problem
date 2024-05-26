function sensorData = genSensorData(type, numSensors, sensorPositions)
    % Generating Sensor Data
    % Assigning Initial Value - Gaussian
    if strcmp(type, "Gaussian")
        sensorData = zeros(numSensors, 1);
        centerX = areaSize / 2;
        centerY = areaSize / 2;
        stdDev = areaSize / 6; % Can be changed
        for i = 1:numSensors
            distance = norm([centerX, centerY] - sensorPositions(i,:)); % Distance from center
            sensorData(i) = exp(-(distance^2) / (2 * stdDev^2)); % Gaussian for this distance
        end
        sensorData = (sensorData / max(sensorData)) * 100; % Normalize
    
    % Assigning Initial Value - Random Uniform
    elseif strcmp(type, "Random Uniform")
        sensorData = rand(numSensors, 1) * 100; % Simulated sensor data
    
    % Assigning Initial Value - Uniform/Gradient
    elseif strcmp(type, "Gradient")
        sensorData = zeros(numSensors, 1);
        % Fill in the sensorData with a linear gradient based on a diagonal from the bottom left to the top right
        for i = 1:numSensors
            % Create a weighted sum of the x and y coordinates
            sensorData(i) = (sensorPositions(i,1) + sensorPositions(i,2)) / (2 * areaSize) * 100;
        end
    else
        error('Wrong name tpye is given. \n The following can be used: Gaussian, Random Uniform or Gradient') 
    end
end
