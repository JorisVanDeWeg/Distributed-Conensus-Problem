clear all;
close all;
clc;

%% Generating Sensor Graph
rng(1096) % Setting Random Seed


areaSize = 100;  % Side Length of the Square Area
numSensors = 200;  % Number of Sensors in the Area

% Communication range of each Sensor in Meters
sensorRange = sqrt(2*log10(numSensors)/numSensors)*areaSize; 

fprintf('Range of Sensors: %.2f \n', sensorRange);
fprintf('Number of Sensors: %d \n', numSensors);

% Randomly place sensors in the area
sensorPositions = areaSize * rand(numSensors, 2);

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

G = graph(adjMatrix);
[bin, binsize] = conncomp(G);
if length(binsize) == 1
    disp('The network is connected.');
else
    disp('The network is not connected.');
end


% Options for type of generation: 
opt = ["Gaussian", "Random Uniform", "Gradient"];

% Generating Sensor Data
sensorData = genSensorData(opt(2), numSensors, sensorPositions);

% Optimal Sensor Values
averageVector = mean(sensorData) * ones(numSensors, 1);
fprintf('The average sensor data before the algorithm is: %.2f\n', mean(sensorData));


% Chosing Type of Simulation
numIterations = 1000000;

% Long term Dropout & One time Addition (consitency)
dropoutInd = randperm(numSensors, 20);
additionPos = areaSize * rand(20, 2);

%% Plot Sensor Graph

% Plotting sensor positions, connections and value
figure;
for i = 1:numSensors
    for j = i+1:numSensors
        distance = norm(sensorPositions(i,:) - sensorPositions(j,:));
        if distance <= sensorRange
            line([sensorPositions(i,1), sensorPositions(j,1)], [sensorPositions(i,2), sensorPositions(j,2)], 'Color', "#4DBEEE");
        end
    end
end
hold on;
scatter(sensorPositions(:,1), sensorPositions(:,2), [], sensorData, 'filled');
colorbar;
axis([0 areaSize 0 areaSize]);
title('Sensor Network Placement and Parameters');
xlabel('X position (m)');
ylabel('Y position (m)');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  PLOTTING
%% RG
alg_name = ["Ideal" , "Transmission Fail", "Dropout", "Random Addition", "Addition"];

Diff_RG = [];

Extras = init_vars("RG");
[dr, ] =  RG(adjMatrix, sensorPositions, sensorData, numIterations, Extras, dropoutInd, additionPos);
Diff_RG = [Diff_RG, dr];
for ex = [1,3,4,5]
    zer = zeros(5,1);
    zer(ex) = 1;
    Extras(:,1) = zer;
    [dr, ] =  RG(adjMatrix, sensorPositions, sensorData, numIterations, Extras, dropoutInd, additionPos);
    Diff_RG = [Diff_RG, dr];
end

totalPlotter(Diff_RG, alg_name, "RG")

%% RGRW

p_restart = 0.5;
alg_name = ["Ideal" , "Transmission Fail", "Dropout", "Random Addition", "Addition"];
Diff_RGRW = [];

Extras = init_vars("RGRW");

[dr, ] =  RGRW(adjMatrix, sensorPositions, sensorData, numIterations, Extras, dropoutInd, additionPos, p_restart);
Diff_RGRW = [Diff_RGRW, dr];
for ex = [1,3,4,5]
    zer = zeros(5,1);
    zer(ex) = 1;
    Extras(:,1) = zer;
    [dr, ] =  RGRW(adjMatrix, sensorPositions, sensorData, numIterations, Extras, dropoutInd, additionPos, p_restart);
    Diff_RGRW = [Diff_RGRW, dr];
end

totalPlotter(Diff_RGRW, alg_name, "RGRW (p_r=" + p_restart + ")")

% Different values
Extras(:,1) = [0,0,0,0,0];
alg_name = [];
Diff_RGRW_pr = [];
for pr = [0.3,0.4, 0.5, 0.6, 0.7]
    alg_name = [alg_name, "pr=" + pr];
    [dr, ] = RGRW(adjMatrix, sensorPositions, sensorData, numIterations, Extras, dropoutInd, additionPos, pr);
    Diff_RGRW_pr = [Diff_RGRW_pr, dr];
end

totalPlotter(Diff_RGRW_pr, alg_name, "RGRW")

%% PDMM

gamma_p = 0.5;
alg_name = ["Ideal" , "Transmission Fail", "Dropout", "Random Addition", "Addition"];
Diff = [];

Extras = init_vars("PDMM");

[dr, ] =  PDMM(adjMatrix, sensorPositions, sensorData, numIterations, Extras, dropoutInd, additionPos, gamma_p);
Diff = [Diff, dr];
for ex = [1,3,4,5]
    zer = zeros(5,1);
    zer(ex) = 1;
    Extras(:,1) = zer;
    [dr, ] =  PDMM(adjMatrix, sensorPositions, sensorData, numIterations, Extras, dropoutInd, additionPos, gamma_p);
    Diff = [Diff, dr];
end
totalPlotter(Diff, alg_name, "PDMM (c=" + gamma_p + ")")
%%
% Different values
Extras = init_vars("PDMM");
Extras(:,1) = [0,0,0,0,0];
alg_name = [];
Diff = [];
for gamma_p = [0.3 0.4, 0.5, 0.6, 0.7, 0.8, 1.0]
    alg_name = [alg_name, "c=" + gamma_p ];
    [dr, ] = PDMM(adjMatrix, sensorPositions, sensorData, numIterations, Extras, dropoutInd, additionPos, gamma_p);
    Diff = [Diff, dr];
end

totalPlotter(Diff, alg_name, "PDMM")

%% ADMM 

rho = 0.5;
alg_name = ["Ideal" , "Transmission Fail", "Dropout", "Random Addition", "Addition"];
Diff = [];

Extras = init_vars("ADMM");

[dr, ] =  ADMM(adjMatrix, sensorPositions, sensorData, numIterations, Extras, dropoutInd, additionPos, rho);
Diff = [Diff, dr];
for ex = [1,3,4,5]
    zer = zeros(5,1);
    zer(ex) = 1;
    Extras(:,1) = zer;
    [dr, ] =  ADMM(adjMatrix, sensorPositions, sensorData, numIterations, Extras, dropoutInd, additionPos, rho);
    Diff = [Diff, dr];
end

totalPlotter(Diff, alg_name, "ADMM (\rho=" + rho + ")")
%%
% Different values
Extras(:,1) = [0,0,0,0,0];
alg_name = [];
Diff = [];
for rho = [0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2]
    alg_name = [alg_name, "\rho=" + rho];
    [dr, ] = ADMM(adjMatrix, sensorPositions, sensorData, numIterations, Extras, dropoutInd, additionPos, rho);
    Diff = [Diff, dr];
end

totalPlotter(Diff, alg_name, "ADMM")

%% Execution of Different Algorithms - END

Diff = [];
alg_name = [];
Extras(:,1) = [0,0,0,0,0];

% RG
alg_name =  [alg_name, "RG"];
[dr, ] = RG(adjMatrix, sensorPositions, sensorData, numIterations, Extras, dropoutInd, additionPos);
Diff = [Diff, dr];
% RGRW
for pr = [0.7, 0.6, 0.5]
    alg_name = [alg_name, "RGRW (p_r=" + pr + ")"];
    [dr, ] = RGRW(adjMatrix, sensorPositions, sensorData, numIterations, Extras, dropoutInd, additionPos, pr);
    Diff = [Diff, dr];
end
% PDMM
for gamma_p = [0.5, 1.0 , 1.5]
    alg_name = [alg_name, "PDMM (c=" + gamma_p + ")"];
    [dr, ] = PDMM(adjMatrix, sensorPositions, sensorData, numIterations, Extras, dropoutInd, additionPos, gamma_p);
    Diff = [Diff, dr];
end
% ADMM
for rho = [0.5, 1, 1.5]
    alg_name = [alg_name, "ADMM (\rho=" + rho + ")"];
    [dr, ] = ADMM(adjMatrix, sensorPositions, sensorData, numIterations, Extras, dropoutInd, additionPos, rho);
    Diff = [Diff, dr];
end

%%
% Plot Combined Results
totalPlotter(Diff, alg_name, "different algorithms")
%%
improvedPlotter(Diff, alg_name, "different algorithms")
