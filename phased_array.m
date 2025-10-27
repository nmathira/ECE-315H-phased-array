% Niranjan Mathirajan, Owen Davis
% October 16, 2025
% ECE 315H MATLAB Phased Array Antenna 

% https://innovationspace.ansys.com/courses/wp-content/uploads/sites/5/2022/07/AntennaArrays_04ManipulatingtheArrayFactor.pdf

f = 5e6; % frequency of antenna 
c = 3e8; % speed of light
lambda = c/f;
dx = 0.25 * lambda; % distance between antennas
N = 5; % number of antennas
nTheta = 180*2;
theta = linspace(0, 2 * pi, nTheta);
% theta = 0; % direction of anetanna figure out this phase stuff later
phase = pi/180 * -25; % phase shift
k = 2 * pi/ lambda;
AF = zeros(1, nTheta);

% Array Factor value, with phase
for t = 1:nTheta
    for n = 1:N
        AF(t) = AF(t) + exp(1j*(n-1)*(k*dx*cos(theta(t))+phase));
    end
    AF(t) = AF(t)/N;
end

% plotting 
figure('Position', [100, 100, 1200, 600]); % Adjusts window size
polarplot(theta, abs(AF), 'LineWidth', 1);
