% Niranjan Mathirajan, Owen Davis
% October 14, 2025
% ECE 315H MATLAB Phased Array Antenna 


% https://innovationspace.ansys.com/courses/wp-content/uploads/sites/5/2022/07/AntennaArrays_04ManipulatingtheArrayFactor.pdf

f = 5e6; % frequency of antenna 
c = 3e8; % speed of light
lambda = c/f;
dx = 0.25 * lambda; % distance between antennas
N = 10; % number of antennas
nTheta = 180;
theta = linspace(-pi/2, pi/2, nTheta);
steerAngle = 0 * pi/180; % direction of anetanna figure out this phase stuff later
phase = -k * dx * cos(steerAngle);
AF = zeros(nTheta,1);

% Array Factor value
for t = 1:nTheta
    for n = 1:N
        AF(t) = AF(t) + exp(1j*(n-1)*(k*dx*cos(theta(t)) + phase));
    end
    AF(t) = AF(t)/N;
end


polarplot(theta, abs(AF),"blue",...
    'LineWidth',2)

ax = gca; 
gca.FontSize = 12;

grid on