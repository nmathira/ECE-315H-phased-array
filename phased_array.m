% Niranjan Mathirajan, Owen Davis
% October 14, 2025
% ECE 315H MATLAB Phased Array Antenna 


% https://innovationspace.ansys.com/courses/wp-content/uploads/sites/5/2022/07/AntennaArrays_04ManipulatingtheArrayFactor.pdf

f = 5e6; % frequency of antenna 
c = 3e8; % speed of light
lambda = c/f;
dx = 0.1 * lambda; % distance between antennas
N = 2; % number of antennas
nTheta = 180;
theta = linspace(0,pi,nTheta);
% theta = 0; % direction of anetanna figure out this phase stuff later
k = 2 * pi/ lambda;
AF = zeros(size(nTheta));

% Array Factor value
for t = 1:length(nTheta)
    for n = 1:N
        AF(t) = AF(t) + exp(1j*(n-1)*(k*dx*cos(theta(t))));
    end
end