% Niranjan Mathirajan, Owen Davis
% November 8, 2025
% ECE 315H MATLAB Phased Array Antenna 
% Source 1: https://innovationspace.ansys.com/courses/wp-content/uploads/sites/5/2022/07/AntennaArrays_04ManipulatingtheArrayFactor.pdf
% Source 2: https://www.antenna-theory.com/arrays/arrayfactor.php
% Source 3: https://www.antenna-theory.com/definitions/wavevector.php

clc;
clf;

% -- constants --
f = 5e6; % frequency of antenna (Hz)
c = 3e8; % speed of light (m/s)
lambda = c/f; % wavelength (m)
k = 2 * pi/lambda; % 1-D wave vector, see source (3) (rad/m)


% -- array parameters --
N = 8; % number of antennas
nTheta = 360; 
theta = linspace(-pi/2, pi/2, nTheta); % range of angles the array "sees"
steerAngle = 0 * pi/180; % steering angle of array
dx_ideal = 0.5 * lambda; % expected distance between antennas
dx_range = .1 * lambda; % max variation in spacing
dx_actual = (dx_ideal-dx_range) + (dx_range + dx_range).*rand(1,N); % compute random spacing
phase = -k * dx_ideal * cos(steerAngle); % term added to each AF to aim beam

% -- compute array factor -- 
AF = zeros(nTheta,1); 
for t = 1:nTheta
    for n = 1:N
        % AF equation from source (1):
        AF(t) = AF(t) + exp(1j*(n-1)*(k*dx_actual(n)*cos(theta(t)) + phase)); 
    end
    AF(t) = AF(t)/N; % normalize AF, so max gain is 0dB
end

% -- plotting --
hps = polaraxes;
polarplot(theta, abs(AF),"blue", 'LineWidth',2)
hps.ThetaZeroLocation = 'top'; % make 0 degrees show up at top of plot
ax = gca; 
gca.FontSize = 12;
title('Normalized Array Factor (dB)');
grid on;
