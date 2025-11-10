% Niranjan Mathirajan, Owen Davis
% October 14, 2025
% ECE 315H MATLAB Phased Array Antenna 


% https://innovationspace.ansys.com/courses/wp-content/uploads/sites/5/2022/07/AntennaArrays_04ManipulatingtheArrayFactor.pdf

f = 5e6; % frequency of antenna 
c = 3e8; % speed of light
lambda = c/f;
k = 2 * pi/ lambda;
nTheta = 180;
N = 16; % number of antennas

% error_mean = 
error_variation = 0.05 * lambda;
dx = 0.25 * lambda;
dx_error = (rand(1,N) - 0.5) * error_variation;

dx_total = (dx + dx_error); % distance between antennas
theta = linspace(-pi/2, pi/2, nTheta);
steerAngle = 30 * pi/180; 
phase = -k * dx * cos(steerAngle);

% absolute antenna positions
% x = zeros(1, N); 
% for n = 1:N
%     x(n) = (n-1) * dx + dx_error(n);
% end


AF_ideal = zeros(nTheta,1);
AF_uneven = zeros(nTheta,1);

% Array Factor value
for t = 1:nTheta
    for n = 1:N
        AF_ideal(t) = AF_ideal(t) + exp(1j*(n-1)*(k*dx*cos(theta(t)) + phase));
        AF_uneven(t) = AF_uneven(t) + exp(1j*(n-1)*(k*dx_total(n)*cos(theta(t)) + phase));
    end
end

AF_ideal(t) = AF_ideal(t)/N;
AF_uneven(t) = AF_uneven(t)/N;
AF_error = AF_ideal - AF_uneven;

polarplot(theta, abs(AF_ideal),"blue",...
    theta, abs(AF_error),"red",...
    theta, abs(AF_uneven), "green",...
    'LineWidth',2)

% 
% polarplot(theta, abs(AF_ideal),"blue",...
%     'LineWidth',2)

ax = gca; 
gca.FontSize = 12;

grid on