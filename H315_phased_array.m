% Niranjan Mathirajan, Owen Davis
% December 7, 2025
% Phased Array Spacing Playground
% Source 1: https://innovationspace.ansys.com/courses/wp-content/uploads/sites/5/2022/07/AntennaArrays_04ManipulatingtheArrayFactor.pdf
% Source 2: https://www.antenna-theory.com/arrays/arrayfactor.php
% Source 3: https://www.antenna-theory.com/definitions/wavevector.php


% Some notes about the simulation:
% the phase factor should always assume ideal spacing, as the computer
% doesn't know about the error in spacing. This is why the simulation and
% calculation's spacing don't incorporate the actual spacing, rather the
% ideal
%
% The expected calculation should never use the actual spacing, just ideal
% and error. The whole point is that the expected value incorporates the
% ideal distribution of dx_error (uniform), so that we can calculate the
% true error.
clc;
clf;

% -- constants --
T = 100; % number of trials performed
f = 5e6; % frequency of antenna (Hz)
c = 3e8; % speed of light (m/s)
lambda = c/f; % wavelength (m)
k = 2 * pi/lambda; % 1-D wave vector, see source (3) (rad/m)

% -- array parameters --
N = 16; % number of antennas
nTheta = 180;
theta = linspace(-pi/2, pi/2, nTheta); % range of angles the array "sees"
steerAngle = 30 * pi/180; % steering angle of array
dx_ideal = 0.25 * lambda; % expected distance between antennas
dx_error = 0.05 * lambda; % max variation in spacing

phase = -k * dx_ideal * cos(steerAngle); % term added to each AF to aim beam
measure_angle_deg = 30; % angle we measure AF at (range from 1 to 90)
measure_angle_rad = measure_angle_deg * 2*pi / 360;

% --- Monte Carlo Simulation ---

simulated_average = zeros(T,1);
for i = 1:T
    AF = zeros(nTheta,1);

    % --- get spacings ---
    dx_actual = zeros(N,1);
    for n = 1:N
        dx_actual(n) = dx_ideal - dx_error + 2*dx_error*rand();
    end

    % --- get Array Factor for every angle
    for t = 1:nTheta
        for n = 1:N
            % AF equation from source (1):
            main_term = k*dx_actual(n)*cos(theta(t));
            phase_term = - k*dx_ideal*cos(steerAngle);

            AF(t) = AF(t) + exp(1j*(n-1)*(main_term + phase_term));
        end
        AF(t) = AF(t)/N; % normalize AF, so max gain is 0dB
    end
    simulated_average(i) = abs(AF(measure_angle_deg+90));
end

simulated_average = mean(simulated_average)

% --- Expected calculation ---

E = 0;
for n = 1:N
    m = n-1;
    % these terms were pulled out of the array factor, should match
    % previous.
    main_term = exp(1j*k*m*cos(measure_angle_rad)*dx_ideal);
    phase_term = exp(-1j*k*m*dx_ideal*cos(steerAngle));

    if m == 0
        % when working at 0th antenna, the random variable doesn't work.
        % I'm guessing it's because there's no randomness in it's
        % placement? The practical reason is that the n-1 term causes a 0
        % in the numerator and denominator.
        random_variable = 1;
    else
        random_variable = (exp(1j*k*m*cos(measure_angle_rad)*dx_error) ...
            - exp(1j*k*m*cos(measure_angle_rad)*-1*dx_error))./(1j*k*m*cos(measure_angle_rad)*2*dx_error);
    end
    E = E + (main_term * phase_term * random_variable);
end

disp(abs(E/N))
% disp(AF(measure_angle_deg +90))


% -- plotting --
hps = polaraxes;
polarplot(theta, abs(AF),"blue", 'LineWidth',2); % plot AF
hold on;
% place dot on AF at our measure angle
polarplot(measure_angle_rad, abs(AF(measure_angle_deg + 90)), 'o', 'MarkerSize', 7, 'MarkerFaceColor', 'r', 'MarkerEdgeColor','black');
hps.ThetaZeroLocation = 'top'; % make 0 degrees show up at top of plot
ax = gca;
hps.ThetaDir = 'clockwise';
gca.FontSize = 12;
title('Normalized Array Factor (Unitless)');
grid on;



