% This code retrieves orbital velocity and period around Itokawa
% as well as propogating the orbit using orbital element inputs

% Gravitational Constants
G = 6.67430e-11; % Gravitational constant (m^3/kg/s^2)
mu_Itokawa = 2.33e-9; % Itokawa's standard gravitational parameter (mu) [km^3/s^2]

% Define Orbital Elements of Spacecraft around Itokawa
a = 1.165;                % Semi-major axis in km
e = 0.00;                 % Eccentricity
i = deg2rad(29);          % Inclination in radians
Omega = deg2rad(00);      % RA of Ascending Node in radians
omega = deg2rad(00);      % Argument of Periapsis in radians
nu = deg2rad(0);          % True Anomaly in radians
mu = mu_Itokawa;          % Gravitational parameter
M = 3.5E10;               % Mass of Itokawa (kg)

% Velocity and Period Calculator
% Orbital Velocity
v_orbital = sqrt(G * M / (a*1000)); % (m/s)

% Orbital Period
T_orbital = 2 * pi * sqrt((a*1000)^3 / (G * M)); % (seconds)

% Display Results
fprintf('Orbital Velocity: %.2f m/s\n', v_orbital); % Orbital Velocity (m/s)
fprintf('Orbital Period: %.2f seconds', T_orbital); % Period in seconds
fprintf(' (%.2f minutes)', T_orbital / 60); % Period in minutes
fprintf(' (%.2f hours)', T_orbital / 60 / 60); % Period in hours
fprintf(' (%.2f days)\n', T_orbital / 60 / 60 / 24); % Period in days

% Orbit Propogation 
% Convert orbital elements to state vectors (position and velocity)
[r0, v0] = orbital_elements_to_state_vector(a, e, i, Omega, omega, nu, mu);

% Set the time span for propagation
T_orbit = 2 * pi * sqrt(a^3 / mu);  % Orbital period based on Kepler's 3rd law
tspan = [0 2*T_orbit];              % Propagate for two orbital periods

% Set ODE options for numerical integration
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-9);

% Propagate the orbit using ode45
[t, rv] = ode45(@(t, rv) orbital_dynamics(t, rv, mu), tspan, [r0; v0], options);

% Extract position and velocity vectors
r = rv(:, 1:3);
v = rv(:, 4:6);

% Plot the propagated orbit
figure;
plot3(r(:,1), r(:,2), r(:,3), 'b');
hold on;
plot3(r0(1), r0(2), r0(3), 'ro'); % Initial position
xlabel('X (km)'); ylabel('Y (km)'); zlabel('Z (km)');

% Plot central body (Itokawa)
[x_sphere, y_sphere, z_sphere] = sphere;
radius_Itokawa = 0.165; % Radius of Itokawa in km
surf(radius_Itokawa * x_sphere, radius_Itokawa * y_sphere, radius_Itokawa * z_sphere, 'FaceColor', [0.5, 0.5, 0.5], 'EdgeColor', 'none');

% Enhance plot appearance
grid on;
title('Spacecraft Orbit Around Itokawa');
legend('Orbit Path', 'Initial Position', 'Central Body (Itokawa)');
axis equal;

%% Function to convert orbital elements to state vector
function [r, v] = orbital_elements_to_state_vector(a, e, i, Omega, omega, nu, mu)
    % Calculate the distance (r) and velocity (v) in the perifocal frame
    p = a * (1 - e^2);  % Semi-latus rectum
    r_pf = (p / (1 + e * cos(nu))) * [cos(nu); sin(nu); 0];
    v_pf = sqrt(mu / p) * [-sin(nu); e + cos(nu); 0];
    
    % Rotation matrices
    R3_W = [cos(-Omega) sin(-Omega) 0;
           -sin(-Omega) cos(-Omega) 0;
            0           0          1];
    R1_i = [1    0           0;
            0  cos(-i)  sin(-i);
            0 -sin(-i)  cos(-i)];
    R3_w = [cos(-omega) sin(-omega) 0;
           -sin(-omega) cos(-omega) 0;
            0           0          1];
    
    % Complete rotation matrix from perifocal to geocentric equatorial frame
    Q_pX = (R3_W * R1_i * R3_w)';
    
    % Calculate position and velocity in geocentric equatorial frame
    r = Q_pX * r_pf;
    v = Q_pX * v_pf;
end

%% Function to define the orbital dynamics
function drv_dt = orbital_dynamics(~, rv, mu)
    % Extract position and velocity vectors
    r = rv(1:3);
    v = rv(4:6);
    
    % Calculate the magnitude of the position vector
    r_norm = norm(r);
    
    % Two-body acceleration
    a = -mu / r_norm^3 * r;
    
    % Return the derivative of the state vector
    drv_dt = [v; a];
end
