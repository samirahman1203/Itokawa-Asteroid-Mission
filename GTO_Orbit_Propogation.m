% Gravitational constant 
mu_earth = 3.99e5;     % Earth mu (mu) [km^3/s^2]

% Define orbital elements 
a = 24367.5;              % Semi-major axis in km
e = 0.728;                % Eccentricity
i = deg2rad(40);          % Inclination in radians
Omega = deg2rad(53.65);   % RA of Ascending Node in radians
omega = deg2rad(270);     % Argument of Periapsis in radians
nu = deg2rad(0);          % True Anomaly in radians
mu = mu_earth;            % Gravitational parameter

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

% Plot central body (Earth)
[x_sphere, y_sphere, z_sphere] = sphere;
radius_earth = 6378; % Radius of Earth in km
surf(radius_earth * x_sphere, radius_earth * y_sphere, radius_earth * z_sphere, 'FaceColor', [0.5, 0.5, 0.5], 'EdgeColor', 'none');

% Enhance plot appearance
grid on;
title('GTO Orbit Propagation');
legend('Orbit Path', 'Initial Position', 'Central Body (Earth)');
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
