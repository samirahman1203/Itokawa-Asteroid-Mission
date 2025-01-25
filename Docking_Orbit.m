% Docking Orbit for the Spacecraft onto Itokawa

% Set up the figure
figure;
axis equal;
hold on;
set(gca, 'Color', [1 1 1]); % Set background to light gray
grid on;

% Define Itokawa's radius 
radius_itokawa = 0.165; % in km 
[X, Y, Z] = sphere(50); % Create data for a sphere to represent Itokawa
surf(X * radius_itokawa, Y * radius_itokawa, Z * radius_itokawa, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none', 'DisplayName', 'Itokawa');

% Light and view setup
light('Position', [1 0 1], 'Style', 'infinite');
material dull;
view(3);

% Define gravitational parameters
G = 6.67430e-20; % Gravitational constant in km^3/kg/s^2
mass_itokawa = 3.5e10; % Mass of Itokawa in kg
mu = G * mass_itokawa; % Standard gravitational parameter for Itokawa in km^3/s^2

% Define spacecraft parameters
initial_orbit_radius = radius_itokawa + 1; % 1 km above Itokawa
initial_velocity = sqrt(mu / initial_orbit_radius); % Circular orbit velocity in km/s

% Define inclination angle
inclination_deg = 29; % Set inclination in degrees
inclination_rad = deg2rad(inclination_deg); % Convert to radians

% Define initial position and velocity with inclination
initial_position = [initial_orbit_radius * cos(inclination_rad); 0; initial_orbit_radius * sin(inclination_rad)]; % Inclined initial position
initial_velocity = [0; initial_velocity * cos(inclination_rad); initial_velocity * sin(inclination_rad)]; % Inclined initial velocity

% Combine into initial conditions for ODE solver
initial_conditions = [initial_position; initial_velocity];

% Define small thrust for docking
thrust_acceleration = -0.00000004; % Small constant thrust to bring the spacecraft closer to Itokawa (in km/s^2)

% Time span for ODE solver
tspan = [0 100000]; % Time span in seconds (increased for docking)

% Define the ODE function with thrust
gravitational_ode = @(t, y) [y(4); y(5); y(6); -mu * y(1) / norm(y(1:3))^3 + thrust_acceleration * y(1) / norm(y(1:3)); -mu * y(2) / norm(y(1:3))^3 + thrust_acceleration * y(2) / norm(y(1:3)); -mu * y(3) / norm(y(1:3))^3 + thrust_acceleration * y(3) / norm(y(1:3))];

% Solve the ODE using ode45
[t, solution] = ode45(gravitational_ode, tspan, initial_conditions);

% Extract positions from the solution
positions = solution(:, 1:3)';

% Find the index where the spacecraft intersects with Itokawa's surface
intersect_index = find(vecnorm(positions) <= radius_itokawa, 1);
if ~isempty(intersect_index)
    docking_position = positions(:, intersect_index);
    positions = positions(:, 1:intersect_index); % Cut off the docking orbit after intersection
else
    docking_position = radius_itokawa * [positions(1, end) / norm(positions(:, end)); positions(2, end) / norm(positions(:, end)); positions(3, end) / norm(positions(:, end))];
end

% Calculate and display descent duration
if ~isempty(intersect_index)
    docking_time = t(intersect_index); % Time at which docking occurs (in seconds)
    disp(['Descent Duration: ', num2str(docking_time), ' seconds (', num2str(docking_time / 60), ' minutes) (', num2str(docking_time / 60 / 60), ' hours)']);
else
    disp('The spacecraft did not intersect with the surface.');
end


% Plot the initial orbit trajectory with inclination
theta = linspace(0, 2*pi, 500);
initial_orbit_x = initial_orbit_radius * cos(theta) * cos(inclination_rad);
initial_orbit_y = initial_orbit_radius * sin(theta);
initial_orbit_z = initial_orbit_radius * cos(theta) * sin(inclination_rad);
plot3(initial_orbit_x, initial_orbit_y, initial_orbit_z, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Initial Orbital Path');

% Plot the spacecraft trajectory
plot3(positions(1, :), positions(2, :), positions(3, :), 'b', 'LineWidth', 1.5, 'DisplayName', 'Docking Orbit Path');

% Final docking representation at the intersection with Itokawa's surface
plot3(docking_position(1), docking_position(2), docking_position(3), 'go', 'MarkerSize', 5, 'MarkerFaceColor', 'g', 'DisplayName', 'Docking Point');

% Title and labels
title('Spacecraft Docking Onto Itokawa');
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
legend('show');

hold off;
