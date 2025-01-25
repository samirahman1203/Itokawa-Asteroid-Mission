% Constants
mu_Itokawa = 2.33e-9; % Gravitational parameter of Itokawa (Km^3/s^2)
au = 1.496e11; % Astronomical Unit in meters

% Orbit parameters
orbit_start_date = datetime(2037, 2, 8);
orbit_end_date = datetime(2039, 8, 6);
time_step = 0.0006944; % a minute denoted as days

% Communication parameters
eclipse_angle = 10; % degrees cone angle for eclipse analysis

% Initialise time array
time = orbit_start_date:days(time_step):orbit_end_date;
num_steps = length(time);

% Obtain ephemeris data for Itokawa
r_itokawa_sun = zeros(3, num_steps);
for i = 1:num_steps
    current_time = time(i);
    % Get Itokawa positions from ephemeris (in heliocentric frame)
    ephemeris_data_itokawa = Itokawa_Ephemeris('Itokawa');
    r_itokawa_sun(:, i) = ephemeris_data_itokawa(1:3); % Use first three elements as position vector
end
    r_itokawa_sun(:, i) = ephemeris_data_itokawa(1:3); % Use first three elements as position vector

% Assuming initial conditions for spacecraft orbiting Itokawa
% Set initial position and velocity around Itokawa
r0_itokawa = [0; 0; 1]; % Initial position 1 km away from Itokawa
v0_itokawa = [0; 0.04 / 1000; 0]; % Initial velocity for a stable orbit in km/s
y0 = double([r0_itokawa; v0_itokawa]); % Initial state vector (position and velocity)

% Integrate the spacecraft trajectory using initial conditions
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
[t, state] = ode45(@(t, y) two_body_equation(t, y, mu_Itokawa), double([0, days(orbit_end_date - orbit_start_date) * 86400]), y0, options); % Adjust time span for orbital period

% Extract spacecraft positions over time
r_spacecraft = state(:, 1:3)';

% Eclipse analysis around Itokawa asteroid
R_itokawa = 0.165; % Radius of Itokawa in km
eclipse_status = ones(1, length(t)); % Initially assume no eclipse

for i = 1:length(t)
    % Find the closest time index in the ephemeris data for current spacecraft time
    ephem_index = min(i, num_steps);
    % Vector from Itokawa to spacecraft
    r_itokawa = r_itokawa_sun(:, ephem_index); % Get Itokawa position at the closest available ephemeris time step
    r_spacecraft_current = r_spacecraft(:, i);
    r_itokawa_spacecraft = r_spacecraft_current - r_itokawa;

    % Vector from Itokawa to Sun (calculated dynamically based on Sun's position)
    r_itokawa_sun_vec = -r_itokawa;

    % Calculate whether the spacecraft falls in Itokawa's shadow
    % Step 1: Check if the spacecraft is behind Itokawa relative to the Sun
    if dot(r_itokawa_sun_vec, r_itokawa_spacecraft) > 0
        % Step 2: Check if the spacecraft is within the shadow cone
        % Calculate the angular radius of Itokawa as seen from the spacecraft
        D = norm(r_itokawa_spacecraft); % Distance between Itokawa and the spacecraft
        angular_radius_itokawa = atan(R_itokawa / D) * (180 / pi); % Angular radius in degrees

        % Calculate the angular separation between Itokawa-Sun vector and Itokawa-Spacecraft vector
        angular_separation = acosd(dot(r_itokawa_spacecraft, r_itokawa_sun_vec) / (norm(r_itokawa_spacecraft) * norm(r_itokawa_sun_vec)));

        % If the angular separation is less than the angular radius, the spacecraft is in shadow
        if angular_separation < angular_radius_itokawa
            eclipse_status(i) = 0; % Eclipse occurs
        end
    end
end

% Plotting eclipse status around Itokawa
time_days = t / 86400; % Convert seconds to days
figure;
plot(time_days, eclipse_status, 'LineWidth', 2);
ylabel('Eclipse Status');
yticklabels({'Eclipse', 'No Eclipse'});
yticks([0 1]);
xlabel('Time (days)');
title('Eclipse Window during Spacecraft Orbit around Itokawa');

grid on;

% Calculate eclipse hours for each day
days_array = unique(floor(time_days)); % Extract unique day numbers
num_days = length(days_array);
eclipse_hours = zeros(1, num_days);

for i = 1:num_days
    % Find indices corresponding to the current day
    day_indices = find(floor(time_days) == days_array(i));
    % Calculate the total time in eclipse for the current day
    eclipse_hours(i) = sum(eclipse_status(day_indices) == 0) * (time_step * 24); % Correct time conversion to hours
end

% Plotting eclipse hours per day
figure;
plot(days_array(1:end-1), eclipse_hours(1:end-1), 'LineWidth', 2);
xlabel('Day Number');
ylabel('Eclipse Duration (hours)');
title('Eclipse Duration per Day for Spacecraft Orbit around Itokawa');
grid on;
grid on;

% Function to define the two-body orbital motion equation
function dydt = two_body_equation(~, y, mu)
    if length(y) ~= 6
        error('State vector must have exactly 6 elements (position and velocity).');
    end
    r = y(1:3);
    v = y(4:6);
    r_norm = norm(r);

    % Derivative of the state vector
    dydt = zeros(6,1);
    dydt(1:3) = v;
    dydt(4:6) = -mu * r / r_norm^3;
end
