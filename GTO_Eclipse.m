% Constants
MuE = 3.986e5; % Gravitational parameter of Earth (Km^3/s^2)
RE = 6378; % Radius of Earth in km

% Orbit parameters
orbit_start_date = datetime(2036, 3, 16);
orbit_end_date = orbit_start_date + days(8); % One week period

% Time parameters
time_step = 0.0006944; % a minute in days
time = orbit_start_date:days(time_step):orbit_end_date;
num_steps = length(time);

% Communication parameters
eclipse_angle = 10; % degrees cone angle for eclipse analysis

% Orbital elements
semi_major_axis = 24367.5; % Semi-major axis in km
eccentricity = 0.728;
inclination_deg = 40; % Inclination in degrees
raan_deg = 53.65; % Right Ascension of Ascending Node (RAAN) in degrees
arg_perigee_deg = 0; % Argument of Perigee in degrees
true_anomaly_deg = 0; % True Anomaly at epoch in degrees

% Convert angles to radians
inclination_rad = deg2rad(inclination_deg);
raan_rad = deg2rad(raan_deg);
arg_perigee_rad = deg2rad(arg_perigee_deg);
true_anomaly_rad = deg2rad(true_anomaly_deg);

% Calculate initial position and velocity in perifocal coordinate system
p = semi_major_axis * (1 - eccentricity^2); % Semi-latus rectum in km
r_perifocal = p / (1 + eccentricity * cos(true_anomaly_rad));

% Position in perifocal frame
r0_perifocal = r_perifocal * [cos(true_anomaly_rad); sin(true_anomaly_rad); 0];

% Velocity in perifocal frame
v0_perifocal = sqrt(MuE / p) * [-sin(true_anomaly_rad); eccentricity + cos(true_anomaly_rad); 0];

% Rotation matrix from perifocal to ECI frame
R_peri_to_eci = [
    cos(raan_rad) * cos(arg_perigee_rad) - sin(raan_rad) * sin(arg_perigee_rad) * cos(inclination_rad), -cos(raan_rad) * sin(arg_perigee_rad) - sin(raan_rad) * cos(arg_perigee_rad) * cos(inclination_rad),  sin(raan_rad) * sin(inclination_rad);
    sin(raan_rad) * cos(arg_perigee_rad) + cos(raan_rad) * sin(arg_perigee_rad) * cos(inclination_rad), -sin(raan_rad) * sin(arg_perigee_rad) + cos(raan_rad) * cos(arg_perigee_rad) * cos(inclination_rad), -cos(raan_rad) * sin(inclination_rad);
    sin(arg_perigee_rad) * sin(inclination_rad),                                             cos(arg_perigee_rad) * sin(inclination_rad),                              cos(inclination_rad)
    ];

% Transform initial position and velocity to ECI frame
r0_earth = R_peri_to_eci * r0_perifocal;
v0_earth = R_peri_to_eci * v0_perifocal;

y0 = double([r0_earth; v0_earth]); % Initial state vector (position and velocity)

% Integrate the spacecraft trajectory using initial conditions
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
[t, state] = ode45(@(t, y) two_body_equation(t, y, MuE), double([0, days(orbit_end_date - orbit_start_date) * 86400]), y0, options); % Adjust time span for orbital period

% Extract spacecraft positions over time
r_spacecraft = state(:, 1:3)';

% Eclipse analysis around Earth
R_earth = RE; % Radius of Earth in km
eclipse_status = ones(1, length(t)); % Initially assume no eclipse

for i = 1:length(t)
    % Vector from Earth to spacecraft
    r_spacecraft_current = r_spacecraft(:, i);
    D = norm(r_spacecraft_current); % Distance between Earth and the spacecraft

    % Update Sun position over time
    orbital_period_earth = 365.25; % Earth orbital period in days
    angular_velocity_earth = 2 * pi / orbital_period_earth; % Angular velocity of Earth's orbit around Sun (rad/day)
    time_since_start = days(time(i) - orbit_start_date); % Time since start date in days
    r_earth_sun_vec = 1.496e8 * [cos(angular_velocity_earth * time_since_start); sin(angular_velocity_earth * time_since_start); 0]; % Sun's position changes over time

    % Step 1: Check if the spacecraft is behind Earth relative to the Sun
    if dot(r_earth_sun_vec, r_spacecraft_current) > 0
        % Step 2: Check if the spacecraft is within the shadow cone
        % Calculate the angular radius of Earth as seen from the spacecraft
        angular_radius_earth = atan(R_earth / D) * (180 / pi); % Angular radius in degrees

        % Calculate the angular separation between Earth-Sun vector and Earth-Spacecraft vector
        angular_separation = acosd(dot(r_spacecraft_current, r_earth_sun_vec) / (norm(r_spacecraft_current) * norm(r_earth_sun_vec)));

        % If the angular separation is less than the angular radius, the spacecraft is in shadow
        if angular_separation < angular_radius_earth
            eclipse_status(i) = 0; % Eclipse occurs
        end
    end
end

% Plotting eclipse status around Earth
time_days = t / 86400; % Convert seconds to days
figure;
plot(time_days, eclipse_status, 'LineWidth', 2);
ylabel('Eclipse Status');
yticklabels({'Eclipse', 'No Eclipse'});
yticks([0 1]);
xlabel('Time (days)');
title('Eclipse Window during Spacecraft Orbit around Earth (GTO)');
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
plot(days_array(1:8), eclipse_hours(1:8), '-o', 'LineWidth', 2);
xlabel('Day Number');
ylabel('Eclipse Duration (hours)');
title('Eclipse Duration per Day for Spacecraft Orbit around Earth (GTO)');
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
