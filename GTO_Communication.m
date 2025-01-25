% Constants
MuE = 3.986e5;  % Gravitational parameter of Earth (Km^3/s^2)
RE = 6378;  % Radius of Earth in km
omega_earth = 7.2921159e-5;  % Earth's angular velocity in rad/s

% Orbit parameters
orbit_start_date = datetime(2036, 3, 16);
orbit_end_date = orbit_start_date + days(8);  % One week period

% Time parameters
time_step = 0.0006944;  % a minute in days
time = orbit_start_date:days(time_step):orbit_end_date;
num_steps = length(time);

% Orbital Elements - GTO
semi_major_axis = 24367.5;  % Semi-major axis in km
eccentricity = 0.728;       % Eccentricity
inclination_deg = 40;       % Inclination in degrees
RAAN_deg = 53.65;           % Right Ascension of Ascending Node in degrees
arg_perigee_deg = 0;        % Argument of Perigee in degrees
true_anomaly_deg = 0;       % True Anomaly at start in degrees

% Convert angles to radians
inclination_rad = deg2rad(inclination_deg);
RAAN_rad = deg2rad(RAAN_deg);
arg_perigee_rad = deg2rad(arg_perigee_deg);
true_anomaly_rad = deg2rad(true_anomaly_deg);

% Calculate initial position and velocity in perifocal coordinates
p = semi_major_axis * (1 - eccentricity^2);  % Semi-latus rectum in km
r_perifocal = p / (1 + eccentricity * cos(true_anomaly_rad));

% Position in perifocal coordinates
r0_perifocal = r_perifocal * [cos(true_anomaly_rad); sin(true_anomaly_rad); 0];

% Velocity in perifocal coordinates
v0_perifocal = sqrt(MuE / p) * [-sin(true_anomaly_rad); eccentricity + cos(true_anomaly_rad); 0];

% Rotation matrices to convert from perifocal to ECI
R3_W = [cos(RAAN_rad), -sin(RAAN_rad), 0; sin(RAAN_rad), cos(RAAN_rad), 0; 0, 0, 1];  % Rotation by RAAN
R1_i = [1, 0, 0; 0, cos(inclination_rad), -sin(inclination_rad); 0, sin(inclination_rad), cos(inclination_rad)];  % Rotation by inclination
R3_w = [cos(arg_perigee_rad), -sin(arg_perigee_rad), 0; sin(arg_perigee_rad), cos(arg_perigee_rad), 0; 0, 0, 1];  % Rotation by argument of perigee

% Complete rotation matrix from perifocal to ECI
Q_p2eci = R3_W * R1_i * R3_w;

% Convert initial position and velocity to ECI coordinates
r0_eci = Q_p2eci * r0_perifocal;
v0_eci = Q_p2eci * v0_perifocal;

% Initial state vector
y0 = double([r0_eci; v0_eci]);  % Initial state vector (position and velocity)

% Integrate the spacecraft trajectory using initial conditions
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
[t, state] = ode45(@(t, y) two_body_equation(t, y, MuE), double([0, days(orbit_end_date - orbit_start_date) * 86400]), y0, options);  % Adjust time span for orbital period

% Extract spacecraft positions over time
r_spacecraft = state(:, 1:3)';

% Ground station coordinates (James Weir Building, University of Strathclyde)
lat_gs = 55.8608;  % Latitude in degrees
lon_gs = -4.2442;  % Longitude in degrees
h_gs = 0.1;  % Altitude in km (approximate)

% Convert ground station coordinates to ECEF (Earth-Centered Earth-Fixed)
lat_gs_rad = deg2rad(lat_gs);
lon_gs_rad = deg2rad(lon_gs);
R_gs_ecef = (RE + h_gs) * [cos(lat_gs_rad) * cos(lon_gs_rad); cos(lat_gs_rad) * sin(lon_gs_rad); sin(lat_gs_rad)];

% Communication window analysis
comm_status = zeros(1, length(t));  % Initially assume no communication

for i = 1:length(t)
    % Update ground station position in ECI frame due to Earth's rotation
    theta = omega_earth * t(i);  % Rotation angle in radians
    R_gs_eci = [cos(theta), -sin(theta), 0; sin(theta), cos(theta), 0; 0, 0, 1] * R_gs_ecef;
    
    % Vector from ground station to spacecraft
    r_gs_sc = r_spacecraft(:, i) - R_gs_eci;
    
    % Check line-of-sight (dot product between ground station position and vector to spacecraft)
    if dot(R_gs_eci, r_gs_sc) > 0
        comm_status(i) = 1;  % Communication possible
    end
end

% Plotting communication status over time
time_days = t / 86400;  % Convert seconds to days
figure;
plot(time_days, comm_status, 'LineWidth', 2);
ylabel('Communication Status');
yticklabels({'No Communication', 'Communication'});
yticks([0 1]);
xlabel('Time (days)');
title('Communication Window between Spacecraft and JWGS');
grid on;

% Calculate communication hours for each day
days_array = unique(floor(time_days));  % Extract unique day numbers
num_days = length(days_array);
comm_hours = zeros(1, num_days);

for i = 1:num_days
    % Find indices corresponding to the current day
    day_indices = find(floor(time_days) == days_array(i));
    % Calculate the total time in communication for the current day
    comm_hours(i) = sum(comm_status(day_indices) == 1) * (time_step * 24);  % Correct time conversion to hours
end

% Plotting communication hours per day
figure;
plot(days_array(1:8), comm_hours(1:8), '-o', 'LineWidth', 2);
xlabel('Day Number');
ylabel('Communication Hours');
title('Communication Hours per Day between Spacecraft and JWGS');
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
