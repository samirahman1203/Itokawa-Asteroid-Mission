% MATLAB code to perform refined eclipse analysis for spacecraft eclipse
% Including realistic trajectory using Lambert's solver and Ephemeris data

% Constants
MuS = 1.32712440018e11; % Gravitational parameter of the Sun (Km^3/s^2)
au = 1.496e11; % Astronomical Unit in meters

% Launch parameters
launch_date = datetime(2039, 8, 7);
transfer_time_days = 228; % 
time_step = 0.01; % in days
end_date = launch_date + days(transfer_time_days);

% Eclipse parameters
eclipse_angle = 10; % degrees cone angle for eclipse

% Initialise time array
time = launch_date:days(time_step):end_date;
num_steps = length(time);

% Obtain ephemeris data for Earth and Itokawa
r_earth_sun = zeros(3, num_steps);
r_itokawa_sun = zeros(3, num_steps);
for i = 1:num_steps
    current_time = time(i);
    % Get Earth and Itokawa positions from ephemeris (in heliocentric frame)
    ephemeris_data_earth = Earth_Ephemeris('Earth');
    ephemeris_data_itokawa = Itokawa_Ephemeris('Itokawa');
    r_earth_sun(:, i) = ephemeris_data_earth(1:3); % Use first three elements as position vector
    r_itokawa_sun(:, i) = ephemeris_data_itokawa(1:3); % Use first three elements as position vector
end

% Using Lambert's method to calculate the spacecraft trajectory
% Initial and final positions from Earth and Itokawa at start and end dates
r1 = r_itokawa_sun(:, 1);
r2 = r_earth_sun(:, end);
[vi_scalar, vf_scalar] = lambertMR(r1, r2, transfer_time_days * 86400, MuS); % Time in seconds

% Adjust for scalar velocity values if returned by lambertMR
if isscalar(vi_scalar)
    % Assuming scalar represents magnitude and direction can be inferred from r1 to r2
    direction_vector = (r2 - r1) / norm(r2 - r1); % Unit vector in the direction of travel
    vi = vi_scalar * direction_vector; % Assign velocity vector based on direction
else
    vi = vi_scalar; % If not scalar, assume proper vector returned
end

% Ensure velocity vector is properly shaped
y0 = [r1; vi]; % Initial state vector (position and velocity)

if length(y0) ~= 6
    error('State vector must have exactly 6 elements (position and velocity). Please check initial conditions.');
end

% Integrate the spacecraft trajectory using initial conditions from Lambert's solution
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
[t, state] = ode45(@(t, y) two_body_equation(t, y, MuS), [0, transfer_time_days * 86400], y0, options); % Adjust time span for 228-day transfer

% Extract spacecraft positions over time
r_spacecraft = state(:, 1:3)';

% Interpolate ephemeris data to match spacecraft time points
r_earth_sun_interp = interp1(linspace(0, transfer_time_days * 86400, num_steps), r_earth_sun', t)';

% Eclipse analysis
eclipse_status = ones(1, length(t)); % Initially assume eclipse is not occurring

for i = 1:length(t)
    % Vector from Earth to spacecraft
    r_earth = r_earth_sun_interp(:, i); % Get Earth position at the closest available ephemeris time step
    r_spacecraft_current = r_spacecraft(:, i);
    r_earth_spacecraft = r_spacecraft_current - r_earth;

    % Vector from Earth to Sun (negative of Earth's position in heliocentric frame)
    r_earth_sun_vec = -r_earth;

    % Calculate angle between vectors (in degrees)
    cos_theta = dot(r_earth_spacecraft, r_earth_sun_vec) / (norm(r_earth_spacecraft) * norm(r_earth_sun_vec));
    theta = acosd(cos_theta);

    % Additional check to see if spacecraft is in the shadow cone (i.e., properly aligned)
    % Check if spacecraft lies between Earth and Sun
    if cos_theta > 0 && theta < eclipse_angle && norm(r_earth) > norm(r_spacecraft_current) && norm(r_earth_spacecraft) < norm(r_earth)
        eclipse_status(i) = 0; % Eclipse occurring
    end
end

% Plotting eclipse status
time_days = t / 86400; % Convert seconds to days
figure;
plot(time_days, eclipse_status, 'LineWidth', 2);
ylabel('Eclipse Status');
yticklabels({'In Eclipse', 'Not in Eclipse'});
yticks([0 1]);
xlabel('Time (days)');
title('Eclipse Status during Spacecraft Transfer to Earth');

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
plot(days_array, eclipse_hours, 'LineWidth', 2);
xlabel('Day Number');
ylabel('Eclipse Duration (hours)');
title('Eclipse Duration per Day for Spacecraft Transfer to Earth');
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
