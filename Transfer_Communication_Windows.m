%% MATLAB Script for Communication Window Analysis between JWGS and Spacecraft to Itokawa Asteroid and Return Journey

% Constants
mu_sun = 1.33e+11; % Gravitational parameter of the Sun in km^3/s^2
launch_date = datetime(2036, 3, 16, 0, 0, 0); % Launch date
transfer_duration = 370; % Transfer duration in days
time_step = 60; % Time step in seconds (1 minute)
total_steps = transfer_duration * 24 * 60; % Total time steps (minutes in 370 days)
Earth_rot_period = 86164; % Earth's rotation period in seconds
radius_earth = 6378; % Radius of Earth in km

% Ground station location (James Weir building, University of Strathclyde)
lat_JWGS = 55.8607; % Latitude in degrees
lon_JWGS = -4.2444; % Longitude in degrees
alt_JWGS = 0.1; % Altitude in km (rough estimate)

% Convert latitude and longitude to radians
lat_JWGS_rad = deg2rad(lat_JWGS);
lon_JWGS_rad = deg2rad(lon_JWGS);

% Preallocate array for communication hours per day (Outbound)
comm_hours_per_day_outbound = zeros(1, transfer_duration);

% Define Lambert solver options
optionsLMR = 0;

% Calculate the initial and final position of Earth and Itokawa (Outbound)
initial_time = datenum(launch_date);
final_time = initial_time + transfer_duration;
kep_earth_initial = Earth_Ephemeris(initial_time);
kep_itokawa_final = Itokawa_Ephemeris(final_time);
[r_earth_initial, ~] = kep2cart(kep_earth_initial, mu_sun);
[r_itokawa_final, ~] = kep2cart(kep_itokawa_final, mu_sun);

% Solve Lambert's problem to get the trajectory of the spacecraft (Outbound)
[~, ~, ~, ~, v_earth_initial, v_sc_final, ~, ~] = lambertMR(r_earth_initial, r_itokawa_final, transfer_duration * 86400, mu_sun, 0, 0, 0, optionsLMR);

% Ensure initial position and velocity are column vectors
r_sc = r_earth_initial(:); % Initial spacecraft position as column vector
v_sc = v_earth_initial(:); % Initial velocity from Lambert's solution as column vector
dt = time_step; % Time step in seconds

for step = 1:total_steps
    % Update position and velocity using Runge-Kutta 4th Order Method (RK4)
    k1_r = v_sc;
    k1_v = -mu_sun * r_sc / norm(r_sc)^3;

    k2_r = v_sc + 0.5 * dt * k1_v;
    k2_v = -mu_sun * (r_sc + 0.5 * dt * k1_r) / norm(r_sc + 0.5 * dt * k1_r)^3;

    k3_r = v_sc + 0.5 * dt * k2_v;
    k3_v = -mu_sun * (r_sc + 0.5 * dt * k2_r) / norm(r_sc + 0.5 * dt * k2_r)^3;

    k4_r = v_sc + dt * k3_v;
    k4_v = -mu_sun * (r_sc + dt * k3_r) / norm(r_sc + dt * k3_r)^3;

    r_sc = r_sc + (dt / 6) * (k1_r + 2 * k2_r + 2 * k3_r + k4_r);
    v_sc = v_sc + (dt / 6) * (k1_v + 2 * k2_v + 2 * k3_v + k4_v);
    
    % Ensure r_sc remains a column vector
    r_sc = r_sc(:);
    v_sc = v_sc(:);
    
    % Calculate JWGS position in ECEF
    theta_G = 2 * pi * mod((step * time_step), Earth_rot_period) / Earth_rot_period;
    x_JWGS = (alt_JWGS + radius_earth) * cos(lat_JWGS_rad) * cos(lon_JWGS_rad + theta_G);
    y_JWGS = (alt_JWGS + radius_earth) * cos(lat_JWGS_rad) * sin(lon_JWGS_rad + theta_G);
    z_JWGS = (alt_JWGS + radius_earth) * sin(lat_JWGS_rad);
    r_JWGS = [x_JWGS; y_JWGS; z_JWGS];
    
    % Ensure r_JWGS is a column vector
    r_JWGS = r_JWGS(:);
    
    % Check line-of-sight condition at every timestep
    r_sc_JWGS = r_sc - r_JWGS; % Vector from JWGS to spacecraft

    % Check if the spacecraft is above the horizon
    % The condition for line-of-sight is that the dot product between r_sc_JWGS
    % and r_JWGS should be greater than zero (i.e., the spacecraft is above the horizon)
    if dot(r_sc_JWGS, r_JWGS) > 0
        % Determine which day of the transfer this step belongs to
        day_idx = ceil(step / (24 * 60));
        % Increment communication hours for the respective day
        comm_hours_per_day_outbound(day_idx) = comm_hours_per_day_outbound(day_idx) + (time_step / 3600);
    end
end

% Return journey parameters
departure_return_date = datetime(2039, 8, 7, 0, 0, 0); % Departure from Itokawa
duration_return = 240; % Return duration in days
total_steps_return = duration_return * 24 * 60; % Total time steps (minutes in 240 days)

% Preallocate array for communication hours per day (Return)
comm_hours_per_day_return = zeros(1, duration_return);

% Calculate the initial and final position of Itokawa and Earth (Return)
initial_time_return = datenum(departure_return_date);
final_time_return = initial_time_return + duration_return;
kep_itokawa_initial = Itokawa_Ephemeris(initial_time_return);
kep_earth_final = Earth_Ephemeris(final_time_return);
[r_itokawa_initial, ~] = kep2cart(kep_itokawa_initial, mu_sun);
[r_earth_final, ~] = kep2cart(kep_earth_final, mu_sun);

% Solve Lambert's problem to get the trajectory of the spacecraft (Return)
[~, ~, ~, ~, v_itokawa_initial, v_earth_final, ~, ~] = lambertMR(r_itokawa_initial, r_earth_final, duration_return * 86400, mu_sun, 0, 0, 0, optionsLMR);

% Ensure initial position and velocity are column vectors
r_sc_return = r_itokawa_initial(:); % Initial spacecraft position as column vector
v_sc_return = v_itokawa_initial(:); % Initial velocity from Lambert's solution as column vector
dt = time_step; % Time step in seconds

for step = 1:total_steps_return
    % Update position and velocity using Runge-Kutta 4th Order Method (RK4)
    k1_r = v_sc_return;
    k1_v = -mu_sun * r_sc_return / norm(r_sc_return)^3;

    k2_r = v_sc_return + 0.5 * dt * k1_v;
    k2_v = -mu_sun * (r_sc_return + 0.5 * dt * k1_r) / norm(r_sc_return + 0.5 * dt * k1_r)^3;

    k3_r = v_sc_return + 0.5 * dt * k2_v;
    k3_v = -mu_sun * (r_sc_return + 0.5 * dt * k2_r) / norm(r_sc_return + 0.5 * dt * k2_r)^3;

    k4_r = v_sc_return + dt * k3_v;
    k4_v = -mu_sun * (r_sc_return + dt * k3_r) / norm(r_sc_return + dt * k3_r)^3;

    r_sc_return = r_sc_return + (dt / 6) * (k1_r + 2 * k2_r + 2 * k3_r + k4_r);
    v_sc_return = v_sc_return + (dt / 6) * (k1_v + 2 * k2_v + 2 * k3_v + k4_v);
    
    % Ensure r_sc_return remains a column vector
    r_sc_return = r_sc_return(:);
    v_sc_return = v_sc_return(:);
    
    % Calculate JWGS position in ECEF
    theta_G = 2 * pi * mod((step * time_step), Earth_rot_period) / Earth_rot_period;
    x_JWGS = (alt_JWGS + radius_earth) * cos(lat_JWGS_rad) * cos(lon_JWGS_rad + theta_G);
    y_JWGS = (alt_JWGS + radius_earth) * cos(lat_JWGS_rad) * sin(lon_JWGS_rad + theta_G);
    z_JWGS = (alt_JWGS + radius_earth) * sin(lat_JWGS_rad);
    r_JWGS = [x_JWGS; y_JWGS; z_JWGS];
    
    % Ensure r_JWGS is a column vector
    r_JWGS = r_JWGS(:);
    
    % Check line-of-sight condition at every timestep
    r_sc_JWGS = r_sc_return - r_JWGS; % Vector from JWGS to spacecraft

    % Check if the spacecraft is above the horizon
    % The condition for line-of-sight is that the dot product between r_sc_JWGS
    % and r_JWGS should be greater than zero (i.e., the spacecraft is above the horizon)
    if dot(r_sc_JWGS, r_JWGS) > 0
        % Determine which day of the transfer this step belongs to
        day_idx = ceil(step / (24 * 60));
        % Increment communication hours for the respective day
        comm_hours_per_day_return(day_idx) = comm_hours_per_day_return(day_idx) + (time_step / 3600);
    end
end

% Plot the number of communication hours per day for outbound journey
figure;
plot(1:transfer_duration, comm_hours_per_day_outbound, 'LineWidth', 2);
xlabel('Day of Transfer');
ylabel('Communication Hours');
title('Communication Hours between Spacecraft and JWGS (Outbound)');
grid on;

% Plot the number of communication hours per day for return journey
figure;
plot(1:duration_return, comm_hours_per_day_return, 'LineWidth', 2);
xlabel('Day of Return Transfer');
ylabel('Communication Hours');
title('Communication Hours between Spacecraft and JWGS (Return)');
grid on;
