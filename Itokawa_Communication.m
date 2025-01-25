clc; clear;

%% Constants and Initial Parameters
mu_earth = 3.986004418e5; % Gravitational parameter of Earth [km^3/s^2]
mu_sun = 1.33e11; % Gravitational parameter of the Sun [km^3/s^2]
mu_itokawa = 2.33e-9; % Gravitational parameter of Itokawa [km^3/s^2]
conv_rads = pi / 180; % Degrees to radians
start_date = datetime(2037, 2, 8);
end_date = datetime(2039, 8, 7);
timestep = 600; % Ten minute step size [s]
time_vector = start_date:timestep/seconds(days(1)):end_date;

% Orbital Elements of Spacecraft (Starting at Perigee)
kep_spacecraft = [1.165, 0, 29 * conv_rads, 0, 0, 0]; % Orbital elements [a, e, i, Om, om, theta]
[r_spacecraft, v_spacecraft] = kep2cart(kep_spacecraft, mu_itokawa); % Initial Cartesian state of spacecraft

if any(isnan(r_spacecraft)) || any(isnan(v_spacecraft))
    error('Initial spacecraft position or velocity contains NaN values. Check the keplerian conversion.');
end

% Ground station location (James Weir Building, University of Strathclyde)
lat_JWGS = 55.8614 * conv_rads; % Latitude in radians
long_JWGS = -4.2446 * conv_rads; % Longitude in radians
alt_JWGS = 0.1; % Altitude above sea level [km]
JWGS = [lat_JWGS, long_JWGS, alt_JWGS];

% Initialize results storage
comm_window_hours = zeros(length(time_vector), 1);

%% Calculate Communication Windows
for t_idx = 1:length(time_vector)
    current_time = time_vector(t_idx);
    mjd_time = juliandate(current_time) - 2400000.5; % Convert datetime to Modified Julian Date (MJD)
    % Get Keplerian elements of Itokawa and Earth at current time
    kep_itokawa = Itokawa_Ephemeris(mjd_time);
    kep_earth = Earth_Ephemeris(mjd_time);

    % Convert Itokawa's and Earth's Keplerian elements to Cartesian coordinates
    [r_itokawa, ~] = kep2cart(kep_itokawa, mu_sun);
    [r_earth, ~] = kep2cart(kep_earth, mu_sun);

    % Propagate Spacecraft's orbit around Itokawa using Runge-Kutta integration
    [r_spacecraft, v_spacecraft] = runge_kutta_4(@(t, y) orbit_dynamics(t, y, mu_itokawa), timestep, [r_spacecraft; v_spacecraft]);

    % Check if propagation generated NaN values
    if any(isnan(r_spacecraft)) || any(isnan(v_spacecraft))
        % warning('Propagation generated NaN values at time %s. Resetting to initial conditions.', datestr(current_time));
        [r_spacecraft, v_spacecraft] = kep2cart(kep_spacecraft, mu_itokawa);
        continue;
    end

    r_spacecraft_abs = r_itokawa + r_spacecraft;

    % Calculate JWGS position in Earth's rotating frame
    theta_earth = 2 * pi * (seconds(current_time - start_date) / 86400); % Earth's rotation angle
    theta_earth = mod(theta_earth, 2 * pi); % Keep within 0 to 2*pi
    JWGS_rotated = JWGS_Rotation(JWGS, theta_earth);
    [x_JWGS, y_JWGS, z_JWGS] = latlong2cart(JWGS_rotated);
    r_JWGS = [x_JWGS, y_JWGS, z_JWGS];

    % Check Line of Sight (LOS) by verifying that Itokawa and the Sun do not block the signal
    los_vector = r_spacecraft_abs - r_JWGS;
    los_distance = norm(los_vector);

    % Check for Blockages (Itokawa, Sun, Earth)
    sun_blockage = Check_Blockage_Sun(r_spacecraft_abs, r_earth);
    itokawa_blockage = Check_Blockage_Itokawa(r_spacecraft_abs, r_itokawa, r_JWGS);
    earth_blockage = Check_Blockage_Earth(r_spacecraft_abs, r_earth, r_JWGS);

    % Determine if there is LOS and no blockages
    if (~sun_blockage && ~itokawa_blockage && ~earth_blockage && los_distance > 0)
        comm_window_hours(t_idx) = 1; % Communication is possible in this timestep
    end
end

% Convert to daily communication hours
total_timesteps_per_day = floor(24 * 3600 / timestep);
num_days = floor(length(comm_window_hours) / total_timesteps_per_day);
comm_window_hours_daily = reshape(comm_window_hours(1:num_days * total_timesteps_per_day), total_timesteps_per_day, num_days);
comm_window_hours_per_day = sum(comm_window_hours_daily, 1) * (timestep / 3600);

%% Plot Communication Windows
figure;
days_in_orbit = 0:(length(comm_window_hours_per_day) - 1);
plot(days_in_orbit, comm_window_hours_per_day, 'b');
xlabel('Days in Orbit');
ylabel('Communication Hours per Day');
title('Communication Windows During Orbit Around Itokawa');
grid on;

%% Functions
function JWGS_rotated = JWGS_Rotation(JWGS, theta)
    % Convert latitude, longitude, and altitude to Cartesian coordinates and rotate
    lat = JWGS(1);
    lon = JWGS(2) + theta;
    JWGS_rotated = [lat, lon, JWGS(3)];
end

function [x, y, z] = latlong2cart(JWGS)
    % Converts latitude, longitude, and altitude to Cartesian coordinates
    R_earth = 6378.1; % Earth's radius [km]
    lat = JWGS(1);
    lon = JWGS(2);
    alt = JWGS(3);
    x = (R_earth + alt) * cos(lat) * cos(lon);
    y = (R_earth + alt) * cos(lat) * sin(lon);
    z = (R_earth + alt) * sin(lat);
end

function blockage = Check_Blockage_Sun(spacecraft_pos, earth_pos)
    % Check if the Sun is blocking the line of sight between the spacecraft and the Earth
    sun_pos = [0, 0, 0]; % Assuming Sun is at origin
    sun_direction = sun_pos - spacecraft_pos;
    earth_direction = earth_pos - spacecraft_pos;
    cos_angle = dot(sun_direction, earth_direction) / (norm(sun_direction) * norm(earth_direction));
    angle_between = acos(cos_angle);
    sun_radius = 696340; % Sun radius [km]
    earth_sun_distance = norm(earth_pos);
    sun_apparent_radius = atan(sun_radius / earth_sun_distance);
    blockage = (cos_angle > 0) && (angle_between < sun_apparent_radius);
end

function blockage = Check_Blockage_Itokawa(spacecraft_pos, itokawa_pos, jwgs_pos)
    % Check if Itokawa is blocking the line of sight between the spacecraft and JWGS
    los_vector = jwgs_pos - spacecraft_pos;
    itokawa_vector = itokawa_pos - spacecraft_pos;
    cos_angle = dot(los_vector, itokawa_vector) / (norm(los_vector) * norm(itokawa_vector));
    angle_between = acos(cos_angle);
    itokawa_radius = 0.165; % Approximate radius of Itokawa [km]
    distance_to_itokawa = norm(itokawa_vector);
    itokawa_apparent_radius = atan(itokawa_radius / distance_to_itokawa);
    % Blockage occurs if Itokawa is along the LOS and its apparent size covers the angle
    blockage = (cos_angle > 0) && (angle_between < itokawa_apparent_radius);
end

function blockage = Check_Blockage_Earth(spacecraft_pos, earth_pos, jwgs_pos)
    % Check if the Earth is blocking the line of sight between the spacecraft and JWGS
    los_vector = jwgs_pos - spacecraft_pos;
    earth_vector = earth_pos - spacecraft_pos;
    cos_angle = dot(los_vector, earth_vector) / (norm(los_vector) * norm(earth_vector));
    angle_between = acos(cos_angle);
    earth_radius = 6378.1; % Earth's radius [km]
    distance_to_earth = norm(earth_vector);
    earth_apparent_radius = atan(earth_radius / distance_to_earth);
    % Blockage occurs if Earth is along the LOS and its apparent size covers the angle
    blockage = (cos_angle > 0) && (angle_between < earth_apparent_radius);
end

function [r_next, v_next] = runge_kutta_4(odefun, dt, y0)
    % Runge-Kutta 4th order integration for numerical propagation
    k1 = odefun(0, y0);
    k2 = odefun(0, y0 + 0.5 * dt * k1);
    k3 = odefun(0, y0 + 0.5 * dt * k2);
    k4 = odefun(0, y0 + dt * k3);
    y_next = y0 + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
    r_next = y_next(1:3);
    v_next = y_next(4:6);
end

function dydt = orbit_dynamics(~, y, mu)
    % Orbital dynamics for Runge-Kutta propagation
    r = y(1:3);
    v = y(4:6);
    r_norm = norm(r);
    a = -mu * r / r_norm^3;
    dydt = [v; a];
end
