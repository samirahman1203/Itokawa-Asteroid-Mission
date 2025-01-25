% Constants
mu_sun = 1.33e+11; % Gravitational parameter of the Sun (km^3/s^2)
au_to_km = 1.496e+8; % Conversion factor from au to km

%% Outward Journey

% Define Launch Date - March 16, 2036
target_date = datetime(2036, 3, 16);

% Convert target date to Modified Julian Date (MJD)
target_mjd = juliandate(target_date) - 2400000.5;

% Get Keplerian Elements for Earth and Itokawa using ephemeris functions
kep_earth = Earth_Ephemeris(target_mjd);
kep_itokawa = Itokawa_Ephemeris(target_mjd);

% Convert Keplerian elements to Cartesian coordinates for current positions
[r_earth, ~] = kep2cart(kep_earth, mu_sun); % Earth's position and velocity in Cartesian coordinates
[r_itokawa, ~] = kep2cart(kep_itokawa, mu_sun); % Itokawa's position and velocity in Cartesian coordinates

% Transfer Calculation using current positions
r1 = norm(r_earth); % Radius of Earth's orbit (km)
r2 = norm(r_itokawa); % Radius of Itokawa's orbit (km)

% Semi-major axis of the transfer orbit
a_transfer = (r1 + r2) / 2; 

% Eccentricity of the transfer orbit
e_transfer = (r2 - r1) / (r1 + r2);

% Update Itokawa's Position After Transfer Time
elapsed_days = 370; % Days taken for the transfer
elapsed_seconds = elapsed_days * 86400; % Convert to seconds

% Update the Mean Anomaly for Itokawa after elapsed time
kep_itokawa_arrival = kep_itokawa;
n_itokawa = sqrt(mu_sun / (kep_itokawa(1)^3)); % Mean motion of Itokawa
M0_itokawa = kep_itokawa(6);
M_arrival = M0_itokawa + n_itokawa * elapsed_seconds; % Mean anomaly at arrival

% Solve Kepler's Equation for the new position of Itokawa
eccAnomaly = M_arrival; % Initial guess for Eccentric Anomaly
tolerance = 1e-8;
while true
    E_next = M_arrival + kep_itokawa(2) * sin(eccAnomaly);
    if abs(E_next - eccAnomaly) < tolerance
        break;
    end
    eccAnomaly = E_next;
end

% Calculate True Anomaly
wom_arrival = 2 * atan2(sqrt(1 + kep_itokawa(2)) * sin(E_next / 2), sqrt(1 - kep_itokawa(2)) * cos(E_next / 2));
kep_itokawa_arrival(6) = wom_arrival;

% Get Itokawa's position at arrival in Cartesian coordinates
[r_itokawa_arrival, ~] = kep2cart(kep_itokawa_arrival, mu_sun);

% Angle from the Sun to Earth's initial position
initial_angle = atan2(r_earth(2), r_earth(1));

% Angle from the Sun to the position of Itokawa at arrival
arrival_angle = atan2(r_itokawa_arrival(2), r_itokawa_arrival(1));

% Calculate true anomaly at which the transfer orbit intersects with Itokawa's arrival
theta_transfer_end = arrival_angle - initial_angle;

% Create an array of true anomalies for the transfer orbit, from 0 to theta_transfer_end
theta_transfer = linspace(0, theta_transfer_end, 500);
r_transfer = (a_transfer * (1 - e_transfer^2)) ./ (1 + e_transfer * cos(theta_transfer));

% Generate transfer orbit in polar coordinates, then rotate to align with Earth's current position
x_transfer = r_transfer .* cos(theta_transfer + initial_angle);
y_transfer = r_transfer .* sin(theta_transfer + initial_angle);
z_transfer = zeros(size(x_transfer)); % Assume transfer orbit is in the ecliptic plane

% Define number of points for plotting Earth's and Itokawa's orbits
theta_values = linspace(0, 2*pi, 1000);

% Preallocate arrays for Earth's Cartesian coordinates
x_earth_orbit = zeros(size(theta_values));
y_earth_orbit = zeros(size(theta_values));
z_earth_orbit = zeros(size(theta_values));

% Preallocate arrays for Itokawa's Cartesian coordinates
x_itokawa_orbit = zeros(size(theta_values));
y_itokawa_orbit = zeros(size(theta_values));
z_itokawa_orbit = zeros(size(theta_values));

% Loop through all true anomalies for Earth and Itokawa
for k = 1:length(theta_values)
    % Set current true anomaly for plotting full orbits
    true_anomaly = theta_values(k);
    
    % Update Earth Keplerian elements with current true anomaly
    kep_earth_orbit = kep_earth;
    kep_earth_orbit(6) = true_anomaly;
    
    % Update Itokawa Keplerian elements with current true anomaly
    kep_itokawa_orbit = kep_itokawa;
    kep_itokawa_orbit(6) = true_anomaly;
    
    % Convert to Cartesian coordinates
    [r_earth_orbit, ~] = kep2cart(kep_earth_orbit, mu_sun);
    [r_itokawa_orbit, ~] = kep2cart(kep_itokawa_orbit, mu_sun);
    
    % Store Cartesian coordinates for Earth and Itokawa
    x_earth_orbit(k) = r_earth_orbit(1);
    y_earth_orbit(k) = r_earth_orbit(2);
    z_earth_orbit(k) = r_earth_orbit(3);
    
    x_itokawa_orbit(k) = r_itokawa_orbit(1);
    y_itokawa_orbit(k) = r_itokawa_orbit(2);
    z_itokawa_orbit(k) = r_itokawa_orbit(3);
end

% Plotting
figure;
hold on;
axis equal;

% Plot the Sun
plot3(0, 0, 0, 'yo', 'MarkerFaceColor', 'y', 'MarkerSize', 12);
text(0, 0, 0, ' Sun', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

% Plot Earth's orbit
plot3(x_earth_orbit, y_earth_orbit, z_earth_orbit, 'b');
text(x_earth_orbit(end), y_earth_orbit(end), z_earth_orbit(end), ' Earth Orbit', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

% Plot Itokawa's orbit
plot3(x_itokawa_orbit, y_itokawa_orbit, z_itokawa_orbit, 'g');
text(x_itokawa_orbit(end), y_itokawa_orbit(end), z_itokawa_orbit(end), ' Itokawa Orbit', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

% Plot Transfer orbit
plot3(x_transfer, y_transfer, z_transfer, 'r', 'LineWidth', 1.5);
text(x_transfer(end), y_transfer(end), z_transfer(end), ' Lambert Transfer', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

% Add markers for Earth and Itokawa at the start and arrival of the transfer
plot3(r_earth(1), r_earth(2), r_earth(3), 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 8);
text(r_earth(1), r_earth(2), r_earth(3), ' Earth (March 16, 2036)', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

plot3(r_itokawa_arrival(1), r_itokawa_arrival(2), r_itokawa_arrival(3), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
text(r_itokawa_arrival(1), r_itokawa_arrival(2), r_itokawa_arrival(3), ' Itokawa (Arrival)', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

% Labels and title
xlabel('X (10^6 km)');
ylabel('Y (10^6 km)');
zlabel('Z (10^6 km)');
title('Lambert Transfer from Earth to Itokawa');

grid on;
hold off;

%% Return Journey 

% Define Launch Date - August 7, 2039
departure_date = datetime(2039, 8, 7);

% Convert departure date to Modified Julian Date (MJD)
departure_mjd = juliandate(departure_date) - 2400000.5;

% Get Keplerian Elements for Itokawa and Earth using ephemeris functions
kep_itokawa = Itokawa_Ephemeris(departure_mjd);
kep_earth = Earth_Ephemeris(departure_mjd);

% Convert Keplerian elements to Cartesian coordinates for current positions
[r_itokawa, v_itokawa] = kep2cart(kep_itokawa, mu_sun); % Itokawa's position and velocity in Cartesian coordinates
[r_earth, v_earth] = kep2cart(kep_earth, mu_sun); % Earth's position and velocity in Cartesian coordinates

% Transfer Calculation using current positions
r1 = norm(r_itokawa); % Radius of Itokawa's orbit (km)
r2 = norm(r_earth); % Radius of Earth's orbit (km)

% Semi-major axis of the transfer orbit
a_transfer = (r1 + r2) / 2; 

% Eccentricity of the transfer orbit
e_transfer = (r2 - r1) / (r1 + r2);

% Update Earth's Position After Transfer Time
elapsed_days = 240; % Days taken for the transfer
elapsed_seconds = elapsed_days * 86400; % Convert to seconds

% Update the Mean Anomaly for Earth after elapsed time
kep_earth_arrival = kep_earth;
n_earth = sqrt(mu_sun / (kep_earth(1)^3)); % Mean motion of Earth
M0_earth = kep_earth(6);
M_arrival = M0_earth + n_earth * elapsed_seconds; % Mean anomaly at arrival

% Solve Kepler's Equation for the new position of Earth
eccAnomaly = M_arrival; % Initial guess for Eccentric Anomaly
tolerance = 1e-8;
while true
    E_next = M_arrival + kep_earth(2) * sin(eccAnomaly);
    if abs(E_next - eccAnomaly) < tolerance
        break;
    end
    eccAnomaly = E_next;
end

% Calculate True Anomaly
wom_arrival = 2 * atan2(sqrt(1 + kep_earth(2)) * sin(E_next / 2), sqrt(1 - kep_earth(2)) * cos(E_next / 2));
kep_earth_arrival(6) = wom_arrival;

% Get Earth's position at arrival in Cartesian coordinates
[r_earth_arrival, ~] = kep2cart(kep_earth_arrival, mu_sun);

% Angle from the Sun to Itokawa's initial position
initial_angle = atan2(r_itokawa(2), r_itokawa(1));

% Angle from the Sun to the position of Earth at arrival
arrival_angle = atan2(r_earth_arrival(2), r_earth_arrival(1));

% Calculate true anomaly at which the transfer orbit intersects with Earth's arrival
theta_transfer_end = arrival_angle - initial_angle;

% Create an array of true anomalies for the transfer orbit, from 0 to theta_transfer_end
theta_transfer = linspace(0, theta_transfer_end, 500);
r_transfer = (a_transfer * (1 - e_transfer^2)) ./ (1 + e_transfer * cos(theta_transfer));

% Generate transfer orbit in polar coordinates, then rotate to align with Itokawa's current position
x_transfer = r_transfer .* cos(theta_transfer + initial_angle);
y_transfer = r_transfer .* sin(theta_transfer + initial_angle);
z_transfer = zeros(size(x_transfer)); % Assume transfer orbit is in the ecliptic plane

% Define number of points for plotting Itokawa's and Earth's orbits
theta_values = linspace(0, 2*pi, 1000);

% Preallocate arrays for Itokawa's Cartesian coordinates
x_itokawa_orbit = zeros(size(theta_values));
y_itokawa_orbit = zeros(size(theta_values));
z_itokawa_orbit = zeros(size(theta_values));

% Preallocate arrays for Earth's Cartesian coordinates
x_earth_orbit = zeros(size(theta_values));
y_earth_orbit = zeros(size(theta_values));
z_earth_orbit = zeros(size(theta_values));

% Loop through all true anomalies for Itokawa and Earth
for k = 1:length(theta_values)
    % Set current true anomaly for plotting full orbits
    true_anomaly = theta_values(k);
    
    % Update Itokawa Keplerian elements with current true anomaly
    kep_itokawa_orbit = kep_itokawa;
    kep_itokawa_orbit(6) = true_anomaly;
    
    % Update Earth Keplerian elements with current true anomaly
    kep_earth_orbit = kep_earth;
    kep_earth_orbit(6) = true_anomaly;
    
    % Convert to Cartesian coordinates
    [r_itokawa_orbit, ~] = kep2cart(kep_itokawa_orbit, mu_sun);
    [r_earth_orbit, ~] = kep2cart(kep_earth_orbit, mu_sun);
    
    % Store Cartesian coordinates for Itokawa and Earth
    x_itokawa_orbit(k) = r_itokawa_orbit(1);
    y_itokawa_orbit(k) = r_itokawa_orbit(2);
    z_itokawa_orbit(k) = r_itokawa_orbit(3);
    
    x_earth_orbit(k) = r_earth_orbit(1);
    y_earth_orbit(k) = r_earth_orbit(2);
    z_earth_orbit(k) = r_earth_orbit(3);
end

% Plotting
figure;
hold on;
axis equal;

% Plot the Sun
plot3(0, 0, 0, 'yo', 'MarkerFaceColor', 'y', 'MarkerSize', 12);
text(0, 0, 0, ' Sun', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

% Plot Itokawa's orbit
plot3(x_itokawa_orbit, y_itokawa_orbit, z_itokawa_orbit, 'g');
text(x_itokawa_orbit(end), y_itokawa_orbit(end), z_itokawa_orbit(end), ' Itokawa Orbit', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

% Plot Earth's orbit
plot3(x_earth_orbit, y_earth_orbit, z_earth_orbit, 'b');
text(x_earth_orbit(end), y_earth_orbit(end), z_earth_orbit(end), ' Earth Orbit', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

% Plot Transfer orbit
plot3(x_transfer, y_transfer, z_transfer, 'r', 'LineWidth', 1.5);
text(x_transfer(end), y_transfer(end), z_transfer(end), ' Lambert Transfer', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

% Add markers for Itokawa and Earth at the start and arrival of the transfer
plot3(r_itokawa(1), r_itokawa(2), r_itokawa(3), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
text(r_itokawa(1), r_itokawa(2), r_itokawa(3), ' Itokawa (August 7, 2039)', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

plot3(r_earth_arrival(1), r_earth_arrival(2), r_earth_arrival(3), 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 8);
text(r_earth_arrival(1), r_earth_arrival(2), r_earth_arrival(3), ' Earth (Arrival)', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

% Labels and title
xlabel('X (10^6 km)');
ylabel('Y (10^6 km)');
zlabel('Z (10^6 km)');
title('Lambert Transfer from Itokawa to Earth');

grid on;
hold off;
