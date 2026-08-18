% sim_two_body.m
% Simulates the Keplaria orbit by calling the dynamics of the
% 2-body-problem together with the RK4-integrator.

clear; clc; close all;

%% Load Constants
constants;

%% Simulation Parameters
dt = 10;                % timestep [s]
n_steps = round(T/dt);  % one orbit [-]
t = (0:n_steps-1)*dt;   % time vector [s]

%% Initial conditions (circular robit in equatorial plane)

% position [m]
rx0 = a;                % perigee on x-axis
ry0 = 0;
rz0 = 0;

% velocity [m/s]
v_c = sqrt(mu/a);       % circular velocity [m/s]

vx0 = 0;
vy0 = v_c;
vz0 = 0;

%% Initial state as 1x6 vector
x0 = [rx0; ry0; rz0; vx0; vy0; vz0];

%% Integration

f = @(X) two_body(X, mu);
X = rk4_integrator(f, x0, dt, n_steps);

%% Plot
%figure;
%plot3(X(1,:)/1000, X(2,:)/1000, X(3,:)/1000, 'b', 'LineWidth', 1.5);
%hold on;

% Draw Earth as a sphere
%[xs, ys, zs] = sphere(50);
%surf(xs*R_earth/1000, ys*R_earth/1000, zs*R_earth/1000, ...
%    'FaceColor', [0.2 0.5 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5);

%axis equal; grid on;
%xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]');
%title('Two-Body Orbit (ECI Frame)');
%legend('Orbit', 'Earth');

%% Export to JSON
metadata = struct(...
    'a_km',         a / 1000, ...
    'e',            e, ...
    'i_deg',        rad2deg(i), ...
    'period_min',   T / 60, ...
    'altitude_km',  alt / 1000, ...
    'dt',           dt, ...
    'n_steps',      n_steps);

data = struct(...
    'metadata', metadata, ...
    't',        t, ...
    'X',        X');      % transpose: N x 6

json_str = jsonencode(data, 'PrettyPrint', true);

fid = fopen('exports/scenarios/two_body_example.json', 'w');
fprintf(fid, '%s', json_str);
fclose(fid);

disp('Simulation complete. JSON exported.');