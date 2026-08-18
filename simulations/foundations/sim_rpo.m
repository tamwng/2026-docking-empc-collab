% sim_rpo.m
% Simulates relative motion between target and chaser spacecraft
% in the Hill frame (LVLH) using selectable dynamics models.
%
% DYNAMICS MODELS:
%   'nonlinear' - Exact nonlinear EoMs (Schaub 14.13)
%   'general'   - Linearized, elliptic target orbit (Schaub 14.20)
%   'CWH'       - Clohessy-Wiltshire-Hill, circular target orbit
%
% OUTPUT:
%   - Phase plot of relative trajectory in Hill frame
%   - Export to JSON for GUI visualization

clear all; clc;

%% Configuration (choose CWH/general/nonlinear, to be expanded)

dynamics_model = 'CWH';

%% Constants
constants;

%% Simulation Parameters (as in two_body)
dt      = 10;               % [s]   timestep
n_steps = round(T / dt);    % [-]   one full target orbit
t       = (0:n_steps-1)*dt; % [s]   time vector

%% Target initial condition (circular orbit in equatorial plane)

rx0_t = a;
ry0_t = 0;
rz0_t = 0;
vc    = sqrt(mu / a);           % circular velocity
vx0_t = 0;
vy0_t = vc;
vz0_t = 0;

X_t0 = [rx0_t; ry0_t; rz0_t; vx0_t; vy0_t; vz0_t];


%% Chaser initial conditions in Hill frame

x0     =  - 1000e3;   % [m]   radial offset
y0     =  1000e3;   % [m]   along-track offset
z0     =   500e3;   % [m]   cross-track offset
x_dot0 =    0;   % [m/s]
y_dot0 = -2 * n * x0;  % [m/s]  drift-free condition, bound relative orbit constrained
z_dot0 =  n * z0;   % [m/s]

Rho0 = [x0; y0; z0; x_dot0; y_dot0; z_dot0];

%% Integrate Target Orbit (like in two_body.m)

f_target = @(X) two_body(X, mu);
X_t = rk4_integrator(f_target, X_t0, dt, n_steps);

%% Integration of Relative Motion

% Switch Cases between the dynamic models
switch dynamics_model
    case 'nonlinear'
        f_rel = @(Rho, X_t_k) relative_motion_nonlinear(Rho, X_t_k, mu);
    case 'general'
        f_rel = @(Rho, X_t_k) relative_motion_general(Rho, X_t_k, mu);
    case 'CWH'
        f_rel = @(Rho, ~) clohessy_wiltshire(Rho, n);
    otherwise
        error('Unknown dynamics model: %s', dynamics_model);
end

Rho = zeros(6, n_steps);
Rho(:, 1) = Rho0;

for k = 1:n_steps-1
    X_t_k = X_t(:, k);
    f     = @(rho) f_rel(rho, X_t_k);
    xk    = Rho(:, k);
    k1    = f(xk);
    k2    = f(xk + 0.5*dt*k1);
    k3    = f(xk + 0.5*dt*k2);
    k4    = f(xk + dt*k3);
    Rho(:, k+1) = xk + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
end

%% Plots

figure;
subplot(1,2,1);
plot(Rho(2,:)/1000, Rho(1,:)/1000, 'b', 'LineWidth', 1.5);
xlabel('Along-track y [km]');
ylabel('Radial x [km]');
title(sprintf('Relative Trajectory (%s)', dynamics_model));
grid on; axis equal;

subplot(1,2,2);
plot3(Rho(1,:)/1000, Rho(2,:)/1000, Rho(3,:)/1000, 'b', 'LineWidth', 1.5);
xlabel('Radial x [km]');
ylabel('Along-track y [km]');
zlabel('Cross-track z [km]');
title('3D Relative Trajectory');
grid on; axis equal;

%% Export to JSON
metadata = struct(...
    'a_km',           a / 1000, ...
    'e',              e, ...
    'i_deg',          rad2deg(i), ...
    'period_min',     T / 60, ...
    'altitude_km',    (a - R_earth) / 1000, ...
    'dt',             dt, ...
    'n_steps',        n_steps, ...
    'dynamics_model', dynamics_model); 

data = struct(...
    'metadata', metadata, ...
    't',        t, ...
    'X_t',      X_t', ...     % target ECI state, N x 6
    'Rho',      Rho');        % relative state in Hill frame, N x 6

json_str = jsonencode(data, 'PrettyPrint', true);

fid = fopen('exports/scenarios/sim_rpo.json', 'w');
fprintf(fid, '%s', json_str);
fclose(fid);

save('exports/scenarios/sim_rpo.mat', 'Rho', 't', 'n', 'metadata');
disp('Simulation complete. JSON exported.');