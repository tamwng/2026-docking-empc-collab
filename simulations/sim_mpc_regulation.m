function [X_out, U_out, t_out, dv_out] = sim_mpc_regulation(params)
% sim_mpc_regulation.m
% Regulation MPC rendezvous in the Hill frame (CWH equations).
%
% Standalone call:  sim_mpc_regulation()
%   Original scenario with plots and JSON export.
%
% Comparison call:  [X, U, t, dv] = sim_mpc_regulation(params)
%   params fields: x0, x_switch, dt, N, u_max, t_final
%   Returns trajectory without plots; uses comparison-mode cost weights.

standalone = nargin < 1;
if standalone; clc; end


%%  LOAD CONSTANTS
constants;


%%  PARAMETERS
if standalone
    dt       = 60;
    t_final  = 2 * T;
    N        = 20;
    x0       = [500; 10000; 200; 0; 0; 0];
    x_switch = zeros(6,1);
    u_max    = 1e-2;
    conv_tol = 1.0;
else
    dt       = params.dt;
    t_final  = params.t_final;
    N        = params.N;
    x0       = params.x0;
    x_switch = params.x_switch;
    u_max    = params.u_max;
    conv_tol = 50;
end

n_steps = round(t_final / dt);
t_all   = (0:n_steps) * dt;


%%  CWH STATE SPACE
Ac = [0      0     0    1     0    0  ;
      0      0     0    0     1    0  ;
      0      0     0    0     0    1  ;
      3*n^2  0     0    0     2*n  0  ;
      0      0     0   -2*n   0    0  ;
      0      0    -n^2  0     0    0 ];
Bc = [zeros(3,3); eye(3)];
[Ad, Bd] = discretize(Ac, Bc, dt);


%%  COST MATRICES
if standalone
    Q = diag([1e-2, 1e-2, 1e-2, 1e0, 1e0, 1e0]);
    R = diag([1e4,  1e4,  1e4]);
else
    % Comparison mode: match EMPC's R so fuel figures are comparable.
    % Q penalises position error; tune as needed.
    Q = 1e-6 * eye(6);
    R = eye(3);
end
[~, P] = lqr_controller(Ad, Bd, Q, R);


%%  CONSTRAINTS
con.u_max        = u_max;
con.y_min_active = standalone;
con.los_cone.active = standalone;
if standalone
    con.los_cone.half_angle = deg2rad(20);
    con.los_cone.n_faces    = 10;
end


%%  SIMULATION LOOP
X_hist = zeros(6, n_steps + 1);
U_hist = zeros(3, n_steps + 1);
X_hist(:,1) = x0;
e_k   = x0 - x_switch;
n_run = 1;

if ~standalone
    fprintf('  Running MPC (comparison)... ');
else
    fprintf('Running Regulation MPC (%d steps, N=%d)...\n', n_steps, N);
end

for k = 1:n_steps
    if norm(e_k) < conv_tol
        n_run = k;
        break;
    end
    [u_opt, ~, ~, ~, ~] = mpc_regulation(e_k, Ad, Bd, Q, R, P, N, con);
    U_hist(:,k)   = u_opt;
    X_hist(:,k+1) = Ad * X_hist(:,k) + Bd * u_opt;
    e_k           = X_hist(:,k+1) - x_switch;
    n_run         = k + 1;
end

X_hist   = X_hist(:, 1:n_run);
U_hist   = U_hist(:, 1:n_run);
dv_total = sum(vecnorm(U_hist(:,1:end-1), 2, 1)) * dt;

if ~standalone
    fprintf('done. Δv = %.3f m/s\n', dv_total);
end

X_out  = X_hist;
U_out  = U_hist;
t_out  = t_all(1:n_run);
dv_out = dv_total;

if nargout > 0; return; end


%%  STANDALONE: POST-PROCESSING
fprintf('Simulation complete.\n');

t_x      = t_all(1:n_run)   / 60;   % [min] for state plots
t_u      = t_all(1:n_run-1) / 60;   % [min] for control plots
rel_dist = vecnorm(X_hist(1:3,:), 2, 1);

conv_idx = find(rel_dist < conv_tol, 1, 'first');
if ~isempty(conv_idx)
    fprintf('Converged at t = %.1f min (step %d)\n', t_x(conv_idx), conv_idx);
else
    fprintf('Did not converge within simulation time — increase t_final or tune Q/R\n');
end


%%  PLOTS (matlab2tikz compatible)
c_x = [0.00 0.60 0.90];
c_y = [0.90 0.40 0.00];
c_z = [0.20 0.80 0.20];

fig = figure('Name', 'Regulation MPC Rendezvous');

ax1 = subplot(2, 3, 1);
hold(ax1, 'on'); grid(ax1, 'on');
plot(ax1, t_x, X_hist(1,:), 'Color', c_x, 'LineWidth', 1.2, 'DisplayName', 'x (radial)');
plot(ax1, t_x, X_hist(2,:), 'Color', c_y, 'LineWidth', 1.2, 'DisplayName', 'y (along-track)');
plot(ax1, t_x, X_hist(3,:), 'Color', c_z, 'LineWidth', 1.2, 'DisplayName', 'z (cross-track)');
xlabel(ax1, 'Time [min]'); ylabel(ax1, 'Position [m]');
title(ax1, 'Position States'); legend(ax1, 'Location', 'northeast');

ax2 = subplot(2, 3, 2);
hold(ax2, 'on'); grid(ax2, 'on');
plot(ax2, t_x, X_hist(4,:), 'Color', c_x, 'LineWidth', 1.2, 'DisplayName', '$\dot{x}$');
plot(ax2, t_x, X_hist(5,:), 'Color', c_y, 'LineWidth', 1.2, 'DisplayName', '$\dot{y}$');
plot(ax2, t_x, X_hist(6,:), 'Color', c_z, 'LineWidth', 1.2, 'DisplayName', '$\dot{z}$');
xlabel(ax2, 'Time [min]'); ylabel(ax2, 'Velocity [m/s]');
title(ax2, 'Velocity States');
legend(ax2, 'Location', 'northeast', 'Interpreter', 'latex');

ax3 = subplot(2, 3, 3);
hold(ax3, 'on'); grid(ax3, 'on');
plot(ax3, t_x, rel_dist, 'Color', c_x, 'LineWidth', 1.2);
xlabel(ax3, 'Time [min]'); ylabel(ax3, 'Distance [m]');
title(ax3, 'Relative Distance');

ax4 = subplot(2, 3, 4);
hold(ax4, 'on'); grid(ax4, 'on');
plot(ax4, t_u, U_hist(1,1:end-1)*1000, 'Color', c_x, 'LineWidth', 1.2, 'DisplayName', '$u_x$');
plot(ax4, t_u, U_hist(2,1:end-1)*1000, 'Color', c_y, 'LineWidth', 1.2, 'DisplayName', '$u_y$');
plot(ax4, t_u, U_hist(3,1:end-1)*1000, 'Color', c_z, 'LineWidth', 1.2, 'DisplayName', '$u_z$');
xlabel(ax4, 'Time [min]');
ylabel(ax4, 'Acceleration [mm/s$^2$]', 'Interpreter', 'latex');
title(ax4, 'Control Inputs');
legend(ax4, 'Location', 'northeast', 'Interpreter', 'latex');
yline(ax4,  con.u_max*1000, '--k', 'LineWidth', 0.8, 'DisplayName', 'u\_max');
yline(ax4, -con.u_max*1000, '--k', 'LineWidth', 0.8);

ax5 = subplot(2, 3, 5);
hold(ax5, 'on'); grid(ax5, 'on'); axis(ax5, 'equal');
plot(ax5, X_hist(2,:), X_hist(1,:), 'Color', c_x, 'LineWidth', 1.2);
plot(ax5, X_hist(2,1), X_hist(1,1), 'o', 'Color', c_y, ...
     'MarkerFaceColor', c_y, 'MarkerSize', 6);
plot(ax5, 0, 0, '+k', 'MarkerSize', 8, 'LineWidth', 1.5);
xlabel(ax5, 'Along-track $y$ [m]', 'Interpreter', 'latex');
ylabel(ax5, 'Radial $x$ [m]',      'Interpreter', 'latex');
title(ax5, 'Hill Frame Trajectory');

ax6 = subplot(2, 3, 6);
hold(ax6, 'on'); grid(ax6, 'on');
dv_cum = cumsum(vecnorm(U_hist(:,1:end-1), 2, 1)) * dt;
plot(ax6, t_u, dv_cum, 'Color', c_z, 'LineWidth', 1.2);
xlabel(ax6, 'Time [min]');
ylabel(ax6, '$\Delta v$ [m/s]', 'Interpreter', 'latex');
title(ax6, 'Cumulative $\Delta v$', 'Interpreter', 'latex');

fprintf('Total delta-v: %.4f m/s\n', dv_total);


%%  EXPORT JSON
X_t0     = [a; 0; 0; 0; sqrt(mu/a); 0];
f_target = @(Xv) two_body(Xv, mu);
X_t      = rk4_integrator(f_target, X_t0, dt, n_run - 1);

metadata = struct(...
    'a_km',           a / 1000, ...
    'e',              e, ...
    'i_deg',          rad2deg(i), ...
    'period_min',     T / 60, ...
    'altitude_km',    (a - R_earth) / 1000, ...
    'dt',             dt, ...
    'n_steps',        n_steps, ...
    'controller',     'MPC_regulation', ...
    'dynamics_model', 'CWH', ...
    'horizon_N',      N, ...
    'u_max',          con.u_max, ...
    'y_min_active',   con.y_min_active);

data = struct(...
    'metadata', metadata, ...
    't',        t_all(1:n_run), ...
    'X_t',      X_t', ...
    'Rho',      X_hist');

json_str = jsonencode(data, 'PrettyPrint', true);
fid = fopen('exports/scenarios/sim_mpc_regulation.json', 'w');
fprintf(fid, '%s', json_str);
fclose(fid);

fprintf('JSON exported.\n');
