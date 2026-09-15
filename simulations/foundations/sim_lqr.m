% sim_lqr.m
% Simulates LQR-controlled rendezvous in the Hill frame based on CW
% equations

clear; clc;

%% constants
constants;

%% Simulation Parameters
dt = 10;
t_final = T;
n_steps = round(t_final / dt);
t = (0:n_steps-1) * dt;

%% CW Equations on cont. State-space
% with X=[x; y; z; x_dot; y_dot; z_dot] and simple u = [ux; uy; uz]

Ac = [0         0         0         1         0         0  ;
      0         0         0         0         1         0  ;
      0         0         0         0         0         1  ;
      3*n^2     0         0         0         2*n       0  ;
      0         0         0        -2*n       0         0  ;
      0         0        -n^2       0         0         0];

Bc = [zeros(3,3); eye(3)];

%% Discretize (ZOH)

[Ad, Bd] = discretize(Ac, Bc, dt);


%% LQR design

Q = diag([1e-2, 1e-2, 1e-2, ...        % position weights
          1e0, 1e0, 1e0]);          % velocity weights

R = diag([1e4, 1e4, 1e4]);        % control weights

[K, P] = lqr_controller(Ad, Bd, Q, R);

fprintf('LQR gain K computed.\n');
fprintf('Closed-loop eigenvalues:\n');
disp(abs(eig(Ad - Bd*K)));

%% Initial Conditions

x0     =  500;
y0     = 5000;
z0     =  200;
x_dot0 =    0;
y_dot0 =    0;
z_dot0 =    0;

X0 = [x0; y0; z0; x_dot0; y_dot0; z_dot0];

%% Integrate target orbit (for JSON export and visualization)

rx0_t = a; ry0_t = 0; rz0_t = 0;
vc    = sqrt(mu / a);
X_t0  = [rx0_t; ry0_t; rz0_t; 0; vc; 0];
f_target = @(X) two_body(X, mu);
X_t = rk4_integrator(f_target, X_t0, dt, n_steps);

%% Closed Loop Simulation
X_hist = zeros(6, n_steps);   % state history
U_hist = zeros(3, n_steps);   % control history

X_hist(:, 1) = X0;

for k = 1:n_steps-1
    u = -K * X_hist(:, k);
    X_hist(:, k+1) = Ad * X_hist(:, k) + Bd * u;
    U_hist(:, k)   = u;
end
U_hist(:, end) = -K * X_hist(:, end);

%% Plots

t_min = t / 60;   % convert to minutes

c_x   = [0.00 0.60 0.90];   % blue: x / ux
c_y   = [0.90 0.40 0.00];   % orange: y / uy
c_z   = [0.20 0.80 0.20];   % green: z / uz

fig = figure('Name', 'LQR Rendezvous');

% --- Plot 1: Position states ---
ax1 = subplot(2, 3, 1);
hold(ax1, 'on'); grid(ax1, 'on');
plot(ax1, t_min, X_hist(1,:), 'Color', c_x, 'LineWidth', 1.2, 'DisplayName', 'x (radial)');
plot(ax1, t_min, X_hist(2,:), 'Color', c_y, 'LineWidth', 1.2, 'DisplayName', 'y (along-track)');
plot(ax1, t_min, X_hist(3,:), 'Color', c_z, 'LineWidth', 1.2, 'DisplayName', 'z (cross-track)');
xlabel(ax1, 'Time [min]');
ylabel(ax1, 'Position [m]');
title(ax1, 'Position States');
legend(ax1, 'Location', 'northeast');

% --- Plot 2: Velocity states ---
ax2 = subplot(2, 3, 2);
hold(ax2, 'on'); grid(ax2, 'on');
plot(ax2, t_min, X_hist(4,:), 'Color', c_x, 'LineWidth', 1.2, 'DisplayName', '$\dot{x}$');
plot(ax2, t_min, X_hist(5,:), 'Color', c_y, 'LineWidth', 1.2, 'DisplayName', '$\dot{y}$');
plot(ax2, t_min, X_hist(6,:), 'Color', c_z, 'LineWidth', 1.2, 'DisplayName', '$\dot{z}$');
xlabel(ax2, 'Time [min]');
ylabel(ax2, 'Velocity [m/s]');
title(ax2, 'Velocity States');
legend(ax2, 'Location', 'northeast', 'Interpreter', 'latex');

% --- Plot 3: Relative distance ---
ax3 = subplot(2, 3, 3);
hold(ax3, 'on'); grid(ax3, 'on');
rel_dist = vecnorm(X_hist(1:3, :), 2, 1);
plot(ax3, t_min, rel_dist, 'Color', c_x, 'LineWidth', 1.2);
xlabel(ax3, 'Time [min]');
ylabel(ax3, 'Distance [m]');
title(ax3, 'Relative Distance');

% --- Plot 4: Control inputs ---
ax4 = subplot(2, 3, 4);
hold(ax4, 'on'); grid(ax4, 'on');
plot(ax4, t_min, U_hist(1,:)*1000, 'Color', c_x, 'LineWidth', 1.2, 'DisplayName', '$u_x$');
plot(ax4, t_min, U_hist(2,:)*1000, 'Color', c_y, 'LineWidth', 1.2, 'DisplayName', '$u_y$');
plot(ax4, t_min, U_hist(3,:)*1000, 'Color', c_z, 'LineWidth', 1.2, 'DisplayName', '$u_z$');
xlabel(ax4, 'Time [min]');
ylabel(ax4, 'Acceleration [mm/s$^2$]', 'Interpreter', 'latex');
title(ax4, 'Control Inputs');
legend(ax4, 'Location', 'northeast', 'Interpreter', 'latex');

% --- Plot 5: Hill frame trajectory ---
ax5 = subplot(2, 3, 5);
hold(ax5, 'on'); grid(ax5, 'on'); axis(ax5, 'equal');
plot(ax5, X_hist(2,:), X_hist(1,:), 'Color', c_x, 'LineWidth', 1.2);
plot(ax5, X_hist(2,1), X_hist(1,1), 'o', 'Color', c_y, ...
     'MarkerFaceColor', c_y, 'MarkerSize', 6);
plot(ax5, 0, 0, '+k', 'MarkerSize', 8, 'LineWidth', 1.5);
xlabel(ax5, 'Along-track $y$ [m]', 'Interpreter', 'latex');
ylabel(ax5, 'Radial $x$ [m]',      'Interpreter', 'latex');
title(ax5, 'Hill Frame Trajectory');

% --- Plot 6: Cumulative delta-v ---
ax6 = subplot(2, 3, 6);
hold(ax6, 'on'); grid(ax6, 'on');
dv = cumsum(vecnorm(U_hist, 2, 1) * dt);
plot(ax6, t_min, dv, 'Color', c_z, 'LineWidth', 1.2);
xlabel(ax6, 'Time [min]');
ylabel(ax6, '$\Delta v$ [m/s]', 'Interpreter', 'latex');
title(ax6, 'Cumulative $\Delta v$', 'Interpreter', 'latex');

%% To Tikz

addpath('tools/matlab2tikz-master/src');

fig_names = {'lqr_position', 'lqr_velocity', 'lqr_distance', ...
             'lqr_control',  'lqr_trajectory', 'lqr_deltav'};
axes_handles = [ax1, ax2, ax3, ax4, ax5, ax6];

for k = 1:6
    fig_tmp = figure('Visible', 'off');
    ax_tmp  = copyobj(axes_handles(k), fig_tmp);
    ax_tmp.Position = [0.15 0.15 0.75 0.75];
    matlab2tikz(sprintf('results/figures/lqr/%s.tikz', fig_names{k}), ...
        'figurehandle',  fig_tmp, ...
        'width',         '\figurewidth', ...
        'height',        '\figureheight', ...
        'showInfo',      false, ...
        'checkForUpdates', false);
    close(fig_tmp);
end
disp('Tikz files exported.');


%% Performance metrics

% Peak thrust (2-norm across all steps)
thrust_norms          = vecnorm(U_hist, 2, 1);
peak_thrust           = max(thrust_norms);
peak_thrust_step      = find(thrust_norms == peak_thrust, 1);

% Per-axis delta-v
dv_per_axis           = sum(abs(U_hist), 2) * dt;   % [dvx; dvy; dvz]

% Closed-loop eigenvalue
A_cl                  = Ad - Bd * K;
eig_max               = max(abs(eig(A_cl)));

% Convergence (position norm < threshold)
conv_threshold        = 1.0;                         % [m]
rel_dist_full         = vecnorm(X_hist(1:3,:), 2, 1);
conv_idx              = find(rel_dist_full < conv_threshold, 1, 'first');
if ~isempty(conv_idx)
    conv_step         = conv_idx;
    conv_time_s       = t(conv_idx);
    final_pos_error   = X_hist(1:3, conv_idx);
    fprintf('Converged (< %.1fm) at t = %.1f min (step %d)\n', ...
            conv_threshold, conv_time_s/60, conv_step);
else
    conv_step         = NaN;
    conv_time_s       = NaN;
    final_pos_error   = X_hist(1:3, end);
    fprintf('Did not converge within simulation time.\n');
end

fprintf('Peak thrust:    %.6f m/s^2 (step %d, t=%.1f min)\n', ...
        peak_thrust, peak_thrust_step, t(peak_thrust_step)/60);
fprintf('u_max (MPC):    1e-2 m/s^2, LQR exceeds by factor %.1f\n', ...
        peak_thrust / 1e-2);
fprintf('Max CL eigenvalue: %.6f\n', eig_max);
fprintf('Delta-v per axis: dvx=%.2f, dvy=%.2f, dvz=%.2f m/s\n', ...
        dv_per_axis(1), dv_per_axis(2), dv_per_axis(3));

%% JSON export

metadata = struct(...
    'a_km',           a / 1000, ...
    'e',              e, ...
    'i_deg',          rad2deg(i), ...
    'period_min',     T / 60, ...
    'altitude_km',    (a - R_earth) / 1000, ...
    'dt',             dt, ...
    'n_steps',        n_steps, ...
    'controller',     'LQR', ...
    'dynamics_model', 'CWH');

performance = struct(...
    'delta_v_total_m_s',       dv(end), ...
    'delta_v_per_axis_m_s',    dv_per_axis', ...
    'convergence_step',        conv_step, ...
    'convergence_time_s',      conv_time_s, ...
    'convergence_threshold_m', conv_threshold, ...
    'final_position_error_m',  final_pos_error', ...
    'peak_thrust_m_s2',        peak_thrust, ...
    'peak_thrust_step',        peak_thrust_step, ...
    'max_eigenvalue_magnitude', eig_max);

data = struct(...
    'metadata',    metadata, ...
    'performance', performance, ...
    't',           t, ...
    'X_t',         X_t', ...
    'Rho',         X_hist');

json_str = jsonencode(data, 'PrettyPrint', true);
fid = fopen('exports/scenarios/sim_lqr.json', 'w');
fprintf(fid, '%s', json_str);
fclose(fid);

fprintf('Simulation complete. JSON exported.\n');
fprintf('Total delta-v: %.4f m/s\n', dv(end));
