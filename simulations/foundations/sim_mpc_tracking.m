%sim_mpc_tracking.m
% Simulates tracking MPC rednezvous in the Hill frame by extending the
% regulation_mpc with a reference trajectory to track. A hold point is used
% for gnerality. For the start, we do not want to dock but rather track a
% holdpoint behind the target

clear; clc;

%% Constants
constants;


%%  SIMULATION PARAMETERS
dt      = 60;
t_final = 2 * T;            
n_steps = round(t_final / dt);
t       = (0:n_steps-1) * dt;
t_min   = t / 60;


%%  MPC PARAMETERS
N     = 20;       % [-]  prediction horizon

%%  CONSTRAINTS
con.u_max        = 1e-2;   % [m/s^2]  symmetric thrust bound
con.y_min_active = true;   % enforce y >= 0: approach only from behind


%%  CWH STATE SPACE
Ac = [0      0     0    1     0    0  ;
      0      0     0    0     1    0  ;
      0      0     0    0     0    1  ;
      3*n^2  0     0    0     2*n  0  ;
      0      0     0   -2*n   0    0  ;
      0      0    -n^2  0     0    0 ];

Bc = [zeros(3,3); eye(3)];


%%  DISCRETIZE
[Ad, Bd] = discretize(Ac, Bc, dt);


%%  COST MATRICES (as before for comparison)
Q = diag([1e-2, 1e-2, 1e-2, 1e0, 1e0, 1e0]);
R = diag([1e4,  1e4,  1e4]);

[~, P] = lqr_controller(Ad, Bd, Q, R);


%%  INITIAL CONDITIONS

x0_ic = 500;    % [m] radial offset
y0_ic = 10000;   % [m] along-track offset
z0_ic = 0;      % [m] cross-track offset

X0 = [x0_ic; y0_ic; z0_ic; 0; 0; 0];


%%  HOLD POINT
x_hold = 0;     % [m] radial   — on V-bar axis
y_hold = 100;   % [m] along-track — 100m behind target
z_hold = 0;     % [m] cross-track

% LoS cone approach corridor — apex at hold point, docking axis along +y.
% Set active = true to enforce.  Hard constraint: infeasible if initial
% state violates the cone.  Check: sqrt(x0^2+z0^2) <= (y0-y_hold)*tan(half_angle).
con.los_cone.active     = false;
con.los_cone.half_angle = deg2rad(20);                % [rad]  20° sensor FOV
con.los_cone.n_faces    = 10;                         % [-]    inner polygon faces
con.los_cone.apex       = [x_hold; y_hold; z_hold];  % [m]    docking port position


%%  REFERENCE TRAJECTORY

traj_mode = 'rbar_vbar';         % 'vbar' | 'rbar' | 'rbar_vbar'

% Reference ends at 70% of simulation time
n_ref            = round(0.7 * n_steps);

[R_ref_short, phase_switch_idx] = reference_trajectory(x0_ic, y0_ic, z0_ic, ...
                                                        n_ref, dt, traj_mode, ...
                                                        x_hold, y_hold, z_hold);

% Pad remaining steps with hold point — MPC regulates to hold point
R_ref_pad = repmat([x_hold; y_hold; z_hold; 0; 0; 0], 1, n_steps - n_ref);
R_ref     = [R_ref_short, R_ref_pad];

fprintf('Trajectory mode:   %s\n', traj_mode);
fprintf('Hold point:        x=%.0fm, y=%.0fm, z=%.0fm\n', x_hold, y_hold, z_hold);
fprintf('Reference ends at: %.1f min (step %d)\n', t_min(n_ref), n_ref);
fprintf('Hold phase:        %.1f min remaining\n', t_min(n_steps) - t_min(n_ref));
if ~isempty(phase_switch_idx)
    fprintf('Phase switch at:   %.1f min (step %d)\n', ...
            t_min(phase_switch_idx), phase_switch_idx);
end



%%  INTEGRATE TARGET ORBIT
X_t0     = [a; 0; 0; 0; sqrt(mu/a); 0];
f_target = @(X) two_body(X, mu);
X_t      = rk4_integrator(f_target, X_t0, dt, n_steps);


%%  CLOSED-LOOP TRACKING MPC SIMULATION
X_hist = zeros(6, n_steps);
U_hist = zeros(3, n_steps);
E_hist = zeros(6, n_steps);   % tracking error history

X_hist(:, 1) = X0;
E_hist(:, 1) = X0 - R_ref(:, 1);

fprintf('Running Tracking MPC (%d steps, N=%d)...\n', n_steps, N);

for k = 1:n_steps-1

    % N-step reference preview: columns k+1 ... k+N, padded with hold point
    k_end      = min(k + N, n_steps);
    R_preview  = [R_ref(:, k+1:k_end), ...
                  repmat(R_ref(:, end), 1, k + N - k_end)];

    % solve tracking MPC
    [u_opt, ~, ~] = mpc_tracking(X_hist(:,k), R_preview, Ad, Bd, Q, R, P, N, con);

    % apply control and propagate
    U_hist(:, k)   = u_opt;
    X_hist(:, k+1) = Ad * X_hist(:,k) + Bd * u_opt;
    E_hist(:, k+1) = X_hist(:, k+1) - R_ref(:, k+1);
end
R_preview_end = repmat(R_ref(:, end), 1, N);
U_hist(:,end) = mpc_tracking(X_hist(:,end), R_preview_end, Ad, Bd, Q, R, P, N, con);

fprintf('Simulation complete.\n');


%%  POST-PROCESSING

rel_dist  = vecnorm(X_hist(1:3,:), 2, 1);
track_err = vecnorm(E_hist(1:3,:), 2, 1);

% Lyapunov on error state
V_hist = zeros(1, n_steps);
for k = 1:n_steps
    V_hist(k) = E_hist(:,k)' * P * E_hist(:,k);
end

% Convergence to hold point
tol      = 1.0;   % [m]
hold_dist = vecnorm(X_hist(1:3,:) - [x_hold; y_hold; z_hold], 2, 1);
conv_idx  = find(hold_dist < tol, 1, 'first');
if ~isempty(conv_idx)
    fprintf('Reached hold point at t = %.1f min (step %d)\n', ...
            t_min(conv_idx), conv_idx);
else
    fprintf('Did not reach hold point within simulation time.\n');
end

% Delta-v
dv = cumsum(vecnorm(U_hist, 2, 1) * dt);
fprintf('Total delta-v: %.4f m/s\n', dv(end));

% Final state diagnostics
fprintf('Final state:      x=%.4f, y=%.4f, z=%.4f\n', ...
        X_hist(1,end), X_hist(2,end), X_hist(3,end));
fprintf('Hold point:       x=%.4f, y=%.4f, z=%.4f\n', x_hold, y_hold, z_hold);
fprintf('Final error:      ex=%.4f, ey=%.4f, ez=%.4f\n', ...
        E_hist(1,end), E_hist(2,end), E_hist(3,end));
fprintf('Final velocity:   vx=%.4f, vy=%.4f, vz=%.4f\n', ...
        X_hist(4,end), X_hist(5,end), X_hist(6,end));



%%  PLOTS (matlab2tikz compatible)

c_x = [0.00 0.60 0.90];
c_y = [0.90 0.40 0.00];
c_z = [0.20 0.80 0.20];
c_r = [0.80 0.00 0.80];

%% Figure 1 — States and Control
fig1 = figure('Name', 'Tracking MPC — States and Control');

ax1 = subplot(2, 2, 1);
hold(ax1, 'on'); grid(ax1, 'on');
plot(ax1, t_min, X_hist(1,:), 'Color', c_x, 'LineWidth', 1.2, 'DisplayName', 'x');
plot(ax1, t_min, X_hist(2,:), 'Color', c_y, 'LineWidth', 1.2, 'DisplayName', 'y');
plot(ax1, t_min, X_hist(3,:), 'Color', c_z, 'LineWidth', 1.2, 'DisplayName', 'z');
plot(ax1, t_min, R_ref(1,:),  '--', 'Color', c_x, 'LineWidth', 0.8, 'DisplayName', 'x ref');
plot(ax1, t_min, R_ref(2,:),  '--', 'Color', c_y, 'LineWidth', 0.8, 'DisplayName', 'y ref');
if ~isempty(phase_switch_idx)
    xline(ax1, t_min(phase_switch_idx), ':', 'Color', [0.5 0.5 0.5], ...
          'LineWidth', 1.0, 'DisplayName', 'Phase switch');
end
xline(ax1, t_min(n_ref), '--', 'Color', [0.4 0.4 0.4], ...
      'LineWidth', 1.0, 'DisplayName', 'Hold phase');
xlabel(ax1, 'Time [min]'); ylabel(ax1, 'Position [m]');
title(ax1, 'Position States');
legend(ax1, 'Location', 'northeast', 'FontSize', 6);

ax2 = subplot(2, 2, 2);
hold(ax2, 'on'); grid(ax2, 'on');
plot(ax2, t_min, X_hist(4,:), 'Color', c_x, 'LineWidth', 1.2, 'DisplayName', '$\dot{x}$');
plot(ax2, t_min, X_hist(5,:), 'Color', c_y, 'LineWidth', 1.2, 'DisplayName', '$\dot{y}$');
plot(ax2, t_min, X_hist(6,:), 'Color', c_z, 'LineWidth', 1.2, 'DisplayName', '$\dot{z}$');
if ~isempty(phase_switch_idx)
    xline(ax2, t_min(phase_switch_idx), ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.0);
end
xline(ax2, t_min(n_ref), '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.0);
xlabel(ax2, 'Time [min]'); ylabel(ax2, 'Velocity [m/s]');
title(ax2, 'Velocity States');
legend(ax2, 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 6);

ax3 = subplot(2, 2, 3);
hold(ax3, 'on'); grid(ax3, 'on');
plot(ax3, t_min, U_hist(1,:)*1000, 'Color', c_x, 'LineWidth', 1.2, 'DisplayName', '$u_x$');
plot(ax3, t_min, U_hist(2,:)*1000, 'Color', c_y, 'LineWidth', 1.2, 'DisplayName', '$u_y$');
plot(ax3, t_min, U_hist(3,:)*1000, 'Color', c_z, 'LineWidth', 1.2, 'DisplayName', '$u_z$');
yline(ax3,  con.u_max*1000, '--k', 'LineWidth', 0.8, 'DisplayName', '$u_{max}$');
yline(ax3, -con.u_max*1000, '--k', 'LineWidth', 0.8);
if ~isempty(phase_switch_idx)
    xline(ax3, t_min(phase_switch_idx), ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.0);
end
xline(ax3, t_min(n_ref), '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.0);
xlabel(ax3, 'Time [min]');
ylabel(ax3, 'Acceleration [mm/s$^2$]', 'Interpreter', 'latex');
title(ax3, 'Control Inputs');
legend(ax3, 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 6);

ax4 = subplot(2, 2, 4);
hold(ax4, 'on'); grid(ax4, 'on');
plot(ax4, t_min, dv, 'Color', c_z, 'LineWidth', 1.2);
if ~isempty(phase_switch_idx)
    xline(ax4, t_min(phase_switch_idx), ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.0);
end
xline(ax4, t_min(n_ref), '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.0);
xlabel(ax4, 'Time [min]');
ylabel(ax4, '$\Delta v$ [m/s]', 'Interpreter', 'latex');
title(ax4, 'Cumulative $\Delta v$', 'Interpreter', 'latex');

%% Figure 2 — Trajectory and Analysis
fig2 = figure('Name', 'Tracking MPC — Trajectory and Analysis');

ax5 = subplot(2, 2, 1);
hold(ax5, 'on'); grid(ax5, 'on'); axis(ax5, 'equal');
plot(ax5, X_hist(2,:), X_hist(1,:), 'Color', c_x, 'LineWidth', 1.2, ...
     'DisplayName', 'Actual');
plot(ax5, R_ref(2,1:n_ref), R_ref(1,1:n_ref), '--', 'Color', c_r, ...
     'LineWidth', 1.0, 'DisplayName', 'Reference');
plot(ax5, X_hist(2,1), X_hist(1,1), 'o', 'Color', c_y, ...
     'MarkerFaceColor', c_y, 'MarkerSize', 6, 'DisplayName', 'Start');
plot(ax5, y_hold, x_hold, 's', 'Color', [0.4 0.4 0.4], ...
     'MarkerFaceColor', [0.4 0.4 0.4], 'MarkerSize', 6, 'DisplayName', 'Hold point');
if ~isempty(phase_switch_idx)
    plot(ax5, R_ref(2,phase_switch_idx), R_ref(1,phase_switch_idx), 's', ...
         'Color', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5], ...
         'MarkerSize', 6, 'DisplayName', 'Phase switch');
end
plot(ax5, 0, 0, '+k', 'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', 'Target');
xlabel(ax5, 'Along-track $y$ [m]', 'Interpreter', 'latex');
ylabel(ax5, 'Radial $x$ [m]',      'Interpreter', 'latex');
title(ax5, 'Hill Frame Trajectory');
legend(ax5, 'Location', 'northeast', 'FontSize', 6);

ax6 = subplot(2, 2, 2);
hold(ax6, 'on'); grid(ax6, 'on');
plot(ax6, t_min, rel_dist,  'Color', c_x, 'LineWidth', 1.2, 'DisplayName', 'Distance to target');
plot(ax6, t_min, hold_dist, 'Color', c_r, 'LineWidth', 1.2, 'DisplayName', 'Distance to hold point');
plot(ax6, t_min, track_err, '--', 'Color', c_y, 'LineWidth', 1.0, 'DisplayName', 'Tracking error');
if ~isempty(phase_switch_idx)
    xline(ax6, t_min(phase_switch_idx), ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.0);
end
xline(ax6, t_min(n_ref), '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.0);
if ~isempty(conv_idx)
    xline(ax6, t_min(conv_idx), '--', 'Color', c_y, 'LineWidth', 1.0, ...
          'DisplayName', sprintf('Hold at %.1f min', t_min(conv_idx)));
end
xlabel(ax6, 'Time [min]'); ylabel(ax6, 'Distance [m]');
title(ax6, 'Distance to Target and Hold Point');
legend(ax6, 'Location', 'northeast', 'FontSize', 6);

ax7 = subplot(2, 2, 3);
hold(ax7, 'on'); grid(ax7, 'on');
plot(ax7, t_min, E_hist(1,:), 'Color', c_x, 'LineWidth', 1.2, 'DisplayName', '$e_x$');
plot(ax7, t_min, E_hist(2,:), 'Color', c_y, 'LineWidth', 1.2, 'DisplayName', '$e_y$');
plot(ax7, t_min, E_hist(3,:), 'Color', c_z, 'LineWidth', 1.2, 'DisplayName', '$e_z$');
if ~isempty(phase_switch_idx)
    xline(ax7, t_min(phase_switch_idx), ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.0);
end
xline(ax7, t_min(n_ref), '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.0);
xlabel(ax7, 'Time [min]'); ylabel(ax7, 'Error [m]');
title(ax7, 'Tracking Error per Axis');
legend(ax7, 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 6);

ax8 = subplot(2, 2, 4);
hold(ax8, 'on'); grid(ax8, 'on');
plot(ax8, t_min, V_hist, 'Color', c_x, 'LineWidth', 1.2);
if ~isempty(phase_switch_idx)
    xline(ax8, t_min(phase_switch_idx), ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.0);
end
xline(ax8, t_min(n_ref), '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.0);
if ~isempty(conv_idx)
    xline(ax8, t_min(conv_idx), '--', 'Color', c_y, 'LineWidth', 1.0);
end
xlabel(ax8, 'Time [min]');
ylabel(ax8, '$V = \mathbf{e}^\top P \mathbf{e}$', 'Interpreter', 'latex');
title(ax8, 'Lyapunov Function (error)');

%%  EXPORT JSON
metadata = struct(...
    'a_km',            a / 1000, ...
    'e',               e, ...
    'i_deg',           rad2deg(i), ...
    'period_min',      T / 60, ...
    'altitude_km',     (a - R_earth) / 1000, ...
    'dt',              dt, ...
    'n_steps',         n_steps, ...
    'controller',      'MPC_tracking', ...
    'dynamics_model',  'CWH', ...
    'trajectory_mode', traj_mode, ...
    'horizon_N',       N, ...
    'u_max',           con.u_max, ...
    'x_hold',          x_hold, ...
    'y_hold',          y_hold, ...
    'z_hold',          z_hold);

data = struct(...
    'metadata', metadata, ...
    't',        t, ...
    'X_t',      X_t', ...
    'Rho',      X_hist', ...
    'R_ref',    R_ref');

json_str = jsonencode(data, 'PrettyPrint', true);
fid = fopen('exports/scenarios/sim_mpc_tracking.json', 'w');
fprintf(fid, '%s', json_str);
fclose(fid);

fprintf('JSON exported.\n');