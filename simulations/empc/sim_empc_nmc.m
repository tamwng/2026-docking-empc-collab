% sim_empc_nmc.m
% Economic MPC — NMC orbit acquisition and maintenance.
%
% SCENARIO:
% The chaser sits at the V-bar hold point (0, 100, 0) m with zero velocity
% after the tracking MPC phase.  The NMC orbit (rho=50m, 2:1 ellipse) passes
% exactly through (0, 100, 0) at phase phi_0 = -pi/2 — so the chaser is at
% the right position but with zero velocity. The EMPC finds the minimum-fuel injection
% impulse and then maintains the orbit with near-zero thrust (free drift).
%
% ECONOMIC COST:
% Stage cost: ℓ(u) = u'R_eco u 
% Terminal cost: e_N'P e_N  where e_N = x(N) - r_nmc(k+N)
%


clear; clc;


%%  CONSTANTS 
constants;

%% Simulation params
dt      = 60;
t_final = 3 * T;
n_steps = round(t_final / dt);
t       = (0:n_steps-1) * dt;
t_min   = t / 60;

%% MPC params
N = 20;  


%%  NMC ORBIT PARAMETERS

rho   = 50;          % [m]    radial semi-axis  (semi-minor of Hill ellipse)
y_c   = 0;           % [m]    along-track drift centre
phi_0 = -pi/2;       % [rad]  initial phase chosen so r_nmc(0) = (0, 100, 0, rho*n, 0, 0)
                     %        — position matches hold point, only velocity differs

% Full NMC reference
X_nmc = nmc_trajectory(rho, phi_0, n, dt, n_steps + N, y_c);

% Verification: NMC at step 0 should be at (0, 100, 0)
assert(abs(X_nmc(1,1)) < 1e-9 && abs(X_nmc(2,1) - 100) < 1e-9, ...
       'Phase phi_0 does not align NMC with hold point — check phi_0.');


%%  CWH STATE SPACE

Ac = [0      0     0    1     0    0  ;
      0      0     0    0     1    0  ;
      0      0     0    0     0    1  ;
      3*n^2  0     0    0     2*n  0  ;
      0      0     0   -2*n   0    0  ;
      0      0    -n^2  0     0    0 ];
Bc = [zeros(3,3); eye(3)];

[Ad, Bd] = discretize(Ac, Bc, dt);


%%  INITIAL CONDITION  — hold point, zero velocity

X0 = [0; 100; 0; 0; 0; 0];


%%  TERMINAL COST

Q_dare = diag([1e-2, 1e-2, 1e-2, 1e0, 1e0, 1e0]);
R_dare = diag([1e4,  1e4,  1e4]);
[~, P] = lqr_controller(Ad, Bd, Q_dare, R_dare);


%%  COST MATRICES

Q_stage = zeros(6);   % no state penalty at intermediate steps
R_eco   = eye(3);     


%%  CONSTRAINTS

con.u_max        = 1e-2;          % [m/s^2]
con.y_min_active = false;         % no overshoot constraint for NMC
con.los_cone.active = false;
con.los_cone.half_angle = deg2rad(20);
con.los_cone.n_faces    = 10;
con.u_ss         = zeros(3, 1);   % NMC is zero-input — bypass feedforward


%%  CLOSED-LOOP EMPC SIMULATION

X_hist = zeros(6, n_steps);
U_hist = zeros(3, n_steps);
X_hist(:, 1) = X0;

fprintf('Running NMC EMPC (%d steps, N=%d, rho=%.0fm)...\n', n_steps, N, rho);

for k = 1:n_steps-1
    R_preview = X_nmc(:, k+1 : k+N);   % phase-synchronised N-step preview
    [u_opt, ~, ~] = mpc_tracking(X_hist(:,k), R_preview, Ad, Bd, ...
                                  Q_stage, R_eco, P, N, con);
    U_hist(:, k)   = u_opt;
    X_hist(:, k+1) = Ad * X_hist(:,k) + Bd * u_opt;
end
R_preview_end = X_nmc(:, n_steps+1 : n_steps+N);
U_hist(:, end) = mpc_tracking(X_hist(:,end), R_preview_end, Ad, Bd, ...
                               Q_stage, R_eco, P, N, con);

fprintf('Simulation complete.\n');


%%  POST-PROCESSING

% NMC position error (how far off the orbit we are)
E_nmc     = X_hist - X_nmc(:, 1:n_steps);
err_pos   = vecnorm(E_nmc(1:3,:), 2, 1);
err_vel   = vecnorm(E_nmc(4:6,:), 2, 1);

% Drift-free condition: ydot + 2*n*x = 0 on the NMC orbit
drift_cond = abs(X_hist(5,:) + 2*n*X_hist(1,:));

% Cumulative delta-v
dv = cumsum(vecnorm(U_hist, 2, 1) * dt);

% One full NMC ellipse for plotting (reference shape)
phi_plot = linspace(0, 2*pi, 500);
x_ell    =  rho * cos(phi_plot);
y_ell    =  y_c - 2*rho * sin(phi_plot);

fprintf('\nTotal delta-v:          %.4f m/s\n', dv(end));
fprintf('Initial NMC error:      %.2f m  (pos)  %.4f m/s (vel)\n', err_pos(1), err_vel(1));
fprintf('Final NMC error:        %.2f m  (pos)  %.4f m/s (vel)\n', err_pos(end), err_vel(end));
fprintf('Final drift condition:  %.2e m/s\n', drift_cond(end));

% Estimate injection ΔV (burn in first quarter orbit ~ 23 steps)
n_inject = round(0.25 * T / dt);
fprintf('Injection delta-v:      %.4f m/s  (first %.0f min)\n', ...
        dv(n_inject), t_min(n_inject));
fprintf('Maintenance delta-v:    %.4f m/s  (remaining)\n', ...
        dv(end) - dv(n_inject));


%%  PLOTS

c_traj = [0.00 0.45 0.70];   % blue  — actual trajectory
c_nmc  = [0.85 0.33 0.10];   % orange — NMC reference
c_x    = [0.00 0.60 0.90];
c_y    = [0.90 0.40 0.00];
c_z    = [0.20 0.80 0.20];

% Orbital period marker positions
T_marks = t_min(round((1:3) * T/dt));

fig1 = figure('Name', 'EMPC NMC — Orbit Acquisition');

%% Hill frame trajectory
ax1 = subplot(2,3,1);
hold(ax1,'on'); grid(ax1,'on'); axis(ax1,'equal');
plot(ax1, y_ell, x_ell, '--', 'Color', c_nmc, 'LineWidth', 1.0, 'DisplayName', 'NMC orbit');
plot(ax1, X_hist(2,:), X_hist(1,:), 'Color', c_traj, 'LineWidth', 1.4, 'DisplayName', 'Actual');
plot(ax1, X0(2), X0(1), 'o', 'Color', c_y, 'MarkerFaceColor', c_y, ...
     'MarkerSize', 6, 'DisplayName', 'Start');
plot(ax1, 0, 0, '+k', 'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', 'Target');
xlabel(ax1, 'Along-track $y$ [m]', 'Interpreter', 'latex');
ylabel(ax1, 'Radial $x$ [m]',      'Interpreter', 'latex');
title(ax1, 'Hill Frame Trajectory');
legend(ax1, 'Location', 'northeast', 'FontSize', 7);

%% NMC position error
ax2 = subplot(2,3,2);
hold(ax2,'on'); grid(ax2,'on');
semilogy(ax2, t_min, err_pos, 'Color', c_traj, 'LineWidth', 1.4, 'DisplayName', 'Position error');
semilogy(ax2, t_min, err_vel * 1000, '--', 'Color', c_nmc, 'LineWidth', 1.2, ...
         'DisplayName', 'Velocity error ×1000');
for tm = T_marks
    xline(ax2, tm, ':', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8);
end
xlabel(ax2, 'Time [min]'); ylabel(ax2, 'Error [m  or  mm/s]');
title(ax2, 'NMC Error (log scale)');
legend(ax2, 'Location', 'northeast', 'FontSize', 7);

%% Drift-free condition
ax3 = subplot(2,3,3);
hold(ax3,'on'); grid(ax3,'on');
semilogy(ax3, t_min, drift_cond + 1e-12, 'Color', c_z, 'LineWidth', 1.4);
for tm = T_marks
    xline(ax3, tm, ':', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8);
end
xlabel(ax3, 'Time [min]');
ylabel(ax3, '$|\dot{y} + 2nx|$ [m/s]', 'Interpreter', 'latex');
title(ax3, 'Drift-free Condition');

%% Cumulative delta-v
ax4 = subplot(2,3,4);
hold(ax4,'on'); grid(ax4,'on');
plot(ax4, t_min, dv, 'Color', c_traj, 'LineWidth', 1.4);
for tm = T_marks
    xline(ax4, tm, ':', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8);
end
xlabel(ax4, 'Time [min]');
ylabel(ax4, '$\Delta v$ [m/s]', 'Interpreter', 'latex');
title(ax4, 'Cumulative $\Delta v$', 'Interpreter', 'latex');

%% Control inputs
ax5 = subplot(2,3,5);
hold(ax5,'on'); grid(ax5,'on');
plot(ax5, t_min, U_hist(1,:)*1000, 'Color', c_x, 'LineWidth', 1.2, 'DisplayName', '$u_x$');
plot(ax5, t_min, U_hist(2,:)*1000, 'Color', c_y, 'LineWidth', 1.2, 'DisplayName', '$u_y$');
plot(ax5, t_min, U_hist(3,:)*1000, 'Color', c_z, 'LineWidth', 1.2, 'DisplayName', '$u_z$');
yline(ax5,  con.u_max*1000, '--k', 'LineWidth', 0.8);
yline(ax5, -con.u_max*1000, '--k', 'LineWidth', 0.8);
for tm = T_marks
    xline(ax5, tm, ':', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8);
end
xlabel(ax5, 'Time [min]');
ylabel(ax5, 'Accel. [mm/s$^2$]', 'Interpreter', 'latex');
title(ax5, 'Control Inputs');
legend(ax5, 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 7);

%% States vs NMC reference
ax6 = subplot(2,3,6);
hold(ax6,'on'); grid(ax6,'on');
plot(ax6, t_min, X_hist(1,:),          'Color', c_x, 'LineWidth', 1.2, 'DisplayName', 'x');
plot(ax6, t_min, X_hist(2,:),          'Color', c_y, 'LineWidth', 1.2, 'DisplayName', 'y');
plot(ax6, t_min, X_nmc(1,1:n_steps),  '--', 'Color', c_x, 'LineWidth', 0.8, 'DisplayName', 'x_{nmc}');
plot(ax6, t_min, X_nmc(2,1:n_steps),  '--', 'Color', c_y, 'LineWidth', 0.8, 'DisplayName', 'y_{nmc}');
for tm = T_marks
    xline(ax6, tm, ':', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8);
end
xlabel(ax6, 'Time [min]'); ylabel(ax6, 'Position [m]');
title(ax6, 'States vs NMC Reference');
legend(ax6, 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 7);


%%  EXPORT JSON

metadata = struct(...
    'a_km',          a / 1000, ...
    'period_min',    T / 60, ...
    'dt',            dt, ...
    'n_steps',       n_steps, ...
    'controller',    'EMPC_NMC', ...
    'dynamics_model','CWH', ...
    'horizon_N',     N, ...
    'u_max',         con.u_max, ...
    'rho_m',         rho, ...
    'phi_0_rad',     phi_0, ...
    'y_c_m',         y_c, ...
    'dv_total_ms',   dv(end));

data = struct(...
    'metadata', metadata, ...
    't',        t, ...
    'X',        X_hist', ...
    'U',        U_hist', ...
    'X_nmc',    X_nmc(:,1:n_steps)');

json_str = jsonencode(data, 'PrettyPrint', true);
fid = fopen('exports/scenarios/sim_empc_nmc.json', 'w');
fprintf(fid, '%s', json_str);
fclose(fid);

fprintf('JSON exported.\n');
