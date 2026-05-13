% sim_mpc_regulation.m
% Simulates Regulation MPC rendezvous in the Hill frame using CWH equations.
%
% drives the relative state to zero (origin = target) subject to hard input constraints.
clear; clc;


%%  LOAD CONSTANTS
constants;

%%  SIMULATION PARAMETERS (as before)
dt      = 60;
t_final = 2 * T;
n_steps = round(t_final / dt);
t       = (0:n_steps-1) * dt;


%%  MPC PARAMETERS

N     = 20;       % [-]      prediction horizon (initial guess)
u_max = 1e-2;     % [m/s^2]  thrust limit 


%%  CWH SS

Ac = [0      0     0    1     0    0  ;
      0      0     0    0     1    0  ;
      0      0     0    0     0    1  ;
      3*n^2  0     0    0     2*n  0  ;
      0      0     0   -2*n   0    0  ;
      0      0    -n^2  0     0    0 ];

Bc = [zeros(3,3); eye(3)];


%%  DISCRETIZE
[Ad, Bd] = discretize(Ac, Bc, dt);


%%  COST MATRICES (same as LQR for direct comparison)
Q = diag([1e-2, 1e-2, 1e-2, 1e0, 1e0, 1e0]);
R = diag([1e4,  1e4,  1e4]);


%%  TERMINAL COST — P from LQR 
[~, P] = lqr_controller(Ad, Bd, Q, R);


%%  INITIAL CONDITIONS
X0 = [500; 6000; 200; 0; 0; 0];


%%  INTEGRATE TARGET ORBIT (as before)
X_t0     = [a; 0; 0; 0; sqrt(mu/a); 0];
f_target = @(X) two_body(X, mu);
X_t      = rk4_integrator(f_target, X_t0, dt, n_steps);


%%  CLOSED-LOOP MPC SIMULATION
X_hist = zeros(6, n_steps);
U_hist = zeros(3, n_steps);
X_hist(:, 1) = X0;

fprintf('Running Regulation MPC (%d steps, N=%d)...\n', n_steps, N);
for k = 1:n_steps-1
    [u_opt, ~, ~]  = mpc_regulation(X_hist(:,k), Ad, Bd, Q, R, P, N, u_max);
    U_hist(:, k)   = u_opt;
    X_hist(:, k+1) = Ad * X_hist(:,k) + Bd * u_opt;
end
U_hist(:,end) = mpc_regulation(X_hist(:,end), Ad, Bd, Q, R, P, N, u_max);

fprintf('Simulation complete.\n');


%%  POST-PROCESSING

t_min = t / 60;

% Relative distance at each timestep
rel_dist = vecnorm(X_hist(1:3,:), 2, 1);   % [1 x n_steps]

% Lyapunov function V = x'*P*x at each timestep
V_hist = zeros(1, n_steps);
for k = 1:n_steps
    V_hist(k) = X_hist(:,k)' * P * X_hist(:,k);
end

% Convergence criterion
tol      = 0.01;                            % [m] convergence threshold
conv_idx = find(rel_dist < tol, 1, 'first');
if ~isempty(conv_idx)
    fprintf('Converged at t = %.1f min (step %d)\n', t_min(conv_idx), conv_idx);
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
plot(ax1, t_min, X_hist(1,:), 'Color', c_x, 'LineWidth', 1.2, 'DisplayName', 'x (radial)');
plot(ax1, t_min, X_hist(2,:), 'Color', c_y, 'LineWidth', 1.2, 'DisplayName', 'y (along-track)');
plot(ax1, t_min, X_hist(3,:), 'Color', c_z, 'LineWidth', 1.2, 'DisplayName', 'z (cross-track)');
xlabel(ax1, 'Time [min]'); ylabel(ax1, 'Position [m]');
title(ax1, 'Position States'); legend(ax1, 'Location', 'northeast');

ax2 = subplot(2, 3, 2);
hold(ax2, 'on'); grid(ax2, 'on');
plot(ax2, t_min, X_hist(4,:), 'Color', c_x, 'LineWidth', 1.2, 'DisplayName', '$\dot{x}$');
plot(ax2, t_min, X_hist(5,:), 'Color', c_y, 'LineWidth', 1.2, 'DisplayName', '$\dot{y}$');
plot(ax2, t_min, X_hist(6,:), 'Color', c_z, 'LineWidth', 1.2, 'DisplayName', '$\dot{z}$');
xlabel(ax2, 'Time [min]'); ylabel(ax2, 'Velocity [m/s]');
title(ax2, 'Velocity States');
legend(ax2, 'Location', 'northeast', 'Interpreter', 'latex');

ax3 = subplot(2, 3, 3);
hold(ax3, 'on'); grid(ax3, 'on');
plot(ax3, t_min, rel_dist, 'Color', c_x, 'LineWidth', 1.2);
xlabel(ax3, 'Time [min]'); ylabel(ax3, 'Distance [m]');
title(ax3, 'Relative Distance');

ax4 = subplot(2, 3, 4);
hold(ax4, 'on'); grid(ax4, 'on');
plot(ax4, t_min, U_hist(1,:)*1000, 'Color', c_x, 'LineWidth', 1.2, 'DisplayName', '$u_x$');
plot(ax4, t_min, U_hist(2,:)*1000, 'Color', c_y, 'LineWidth', 1.2, 'DisplayName', '$u_y$');
plot(ax4, t_min, U_hist(3,:)*1000, 'Color', c_z, 'LineWidth', 1.2, 'DisplayName', '$u_z$');
xlabel(ax4, 'Time [min]');
ylabel(ax4, 'Acceleration [mm/s$^2$]', 'Interpreter', 'latex');
title(ax4, 'Control Inputs');
legend(ax4, 'Location', 'northeast', 'Interpreter', 'latex');
yline(ax4,  u_max*1000, '--k', 'LineWidth', 0.8, 'DisplayName', 'u\_max');
yline(ax4, -u_max*1000, '--k', 'LineWidth', 0.8);

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
dv = cumsum(vecnorm(U_hist, 2, 1) * dt);
plot(ax6, t_min, dv, 'Color', c_z, 'LineWidth', 1.2);
xlabel(ax6, 'Time [min]');
ylabel(ax6, '$\Delta v$ [m/s]', 'Interpreter', 'latex');
title(ax6, 'Cumulative $\Delta v$', 'Interpreter', 'latex');

fprintf('Total delta-v: %.4f m/s\n', dv(end));

%% -------------------------------------------------------------------------
%  EXPORT JSON
%  -------------------------------------------------------------------------
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
    'u_max',          u_max);

data = struct(...
    'metadata', metadata, ...
    't',        t, ...
    'X_t',      X_t', ...
    'Rho',      X_hist');

json_str = jsonencode(data, 'PrettyPrint', true);
fid = fopen('exports/scenarios/sim_mpc_regulation.json', 'w');
fprintf(fid, '%s', json_str);
fclose(fid);

fprintf('JSON exported.\n');