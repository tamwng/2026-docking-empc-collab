% sim_scenario2_flyaround.m
% Scenario 2 — NMC Fly-Around: Economic MPC, Run 1
%
% One-step receding-horizon EMPC drives the chaser from the V-bar hold point
% (handover from S4, ~300 m) into a fuel-optimal natural motion circumnavigation
% (NMC) orbit inside the keep-in band [rho_min, rho_max]. The controller is given
% only the band — no reference orbit. It discovers the NMC by minimising fuel.
%
% Stage cost:  l(u) = ||u||^2  (no state penalty)
% No terminal cost or terminal constraint.
% Band constraints are linearised per prediction stage around the warm-start.

clear; clc;

%% Paths
script_dir = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(script_dir, '..', 'tools', 'matlab2tikz-master', 'src')));

%% Constants (mu, a, T, n, ...)
constants;

%% Discretisation at P-exact time step (NMT must close in integer steps)
P  = 92;        % integer steps per orbital period
dt = T / P;     % ≈ 60.27 s — do NOT hardcode 60

Ac = [0      0     0    1     0    0  ;
      0      0     0    0     1    0  ;
      0      0     0    0     0    1  ;
      3*n^2  0     0    0     2*n  0  ;
      0      0     0   -2*n   0    0  ;
      0      0    -n^2  0     0    0 ];
Bc = [zeros(3,3); eye(3)];
[A_d, B_d] = discretize(Ac, Bc, dt);

%% NMT closure sanity check — abort if discrete NMT does not close in P steps
x_nmt_test = [100; 0; 0; 0; -2*n*100; 0];   % NMC condition: vy0 = -2*n*x0
residual   = norm(A_d^P * x_nmt_test - x_nmt_test);
fprintf('NMT closure residual: %.2e  (must be < 1e-4)\n', residual);
assert(residual < 1e-4, 'Discrete NMT does not close — check dt and P');

%% Initial condition — NMC injection at radial axis
% Chaser arrives at rho_b = 150 m on the radial axis with NMC injection
% velocity vy = -2*n*rho_b. Free drift traces a 2:1 ellipse with range
% 150-300 m (exceeds outer band), so EMPC must contract onto a smaller NMC.
rho_b_ic = 150;   % [m]  desired injection semi-axis
x0 = [rho_b_ic; 0; 0; 0; -2*n*rho_b_ic; 0];
fprintf('IC: NMC injection  rho_b=%d m  [%.2f; %.2f; %.4f; %.5f; %.5f; %.2e]\n', rho_b_ic, x0);

%% Parameters
rho_min = 50;       % [m]      inner keep-out radius
rho_max = 200;      % [m]      outer keep-in radius
u_max   = 1e-2;     % [m/s^2]  per-axis thrust bound
N       = P;        % horizon  = 1 full orbital period
n_sim   = 6 * P;    % 552 steps ≈ 6 orbital periods

n_x = 6;  n_u = 3;

%% Build condensed prediction matrices Sx (6N×6) and Su (6N×3N)
Sx = zeros(n_x * N, n_x);
Su = zeros(n_x * N, n_u * N);

Ad_pow = eye(n_x);
for i = 1:N
    Ad_pow = Ad_pow * A_d;
    Sx((i-1)*n_x+1 : i*n_x, :) = Ad_pow;
    for j = 1:i
        Su((i-1)*n_x+1 : i*n_x, (j-1)*n_u+1 : j*n_u) = A_d^(i-j) * B_d;
    end
end

%% Pre-extract position rows of Sx and Su to avoid repeated indexing in loop
% Position rows at stage ki occupy rows 6(ki-1)+1 : 6(ki-1)+3
pos_row_idx = zeros(1, 3*N);
for ki = 1:N
    pos_row_idx((ki-1)*3+1 : ki*3) = (ki-1)*6+1 : (ki-1)*6+3;
end
Sx_pos_all = Sx(pos_row_idx, :);   % 3N × 6
Su_pos_all = Su(pos_row_idx, :);   % 3N × 3N

%% NMC terminal manifold equality constraint  (computed once — constant in U)
% Forces the planned trajectory to end on the NMC manifold:
%   vx(N) = 0
%   vy(N) + 2*n*x(N) = 0
% State order [x,y,z,vx,vy,vz] → terminal block rows 1(x), 4(vx), 5(vy).
Sx_term = Sx((N-1)*n_x+1 : N*n_x, :);   % 6 × 6
Su_term = Su((N-1)*n_x+1 : N*n_x, :);   % 6 × 3N

Aeq          = [Su_term(4, :);
                Su_term(5, :) + 2*n * Su_term(1, :)];   % 2 × 3N
Sx_term_nmc  = [Sx_term(4, :);
                Sx_term(5, :) + 2*n * Sx_term(1, :)];   % 2 × 6  (for beq in loop)

%% Cost matrices (Q = 0, fuel-only stage cost)
% H = 2 * kron(eye(N), R_cost),  R_cost = eye(3)
H = 2 * eye(3*N);          % 3N × 3N
H = (H + H') / 2;          % symmetrise for quadprog
f = zeros(3*N, 1);         % constant (no state cost)

%% Static input-bound constraint (box constraint on u)
A_u = [eye(3*N); -eye(3*N)];       % 6N × 3N
b_u = u_max * ones(6*N, 1);

%% Pre-allocate logs
x_log        = zeros(6, n_sim + 1);
u_log        = zeros(3, n_sim);
exitflag_log = zeros(1, n_sim);

x_log(:, 1) = x0;
U_warm       = zeros(3*N, 1);   % warm-start input sequence

%% Simulation loop
opts = optimoptions('quadprog', 'Display', 'off');

fprintf('Running S2 simulation (%d steps, %.1f orbits)...\n', n_sim, n_sim / P);
t_start = tic;

for k = 1:n_sim

    %% 1. Nominal positions from warm-start prediction
    x_pred_warm = Sx * x_log(:, k) + Su * U_warm;   % 6N × 1
    r_nom = reshape(x_pred_warm(1:6*N), 6, N);       % 6 × N
    r_nom = r_nom(1:3, :);                           % 3 × N  (positions only)

    %% 2. Build linearised band constraints
    r_norms   = vecnorm(r_nom, 2, 1);   % 1 × N
    A_band    = zeros(N, 3*N);
    b_outer   = zeros(N, 1);
    b_inner   = zeros(N, 1);

    for ki = 1:N
        idx3 = (ki-1)*3+1 : ki*3;

        if r_norms(ki) < 1e-6
            nk = [0; 1; 0];              % fallback: along-track direction
        else
            nk = r_nom(:, ki) / r_norms(ki);
        end

        row   = nk' * Su_pos_all(idx3, :);                       % 1 × 3N
        b_sx  = nk' * Sx_pos_all(idx3, :) * x_log(:, k);        % scalar

        A_band(ki, :) = row;
        b_outer(ki)   = rho_max - b_sx;   % outer:  n'*r ≤ rho_max
        b_inner(ki)   = rho_min - b_sx;   % inner:  n'*r ≥ rho_min  → -n'*r ≤ -rho_min
    end

    %% 3. Assemble combined inequality system
    % [A_u; A_band; -A_band] * U ≤ [b_u; b_outer; -b_inner]
    A_ineq = [A_u;  A_band; -A_band];
    b_ineq = [b_u; b_outer; -b_inner];

    %% 4. NMC terminal manifold RHS (state-dependent)
    beq = -Sx_term_nmc * x_log(:, k);   % 2 × 1

    %% 5. Solve QP
    [U_opt, ~, exitflag] = quadprog(H, f, A_ineq, b_ineq, Aeq, beq, [], [], U_warm, opts);
    exitflag_log(k) = exitflag;

    %% 6. Extract and apply first input; shift warm-start
    if exitflag > 0
        u_log(:, k) = U_opt(1:3);
        U_warm = [U_opt(4:end); zeros(3, 1)];
    else
        u_log(:, k) = zeros(3, 1);
        U_warm = [U_warm(4:end); zeros(3, 1)];
    end

    %% 7. Simulate one step
    x_log(:, k+1) = A_d * x_log(:, k) + B_d * u_log(:, k);

end

fprintf('Simulation complete (%.1f s wall-clock)\n', toc(t_start));

%% ── POST-PROCESS ──────────────────────────────────────────────────────────────
fprintf('\n--- S2 Run 1 Results ---\n');
fprintf('dt = %.4f s  (P=%d, T=%.2f s)\n', dt, P, T);
fprintf('NMT closure residual: %.2e\n', residual);

range_log = vecnorm(x_log(1:3, :));        % 1 × (n_sim+1)
dv_log    = vecnorm(u_log, 2, 1) * dt;     % per-step ΔV [m/s]

% Injection boundary: first step where rolling mean ‖u‖ over P/4 steps < 1e-4
u_norm       = vecnorm(u_log, 2, 1);       % 1 × n_sim
rolling_mean = movmean(u_norm, P/4);
inject_end   = find(rolling_mean < 1e-4, 1);
if isempty(inject_end); inject_end = n_sim; end

dv_total  = sum(dv_log);
dv_inject = sum(dv_log(1:inject_end));
dv_steady = sum(dv_log(inject_end+1:end));

range_ss    = range_log(end - 2*P : end);   % last 2 orbits
b_est       = min(range_ss);               % NMT radial semi-axis ≈ min range
ratio_2to1  = max(range_ss) / min(range_ss);

u_mean_last = mean(u_norm(end-P:end));
u_max_last  = max(u_norm(end-P:end));

fprintf('Total ΔV:       %.4f m/s\n', dv_total);
fprintf('  Injection ΔV: %.4f m/s  (steps 1–%d)\n', dv_inject, inject_end);
fprintf('  Steady-state: %.4f m/s  (steps %d–%d)\n', dv_steady, inject_end+1, n_sim);
fprintf('Discovered NMT semi-axis b ≈ %.1f m\n', b_est);
fprintf('Range ratio max/min (expect ≈ 2.0): %.3f\n', ratio_2to1);
fprintf('Band violations: %d\n', sum(range_log < rho_min | range_log > rho_max));
fprintf('Infeasible QP steps: %d\n', sum(exitflag_log <= 0));
fprintf('Mean ‖u‖ last orbit: %.2e m/s²\n', u_mean_last);
fprintf('Max  ‖u‖ last orbit: %.2e m/s²\n', u_max_last);

%% ── FIGURES ───────────────────────────────────────────────────────────────────
fig_dir = fullfile(script_dir, '..', 'results', 'figures', 's2_flyaround');
if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end

time_s  = (0:n_sim) * dt;
theta_c = linspace(0, 2*pi, 300);

%% Figure 1 — Hill-frame trajectory (rx vs ry, top view)
fig1 = figure('Name', 'S2 Run 1 — Hill-frame trajectory');
ax1  = axes(fig1);
hold(ax1, 'on'); grid(ax1, 'on'); axis(ax1, 'equal');

plot(ax1, rho_min * cos(theta_c), rho_min * sin(theta_c), '--r', 'LineWidth', 1.0, ...
    'DisplayName', sprintf('\\rho_{min} = %d m', rho_min));
plot(ax1, rho_max * cos(theta_c), rho_max * sin(theta_c), '--g', 'LineWidth', 1.0, ...
    'DisplayName', sprintf('\\rho_{max} = %d m', rho_max));
plot(ax1, x_log(1, :), x_log(2, :), 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8, ...
    'DisplayName', 'Trajectory');

idx_ss = (n_sim - 2*P) : n_sim;
plot(ax1, x_log(1, idx_ss+1), x_log(2, idx_ss+1), ...
    'Color', [0.00 0.45 0.70], 'LineWidth', 1.4, ...
    'DisplayName', 'Last 2 orbits (settled)');

plot(ax1, x0(1), x0(2), 'x', 'Color', [0.3 0.3 0.3], 'MarkerSize', 10, ...
    'LineWidth', 2.0, 'DisplayName', sprintf('IC (%d m)', rho_b_ic));
plot(ax1, 0, 0, '.k', 'MarkerSize', 14, 'DisplayName', 'Client');

xlabel(ax1, 'r_x (radial) [m]');
ylabel(ax1, 'r_y (along-track) [m]');
title(ax1, 'S2 Run 1 — Hill-frame trajectory');
legend(ax1, 'Location', 'northeast');

savefig(fig1, fullfile(fig_dir, 'hill_frame_trajectory.fig'));
exportgraphics(fig1, fullfile(fig_dir, 'hill_frame_trajectory.pdf'), 'ContentType', 'vector');

%% Figure 2 — Range history (‖r‖ vs time)
fig2 = figure('Name', 'S2 Run 1 — Range history');
ax2  = axes(fig2);
hold(ax2, 'on'); grid(ax2, 'on');

plot(ax2, time_s, range_log, 'Color', [0.00 0.45 0.70], 'LineWidth', 1.2);
yline(ax2, rho_min, '--r', 'LineWidth', 1.0, ...
    'DisplayName', sprintf('\\rho_{min} = %d m', rho_min));
yline(ax2, rho_max, '--g', 'LineWidth', 1.0, ...
    'DisplayName', sprintf('\\rho_{max} = %d m', rho_max));
xline(ax2, inject_end * dt, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.9, ...
    'DisplayName', 'Injection end');

xlabel(ax2, 'Time [s]');
ylabel(ax2, 'Range ||r|| [m]');
title(ax2, 'S2 Run 1 — Range history');
legend(ax2, 'Location', 'northeast');

savefig(fig2, fullfile(fig_dir, 'range_history.fig'));
exportgraphics(fig2, fullfile(fig_dir, 'range_history.pdf'), 'ContentType', 'vector');

%% Figure 3 — Control effort (log scale)
fig3 = figure('Name', 'S2 Run 1 — Control effort');
ax3  = axes(fig3);
hold(ax3, 'on'); grid(ax3, 'on');

semilogy(ax3, time_s(1:n_sim), u_norm, 'Color', [0.00 0.45 0.70], 'LineWidth', 1.2);
xline(ax3, inject_end * dt, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.9, ...
    'DisplayName', 'Injection end');

xlabel(ax3, 'Time [s]');
ylabel(ax3, '||u|| [m/s^2]');
title(ax3, 'S2 Run 1 — Control effort (log scale)');
legend(ax3, 'Location', 'northeast');

savefig(fig3, fullfile(fig_dir, 'control_effort.fig'));
exportgraphics(fig3, fullfile(fig_dir, 'control_effort.pdf'), 'ContentType', 'vector');

%% Figure 4 — Cumulative ΔV
fig4 = figure('Name', 'S2 Run 1 — Cumulative DeltaV');
ax4  = axes(fig4);
hold(ax4, 'on'); grid(ax4, 'on');

plot(ax4, time_s(1:n_sim), cumsum(dv_log), 'Color', [0.00 0.45 0.70], 'LineWidth', 1.2);

xlabel(ax4, 'Time [s]');
ylabel(ax4, 'Cumulative \DeltaV [m/s]');
title(ax4, 'S2 Run 1 — Cumulative \DeltaV');

savefig(fig4, fullfile(fig_dir, 'cumulative_dv.fig'));
exportgraphics(fig4, fullfile(fig_dir, 'cumulative_dv.pdf'), 'ContentType', 'vector');

fprintf('Figures saved to %s\n', fig_dir);

%% ── JSON EXPORT ───────────────────────────────────────────────────────────────
metadata_s = struct( ...
    'scenario',       'S2_flyaround', ...
    'run',            1, ...
    'description',    ['One-step EMPC, fuel-only stage cost, no terminal cost, ' ...
                       'outer+inner band constraints (linearised per stage). ' ...
                       'EMPC discovers NMT orbit without being told which orbit to track.'], ...
    'dt_s',           dt, ...
    'P_orbit_steps',  92, ...
    'T_orbit_s',      T, ...
    'N_horizon',      92, ...
    'n_sim_steps',    552, ...
    'u_max_mps2',     0.01, ...
    'rho_min_m',      50, ...
    'rho_max_m',      200, ...
    'x0_m',           x0', ...
    'Q_stage',        'zeros(6)', ...
    'R_cost',         'eye(3)', ...
    'terminal_cost',  'none');

results_s = struct( ...
    'nmt_closure_residual',        residual, ...
    'delta_v_total_m_s',           dv_total, ...
    'delta_v_injection_m_s',       dv_inject, ...
    'delta_v_injection_steps',     inject_end, ...
    'delta_v_steady_state_m_s',    dv_steady, ...
    'discovered_nmt_semiaxis_m',   b_est, ...
    'range_ratio_2to1',            ratio_2to1, ...
    'max_range_steady_m',          max(range_ss), ...
    'min_range_steady_m',          min(range_ss), ...
    'band_violations',             sum(range_log < rho_min | range_log > rho_max), ...
    'infeasible_qp_steps',         sum(exitflag_log <= 0), ...
    'u_mean_last_orbit_mps2',      u_mean_last, ...
    'u_max_last_orbit_mps2',       u_max_last, ...
    'control_converged',           u_mean_last < 1e-4);

data = struct('metadata', metadata_s, 'results', results_s);

json_str    = jsonencode(data, 'PrettyPrint', true);
export_path = fullfile(script_dir, '..', 'exports', 'scenarios', ...
                       'sim_scenario2_flyaround_run1.json');
fid = fopen(export_path, 'w');
fprintf(fid, '%s', json_str);
fclose(fid);

fprintf('JSON exported to %s\n', export_path);
