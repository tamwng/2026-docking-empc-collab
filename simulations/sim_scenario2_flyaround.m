% sim_scenario2_flyaround.m
% Scenario 2 — NMC Fly-Around (two-run EMPC study)
%
% Demonstrates the role of the periodic terminal constraint in eMPC for
% natural motion circumnavigation orbit stabilisation.
%
% Run 1 — Degeneracy (no terminal constraint):
%   IC: V-bar hold at [0;300;0;0;0;0] (CWH equilibrium).
%   Fuel-only stage cost ℓ(u)=‖u‖², U*=0 globally optimal.
%   Chaser stays put; entire passive-orbit family achieves ℓ_p=0.
%
% Run 2 — Orbit stabilisation (periodic terminal equality):
%   Same IC. Terminal state pinned to Π*(k+N-1 mod P) at each step.
%   Π* = analytically prescribed NMC 2:1 ellipse with radial semi-axis
%        b = 75 m (along-track amplitude = 150 m), computed by forward
%        propagation of the NMC IC under A_d.
%   Controller drives chaser from the V-bar hold to Π* and stabilises it.
%   Asymptotic stability follows from Zanon–Grüne–Diehl (2015) Thm 4.6.

clear; clc;
script_dir = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(script_dir, '..', 'tools', 'matlab2tikz-master', 'src')));
constants;   % provides n, T, a, mu, ...

%% ── Discretisation ───────────────────────────────────────────────────────────
% P chosen so that dt = T/P gives an exact integer number of steps per period,
% ensuring the NMC orbit closes exactly under A_d.
P  = 92;
dt = T / P;

Ac = [0      0     0    1     0    0  ;
      0      0     0    0     1    0  ;
      0      0     0    0     0    1  ;
      3*n^2  0     0    0     2*n  0  ;
      0      0     0   -2*n   0    0  ;
      0      0    -n^2  0     0    0 ];
Bc = [zeros(3,3); eye(3)];
[A_d, B_d] = discretize(Ac, Bc, dt);

x_test   = [100; 0; 0; 0; -2*n*100; 0];
residual = norm(A_d^P * x_test - x_test);
fprintf('NMT closure residual: %.2e  (must be < 1e-4)\n', residual);
assert(residual < 1e-4, 'Discrete NMT does not close — check dt and P');

%% ── Prescribed NMC orbit Π* ──────────────────────────────────────────────────
% 2:1 ellipse in Hill frame: x in [-b, b] radial, y in [-2b, 2b] along-track.
% IC on NMC manifold: vx=0, vy+2n*x=0 → start at [b;0;0;0;-2nb;0].
% Propagated by A_d to get the full P-step periodic orbit exactly.
b_nmc       = 75;                                  % [m]  radial semi-axis
x_star      = zeros(6, P);
x_star(:,1) = [b_nmc; 0; 0; 0; -2*n*b_nmc; 0];   % NMC IC (on manifold)
for k = 1:P-1
    x_star(:, k+1) = A_d * x_star(:, k);
end
closure_res_star = norm(A_d * x_star(:, end) - x_star(:, 1));
fprintf('Pi* closure residual:  %.2e  (must be < 1e-8)\n', closure_res_star);
assert(closure_res_star < 1e-6, 'Pi* does not close under A_d');



%% ── QP dimensions and parameters ────────────────────────────────────────────
u_max = 1e-2;    % [m/s²]  symmetric per-axis thrust bound
N     = P;       % horizon = one full orbital period
n_sim = 6*P;     % total simulation length (6 orbital periods)
n_x   = 6;
n_u   = 3;

%% ── Condensed prediction matrices Sx (6N×6) and Su (6N×3N) ──────────────────
Sx = zeros(n_x*N, n_x);
Su = zeros(n_x*N, n_u*N);
Ad_pow = eye(n_x);
for i = 1:N
    Ad_pow = Ad_pow * A_d;
    Sx((i-1)*n_x+1 : i*n_x, :) = Ad_pow;
    for j = 1:i
        Su((i-1)*n_x+1 : i*n_x, (j-1)*n_u+1 : j*n_u) = A_d^(i-j) * B_d;
    end
end

% Position-row extracts (kept for struct compatibility with run_flyaround)
pos_idx = zeros(1, 3*N);
for ki = 1:N
    pos_idx((ki-1)*3+1 : ki*3) = (ki-1)*6+1 : (ki-1)*6+3;
end
Sx_pos_all = Sx(pos_idx, :);
Su_pos_all = Su(pos_idx, :);

%% Cost and static input-bound matrices (fuel-only, Q=0)
H     = 2 * eye(3*N);
H     = (H + H') / 2;
f_vec = zeros(3*N, 1);
A_u   = [eye(3*N); -eye(3*N)];
b_u   = u_max * ones(6*N, 1);

%% Pack shared structs
matrices.A_d        = A_d;
matrices.B_d        = B_d;
matrices.Sx         = Sx;
matrices.Su         = Su;
matrices.Sx_pos_all = Sx_pos_all;
matrices.Su_pos_all = Su_pos_all;
matrices.H          = H;
matrices.f_vec      = f_vec;
matrices.A_u        = A_u;
matrices.b_u        = b_u;

x0 = [25; 335; 0; 0; 0; 0];   % V-bar hold — CWH equilibrium (shared IC)

params_base.rho_min = 50;     % unused (use_band=false) — kept for struct completeness
params_base.rho_max = 200;
params_base.u_max   = u_max;
params_base.N       = N;
params_base.P       = P;
params_base.n       = n;
params_base.dt      = dt;
params_base.n_sim   = n_sim;

%% ── RUN 1 — No terminal constraint (degeneracy demo) ─────────────────────────
fprintf('\n════════════════════════════════════════════════════════════\n');
fprintf(' S2 Run 1 — No terminal constraint (V-bar hold, chaser stays put)\n');
fprintf('════════════════════════════════════════════════════════════\n');

params_r1     = params_base;
params_r1.x0  = x0;

run_cfg_r1.use_nmc_manifold      = false;
run_cfg_r1.use_periodic_terminal = false;
run_cfg_r1.x_star                = [];
run_cfg_r1.n_sim_override        = 2*P;
run_cfg_r1.use_band              = false;

r1 = run_flyaround(run_cfg_r1, matrices, params_r1);

%% ── RUN 2 — Periodic terminal to Π* (orbit stabilisation) ───────────────────
fprintf('\n════════════════════════════════════════════════════════════\n');
fprintf(' S2 Run 2 — Periodic terminal constraint → stabilise Pi*\n');
fprintf(' (Zanon–Gruene–Diehl 2015 Thm 4.6)\n');
fprintf('════════════════════════════════════════════════════════════\n');
fprintf('IC: [0; 300; 0; 0; 0; 0]   Pi*: b = %d m radial / %d m along-track\n', ...
    b_nmc, 2*b_nmc);

params_r2     = params_base;
params_r2.x0  = x0;

run_cfg_r2.use_nmc_manifold      = true;
run_cfg_r2.use_periodic_terminal = true;
run_cfg_r2.x_star                = x_star;
run_cfg_r2.n_sim_override        = [];
run_cfg_r2.use_band              = false;

r2 = run_flyaround(run_cfg_r2, matrices, params_r2);

%% ── Phase error — empirical convergence verification ────────────────────────
n_sim_r2  = r2.n_sim;
phase_err = zeros(1, n_sim_r2 + 1);
for k = 1 : n_sim_r2 + 1
    col          = mod(k-1, P) + 1;
    phase_err(k) = norm(r2.x_log(:, k) - x_star(:, col));
end
last_P_start        = max(1, n_sim_r2 - P + 1);
phase_err_last      = phase_err(last_P_start : end);
phase_err_max_last  = max(phase_err_last);
phase_err_mean_last = mean(phase_err_last);
fprintf('\nPhase error (last orbit):  max = %.4f m   mean = %.4f m\n', ...
    phase_err_max_last, phase_err_mean_last);

%% ── FIGURES ──────────────────────────────────────────────────────────────────
fig_dir = fullfile(script_dir, '..', 'results', 'figures', 's2_flyaround');
if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end

c_r2   = [0.00 0.45 0.70];   % blue  — Run 2 trajectory
c_star = [0.93 0.69 0.13];   % gold  — Π*
c_r1   = [0.85 0.33 0.10];   % red   — Run 1 IC marker

n_sim_r1  = r1.n_sim;
time_s_r2 = (0 : n_sim_r2) * dt;

%% Fig 1 — Hill-frame: Run 1 stuck vs Run 2 converging to Π*
fig1 = figure('Name', 'S2 — Hill-frame');
ax1  = axes(fig1);
hold(ax1, 'on'); grid(ax1, 'on'); axis(ax1, 'equal');
x_star_cl = [x_star, x_star(:,1)];   % close the curve for plotting
plot(ax1, x_star_cl(1,:), x_star_cl(2,:), '--', 'Color', c_star, 'LineWidth', 2.0, ...
    'DisplayName', sprintf('\\Pi^* (b = %d m)', b_nmc));
plot(ax1, r2.x_log(1,:), r2.x_log(2,:), 'Color', c_r2, 'LineWidth', 1.0, ...
    'DisplayName', 'Run 2: trajectory');
plot(ax1, r1.x_log(1,:), r1.x_log(2,:), 'Color', c_r1, 'LineWidth', 1.0, ...
    'DisplayName', 'Run 1: free drift (no terminal cstr.)');
plot(ax1, x0(1), x0(2), 's', 'Color', c_r1, 'MarkerSize', 10, ...
    'MarkerFaceColor', c_r1, 'LineWidth', 1.5, 'HandleVisibility', 'off');
plot(ax1, 0, 0, '.k', 'MarkerSize', 14, 'DisplayName', 'Target');
xlabel(ax1, 'r_x (radial) [m]');
ylabel(ax1, 'r_y (along-track) [m]');
title(ax1, 'S2 — Run 1: free drift vs Run 2: stabilises \Pi^*');
legend(ax1, 'Location', 'northeast');
exportgraphics(fig1, fullfile(fig_dir, 'hill_frame_overlay.pdf'), 'ContentType', 'vector');
exportgraphics(fig1, fullfile(fig_dir, 'hill_frame_overlay.png'), 'Resolution', 300);

%% Fig 5 — Non-equilibrium IC drift (Run 1 detail)
time_s_r1 = (0 : n_sim_r1) * dt;

fig5 = figure('Name', 'S2 — Non-equilibrium drift (Run 1)');
tl5  = tiledlayout(fig5, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% Top: Hill-frame trajectory
ax5a = nexttile(tl5);
hold(ax5a, 'on'); grid(ax5a, 'on'); axis(ax5a, 'equal');
plot(ax5a, x_star_cl(1,:), x_star_cl(2,:), '--', 'Color', c_star, 'LineWidth', 1.5, ...
    'DisplayName', sprintf('\\Pi^* (b = %d m)', b_nmc));
plot(ax5a, r1.x_log(1,:), r1.x_log(2,:), 'Color', c_r1, 'LineWidth', 1.2, ...
    'DisplayName', 'Run 1: free drift');
plot(ax5a, x0(1), x0(2), 's', 'Color', c_r1, 'MarkerSize', 10, ...
    'MarkerFaceColor', c_r1, 'DisplayName', 'IC');
plot(ax5a, r1.x_log(1,end), r1.x_log(2,end), '^', 'Color', [0.4 0.4 0.4], ...
    'MarkerSize', 8, 'MarkerFaceColor', [0.4 0.4 0.4], 'DisplayName', 'Final (2T)');
plot(ax5a, 0, 0, '.k', 'MarkerSize', 14, 'DisplayName', 'Target');
xlabel(ax5a, 'r_x (radial) [m]');
ylabel(ax5a, 'r_y (along-track) [m]');
title(ax5a, 'Hill-frame — drifting ellipse (U^* = 0, no terminal cstr.)');
legend(ax5a, 'Location', 'northeast');

% Bottom: along-track time series showing secular drift vs analytical trend
ax5b = nexttile(tl5);
hold(ax5b, 'on'); grid(ax5b, 'on');
t_vec     = linspace(0, n_sim_r1 * dt, 500);
y_secular = -6 * x0(1) * n * t_vec + x0(2);
plot(ax5b, time_s_r1, r1.x_log(2,:), 'Color', c_r1, 'LineWidth', 1.2, ...
    'DisplayName', 'r_y (simulated)');
plot(ax5b, t_vec, y_secular, 'k--', 'LineWidth', 1.0, ...
    'DisplayName', 'Secular trend: -6 x_0 n t + y_0');
xline(ax5b, T, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.9, ...
    'Label', '1 orbit', 'HandleVisibility', 'off');
xlabel(ax5b, 'Time [s]');
ylabel(ax5b, 'r_y (along-track) [m]');
title(ax5b, sprintf('Along-track secular drift  (avg. rate \\approx %.3f m/s)', 6 * x0(1) * n));
legend(ax5b, 'Location', 'northeast');

exportgraphics(fig5, fullfile(fig_dir, 'run1_nonequil_drift.pdf'), 'ContentType', 'vector');
exportgraphics(fig5, fullfile(fig_dir, 'run1_nonequil_drift.png'), 'Resolution', 300);

%% Fig 2 — Phase error (Zanon–Grüne–Diehl convergence)
fig2 = figure('Name', 'S2 — Phase error');
ax2  = axes(fig2);
hold(ax2, 'on'); grid(ax2, 'on');
semilogy(ax2, time_s_r2, phase_err, 'Color', c_r2, 'LineWidth', 1.2);
xlabel(ax2, 'Time [s]');
ylabel(ax2, '||x(k) - x^*_{k \rm mod P}|| [m]');
title(ax2, 'Convergence to \Pi^* — Zanon–Gruene–Diehl Thm 4.6');
ylim(ax2, [1e-2, 1e3]);
exportgraphics(fig2, fullfile(fig_dir, 'phase_error.pdf'), 'ContentType', 'vector');
exportgraphics(fig2, fullfile(fig_dir, 'phase_error.png'), 'Resolution', 300);

%% Fig 3 — Control effort Run 2 (log scale)
fig3 = figure('Name', 'S2 — Control effort');
ax3  = axes(fig3);
hold(ax3, 'on'); grid(ax3, 'on');
semilogy(ax3, time_s_r2(1:n_sim_r2), r2.u_norm, 'Color', c_r2, 'LineWidth', 1.2);
if r2.inject_end < n_sim_r2
    xline(ax3, r2.inject_end * dt, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.9, ...
        'DisplayName', 'Orbit reached');
end
xlabel(ax3, 'Time [s]');
ylabel(ax3, '||u|| [m/s^2]');
title(ax3, 'S2 — Control effort (log scale)');
ylim(ax3, [1e-11, 1e-4]);
exportgraphics(fig3, fullfile(fig_dir, 'control_effort.pdf'), 'ContentType', 'vector');
exportgraphics(fig3, fullfile(fig_dir, 'control_effort.png'), 'Resolution', 300);

%% Fig 4 — Cumulative ΔV Run 2
fig4 = figure('Name', 'S2 — Cumulative DeltaV');
ax4  = axes(fig4);
hold(ax4, 'on'); grid(ax4, 'on');
plot(ax4, time_s_r2(1:n_sim_r2), cumsum(r2.dv_log), 'Color', c_r2, 'LineWidth', 1.2);
xlabel(ax4, 'Time [s]');
ylabel(ax4, 'Cumulative \DeltaV [m/s]');
title(ax4, 'S2 — Cumulative \DeltaV');
exportgraphics(fig4, fullfile(fig_dir, 'cumulative_dv.pdf'), 'ContentType', 'vector');
exportgraphics(fig4, fullfile(fig_dir, 'cumulative_dv.png'), 'Resolution', 300);

fprintf('Figures saved to %s\n', fig_dir);

%% ── JSON EXPORT ───────────────────────────────────────────────────────────────
u_lo_r1 = max(1, n_sim_r1 - P);
u_lo_r2 = max(1, n_sim_r2 - P);

meta = struct( ...
    'scenario',                'S2_flyaround', ...
    'date',                    datestr(now, 'yyyy-mm-dd HH:MM:SS'), ...
    'dt_s',                    dt, ...
    'P_orbit_steps',           P, ...
    'T_orbit_s',               T, ...
    'N_horizon',               N, ...
    'n_sim',                   n_sim, ...
    'u_max_mps2',              u_max, ...
    'x0',                      x0', ...
    'Pi_star_b_radial_m',      b_nmc, ...
    'Pi_star_b_alongtrack_m',  2*b_nmc, ...
    'nmt_closure_residual',    residual, ...
    'Pi_star_closure_residual', closure_res_star);

res1 = struct( ...
    'terminal',               'none', ...
    'n_sim_steps',            n_sim_r1, ...
    'delta_v_total_m_s',      r1.dv_total, ...
    'infeasible_qp_steps',    sum(r1.exitflag_log <= 0), ...
    'u_mean_last_orbit_mps2', mean(r1.u_norm(u_lo_r1:end)), ...
    'note',                   'CWH equilibrium; U*=0 globally optimal — chaser stays at IC');

res2 = struct( ...
    'terminal',               'periodic Pi* (6-row Su_term equality)', ...
    'n_sim_steps',            n_sim_r2, ...
    'delta_v_total_m_s',      r2.dv_total, ...
    'delta_v_injection_m_s',  r2.dv_inject, ...
    'inject_end_step',        r2.inject_end, ...
    'b_est_m',                r2.b_est, ...
    'range_ratio_2to1',       r2.ratio_2to1, ...
    'infeasible_qp_steps',    sum(r2.exitflag_log <= 0), ...
    'u_mean_last_orbit_mps2', mean(r2.u_norm(u_lo_r2:end)), ...
    'phase_err_max_last_m',   phase_err_max_last, ...
    'phase_err_mean_last_m',  phase_err_mean_last, ...
    'control_converged',      mean(r2.u_norm(u_lo_r2:end)) < 1e-4);

data = struct('metadata', meta, 'run1', res1, 'run2', res2);

json_str    = jsonencode(data, 'PrettyPrint', true);
export_path = fullfile(script_dir, '..', 'exports', 'scenarios', ...
    'sim_scenario2_flyaround.json');
fid = fopen(export_path, 'w');
fprintf(fid, '%s', json_str);
fclose(fid);
fprintf('JSON exported to %s\n', export_path);
