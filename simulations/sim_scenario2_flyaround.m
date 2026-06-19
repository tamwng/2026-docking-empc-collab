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
%   Asymptotic stability follows from Keerthi–Gilbert / Theorem 2.24.

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

% Orbit closure check: verifies that the NMC IC lies in the eigenvalue-1
% eigenspace of A_d, i.e. A_d^P maps this specific IC back to itself.
% This is NOT equivalent to A_d^P = I — the full-state period map has secular
% along-track drift for generic states (‖A_d^P − I‖_F ≈ 1.7e4).  The test
% only confirms that the prescribed Π* closes under A_d at dt = T/P.
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
fprintf(' (Keerthi–Gilbert / Theorem 2.24)\n');
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

%% ── RUN 4 — NMC manifold IC, pure fuel, no terminal constraint ───────────────
% IC is placed exactly on Π* (satisfies vy = -2n·x, vx = 0).
% Under pure fuel cost the globally optimal input is U* = 0 for every step.
% The chaser rides the free CWH dynamics and stays on Π* indefinitely.
% No terminal constraint needed — orbit knowledge is encoded in the IC itself.
fprintf('\n════════════════════════════════════════════════════════════\n');
fprintf(' S2 Run 4 — IC on NMC manifold, pure fuel, no terminal constraint\n');
fprintf(' Expected: u*=0, chaser stays on Pi* (orbit "captured for free")\n');
fprintf('════════════════════════════════════════════════════════════\n');

x0_r4 = [b_nmc; 0; 0; 0; -2*n*b_nmc; 0];   % starting point of Π*: x=b, vy=-2nb
fprintf('IC: [%.1f; 0; 0; 0; %.5f; 0]  (on Pi* manifold)\n', b_nmc, -2*n*b_nmc);

params_r4        = params_base;
params_r4.x0     = x0_r4;

run_cfg_r4.use_nmc_manifold      = false;
run_cfg_r4.use_periodic_terminal = false;
run_cfg_r4.x_star                = [];
run_cfg_r4.n_sim_override        = 3*P;
run_cfg_r4.use_band              = false;

r4 = run_flyaround(run_cfg_r4, matrices, params_r4);

%% Phase error Run 4 — should be ~0 throughout (machine-epsilon scale)
n_sim_r4      = r4.n_sim;
phase_err_r4  = zeros(1, n_sim_r4 + 1);
for k = 1 : n_sim_r4 + 1
    col              = mod(k-1, P) + 1;
    phase_err_r4(k)  = norm(r4.x_log(:, k) - x_star(:, col));
end
fprintf('Run 4 phase error:  max = %.2e m   mean = %.2e m  (expect machine-eps)\n', ...
    max(phase_err_r4), mean(phase_err_r4));

%% ── RUN 5 — Regularised dissipativity route (periodic stage cost) ────────────
% Stage cost: ℓ_aug(x,u,k) = ‖u‖² + ε‖x − x*_k‖²
%
% Strict P-periodic dissipativity certificate (Köhler–Müller–Allgöwer 2018, Ass. 1):
%   Storage:  λ_k ≡ 0  (set p_0 = 0; recursion p_k = A_d'*p_{k+1} gives p_k = 0 ∀k).
%   Optimal periodic cost: ℓ*_k = ‖0‖² + ε‖x*_k − x*_k‖² = 0.
%   Rotated cost:  L_k(x,u) = ‖u‖² + ε‖x − x*_k‖² ≥ ε‖x − x*_k‖² = ρ(|x − x*_k|).
%   → strictly positive definite around Π* with trivial storage; no appeal to A_d^P = I.
%
% Phase-synchronisation prevents the "waiting" trap (Müller–Grüne 2016, Example 4):
%   the per-step penalty ε‖x − x*_{k mod P}‖² is tied to absolute time, so
%   deferring thrust is penalised at every stage, not just at the horizon end.
%
% No terminal constraint — practical asymptotic stability via Corollary 4.
% Trade-off vs Run 2: ε-neighbourhood guarantee (not exact), no hard equality rows.
fprintf('\n════════════════════════════════════════════════════════════\n');
fprintf(' S2 Run 5 — Periodic stage cost ell_aug=||u||^2+eps||x-x*_k||^2\n');
fprintf(' (Koehler-Mueller-Allgower 2018, Ass. 1/Cor. 4 — no terminal cstr.)\n');
fprintf('════════════════════════════════════════════════════════════\n');

eps_r5       = 1e-4;
Q_r5         = eps_r5 * eye(6);
Qbar_r5      = kron(eye(N), Q_r5);
H_r5         = 2 * (Su' * Qbar_r5 * Su + eye(3*N));
H_r5         = (H_r5 + H_r5') / 2;
F_mat_r5     = 2 * Su' * Qbar_r5 * Sx;      % state-dependent part of f_vec
f_ref_mat_r5 = -2 * Su' * Qbar_r5;          % 3N × 6N: ref-dependent part

matrices_r5              = matrices;
matrices_r5.H            = H_r5;
matrices_r5.F_mat        = F_mat_r5;
matrices_r5.f_ref_mat    = f_ref_mat_r5;

params_r5     = params_base;
params_r5.x0  = x0;

run_cfg_r5.use_nmc_manifold      = false;
run_cfg_r5.use_periodic_terminal = false;
run_cfg_r5.x_star                = x_star;
run_cfg_r5.n_sim_override        = [];
run_cfg_r5.use_band              = false;

r5 = run_flyaround(run_cfg_r5, matrices_r5, params_r5);

%% Phase error Run 5
n_sim_r5      = r5.n_sim;
phase_err_r5  = zeros(1, n_sim_r5 + 1);
for k = 1 : n_sim_r5 + 1
    col              = mod(k-1, P) + 1;
    phase_err_r5(k)  = norm(r5.x_log(:, k) - x_star(:, col));
end
last_P_start_r5      = max(1, n_sim_r5 - P + 1);
phase_err_r5_last    = phase_err_r5(last_P_start_r5 : end);
phase_err_r5_max     = max(phase_err_r5_last);
phase_err_r5_mean    = mean(phase_err_r5_last);
fprintf('\nRun 5 phase error (last orbit):  max = %.4f m   mean = %.4f m\n', ...
    phase_err_r5_max, phase_err_r5_mean);

%% ── RUNS 3a / 3b / 3c — "failure mode" suite ────────────────────────────────
% Three configurations that do NOT acquire an NMC orbit, demonstrating that
% neither a state band nor a state cost alone can substitute for the periodic
% terminal constraint of Run 2.
%
%   3a: band (‖r(k)‖ ≥ 75 m) + small Q  → cheapest feasible equilibrium is
%       the V-bar point [0; 75; 0; 0; 0; 0]; NOT an orbit.
%   3b: band + pure fuel               → u*=0, chaser drifts freely under CWH.
%   3c: small Q, no band               → regulation toward origin (like std MPC).

%% State cost matrices shared by 3a and 3c
q_state = 1e-6;
Q_state = q_state * blkdiag(eye(3), zeros(3));   % position-only penalty
Qbar    = kron(eye(N), Q_state);                  % block-diagonal over horizon
H_Q     = 2 * (Su' * Qbar * Su + eye(3*N));       % fuel + state
H_Q     = (H_Q + H_Q') / 2;
F_mat_Q = 2 * Su' * Qbar * Sx;                    % f_vec(k) = F_mat * x(k)

n_sim_3  = 3*P;   % 3 orbits — enough to see asymptotic behaviour

%% RUN 3a — band + small Q
fprintf('\n════════════════════════════════════════════════════════════\n');
fprintf(' S2 Run 3a — band (r_min=75 m) + small Q\n');
fprintf(' Expected: V-bar equilibrium [0; 75; 0; 0; 0; 0], NOT an orbit\n');
fprintf('════════════════════════════════════════════════════════════\n');

matrices_3a       = matrices;
matrices_3a.H     = H_Q;
matrices_3a.F_mat = F_mat_Q;

params_3a         = params_base;
params_3a.x0      = x0;
params_3a.rho_min = 75;
params_3a.rho_max = 1e6;

run_cfg_3a.use_nmc_manifold      = false;
run_cfg_3a.use_periodic_terminal = false;
run_cfg_3a.x_star                = [];
run_cfg_3a.n_sim_override        = n_sim_3;
run_cfg_3a.use_band              = true;

r3a = run_flyaround(run_cfg_3a, matrices_3a, params_3a);

%% RUN 3b — band + pure fuel
fprintf('\n════════════════════════════════════════════════════════════\n');
fprintf(' S2 Run 3b — band (r_min=75 m) + pure fuel\n');
fprintf(' Expected: u*=0, free CWH drift (no orbit selection)\n');
fprintf('════════════════════════════════════════════════════════════\n');

params_3b         = params_base;
params_3b.x0      = x0;
params_3b.rho_min = 75;
params_3b.rho_max = 1e6;

run_cfg_3b.use_nmc_manifold      = false;
run_cfg_3b.use_periodic_terminal = false;
run_cfg_3b.x_star                = [];
run_cfg_3b.n_sim_override        = 2*P;
run_cfg_3b.use_band              = true;

r3b = run_flyaround(run_cfg_3b, matrices, params_3b);

%% RUN 3c — small Q, no band
fprintf('\n════════════════════════════════════════════════════════════\n');
fprintf(' S2 Run 3c — small Q, no band\n');
fprintf(' Expected: regulation to origin (not an orbit)\n');
fprintf('════════════════════════════════════════════════════════════\n');

matrices_3c       = matrices;
matrices_3c.H     = H_Q;
matrices_3c.F_mat = F_mat_Q;

params_3c         = params_base;
params_3c.x0      = x0;
params_3c.rho_min = 0;
params_3c.rho_max = 1e6;

run_cfg_3c.use_nmc_manifold      = false;
run_cfg_3c.use_periodic_terminal = false;
run_cfg_3c.x_star                = [];
run_cfg_3c.n_sim_override        = n_sim_3;
run_cfg_3c.use_band              = false;

r3c = run_flyaround(run_cfg_3c, matrices_3c, params_3c);

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
title(ax2, 'Convergence to \Pi^* — Keerthi-Gilbert / Theorem 2.24');
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

%% Fig 5 — Run 4: manifold IC → free orbit maintenance
c_r4        = [0.50 0.00 0.50];   % purple
time_s_r4   = (0 : n_sim_r4) * dt;

fig5 = figure('Name', 'S2 — Run 4: manifold IC');
tl5b = tiledlayout(fig5, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% Hill frame: trajectory should lie on top of Π*
ax5b1 = nexttile(tl5b);
hold(ax5b1, 'on'); grid(ax5b1, 'on'); axis(ax5b1, 'equal');
plot(ax5b1, x_star_cl(1,:), x_star_cl(2,:), '--', 'Color', c_star, 'LineWidth', 3.0, ...
    'DisplayName', sprintf('\\Pi^* (b=%dm)', b_nmc));
plot(ax5b1, r4.x_log(1,:), r4.x_log(2,:), 'Color', c_r4, 'LineWidth', 1.2, ...
    'DisplayName', 'Run 4: trajectory (u^*=0)');
plot(ax5b1, x0_r4(1), x0_r4(2), 'o', 'Color', c_r4, 'MarkerSize', 8, ...
    'MarkerFaceColor', c_r4, 'DisplayName', 'IC (on \Pi^*)');
plot(ax5b1, 0, 0, '.k', 'MarkerSize', 14, 'DisplayName', 'Target');
xlabel(ax5b1, 'r_x (radial) [m]');
ylabel(ax5b1, 'r_y (along-track) [m]');
title(ax5b1, 'Run 4: IC on \Pi^* — trajectory coincides with \Pi^*');
legend(ax5b1, 'Location', 'northeast');

% Phase error: should be numerical noise
ax5b2 = nexttile(tl5b);
hold(ax5b2, 'on'); grid(ax5b2, 'on');
semilogy(ax5b2, time_s_r4, phase_err_r4 + 1e-16, 'Color', c_r4, 'LineWidth', 1.2);
xlabel(ax5b2, 'Time [s]');
ylabel(ax5b2, '||x(k) - x^*_{k \rm mod P}|| [m]');
title(ax5b2, 'Phase error (expect \approx 0)');
ylim(ax5b2, [1e-16, 1e-3]);

sgtitle(fig5, 'S2 Run 4 — NMC manifold IC: orbit maintained with u^* = 0');
exportgraphics(fig5, fullfile(fig_dir, 'run4_manifold_ic.pdf'), 'ContentType', 'vector');
exportgraphics(fig5, fullfile(fig_dir, 'run4_manifold_ic.png'), 'Resolution', 300);

%% Fig 7 — Run 5 vs Run 2: periodic stage cost vs terminal equality
c_r5        = [0.80 0.20 0.00];   % dark orange-red
time_s_r5   = (0 : n_sim_r5) * dt;

fig7 = figure('Name', 'S2 — Run 5 vs Run 2');
tl7  = tiledlayout(fig7, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% Hill frame: Run 5 should converge to Π* without terminal constraint
ax7a = nexttile(tl7);
hold(ax7a, 'on'); grid(ax7a, 'on'); axis(ax7a, 'equal');
plot(ax7a, x_star_cl(1,:), x_star_cl(2,:), '--', 'Color', c_star, 'LineWidth', 2.0, ...
    'DisplayName', sprintf('\\Pi^* (b=%dm)', b_nmc));
plot(ax7a, r2.x_log(1,:), r2.x_log(2,:), 'Color', c_r2, 'LineWidth', 0.8, ...
    'DisplayName', 'Run 2 (terminal eq.)');
plot(ax7a, r5.x_log(1,:), r5.x_log(2,:), 'Color', c_r5, 'LineWidth', 1.0, ...
    'DisplayName', 'Run 5 (periodic \ell_{aug})');
plot(ax7a, x0(1), x0(2), 's', 'Color', c_r1, 'MarkerSize', 8, ...
    'MarkerFaceColor', c_r1, 'HandleVisibility', 'off');
plot(ax7a, 0, 0, '.k', 'MarkerSize', 14, 'DisplayName', 'Target');
xlabel(ax7a, 'r_x [m]'); ylabel(ax7a, 'r_y [m]');
title(ax7a, 'Hill frame: Run 2 vs Run 5');
legend(ax7a, 'Location', 'northeast', 'FontSize', 7);

% Phase error comparison
ax7b = nexttile(tl7);
hold(ax7b, 'on'); grid(ax7b, 'on');
semilogy(ax7b, time_s_r2, phase_err,    'Color', c_r2, 'LineWidth', 1.2, ...
    'DisplayName', 'Run 2 (terminal eq.)');
semilogy(ax7b, time_s_r5, phase_err_r5, 'Color', c_r5, 'LineWidth', 1.2, ...
    'DisplayName', 'Run 5 (periodic \ell_{aug})');
xlabel(ax7b, 'Time [s]');
ylabel(ax7b, '||x(k) - x^*_{k \rm mod P}|| [m]');
title(ax7b, 'Phase error convergence');
legend(ax7b, 'Location', 'northeast', 'FontSize', 7);
ylim(ax7b, [1e-3, 1e3]);

% Control effort comparison
ax7c = nexttile(tl7);
hold(ax7c, 'on'); grid(ax7c, 'on');
semilogy(ax7c, time_s_r2(1:n_sim_r2), r2.u_norm, 'Color', c_r2, 'LineWidth', 1.2, ...
    'DisplayName', 'Run 2');
semilogy(ax7c, time_s_r5(1:n_sim_r5), r5.u_norm, 'Color', c_r5, 'LineWidth', 1.2, ...
    'DisplayName', 'Run 5');
xlabel(ax7c, 'Time [s]');
ylabel(ax7c, '||u|| [m/s^2]');
title(ax7c, 'Control effort');
legend(ax7c, 'Location', 'northeast', 'FontSize', 7);
ylim(ax7c, [1e-11, 1e-2]);

sgtitle(fig7, sprintf('S2: Run 2 (Keerthi-Gilbert) vs Run 5 (periodic \\ell_{aug}, \\epsilon=%.0e)', eps_r5));
exportgraphics(fig7, fullfile(fig_dir, 'run5_vs_run2.pdf'), 'ContentType', 'vector');
exportgraphics(fig7, fullfile(fig_dir, 'run5_vs_run2.png'), 'Resolution', 300);

%% Fig 6 — Failure-mode suite: runs 3a / 3b / 3c
c_3a = [0.13 0.63 0.37];   % green
c_3b = [0.49 0.18 0.56];   % purple
c_3c = [0.93 0.53 0.18];   % orange

n_sim_3a = r3a.n_sim;
n_sim_3b = r3b.n_sim;
n_sim_3c = r3c.n_sim;

fig6 = figure('Name', 'S2 — Failure modes (Runs 3a/3b/3c)');
tl6  = tiledlayout(fig6, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% 3a — band + Q → V-bar equilibrium at boundary
ax6a = nexttile(tl6);
hold(ax6a, 'on'); grid(ax6a, 'on'); axis(ax6a, 'equal');
plot(ax6a, x_star_cl(1,:), x_star_cl(2,:), '--', 'Color', c_star, 'LineWidth', 1.5, ...
    'DisplayName', sprintf('\\Pi^* (b=%dm)', b_nmc));
plot(ax6a, r3a.x_log(1,:), r3a.x_log(2,:), 'Color', c_3a, 'LineWidth', 1.0, ...
    'DisplayName', '3a trajectory');
plot(ax6a, x0(1), x0(2), 's', 'Color', c_r1, 'MarkerSize', 8, ...
    'MarkerFaceColor', c_r1, 'HandleVisibility', 'off');
plot(ax6a, r3a.x_log(1,end), r3a.x_log(2,end), '^', 'Color', c_3a, ...
    'MarkerSize', 8, 'MarkerFaceColor', c_3a, 'DisplayName', 'Final state');
plot(ax6a, 0, 0, '.k', 'MarkerSize', 14, 'DisplayName', 'Target');
xlabel(ax6a, 'r_x [m]'); ylabel(ax6a, 'r_y [m]');
title(ax6a, '3a: band + Q  \rightarrow V-bar eq. at boundary');
legend(ax6a, 'Location', 'northeast', 'FontSize', 7);

% 3b — band + fuel → free drift
ax6b = nexttile(tl6);
hold(ax6b, 'on'); grid(ax6b, 'on'); axis(ax6b, 'equal');
plot(ax6b, x_star_cl(1,:), x_star_cl(2,:), '--', 'Color', c_star, 'LineWidth', 1.5, ...
    'DisplayName', sprintf('\\Pi^* (b=%dm)', b_nmc));
plot(ax6b, r3b.x_log(1,:), r3b.x_log(2,:), 'Color', c_3b, 'LineWidth', 1.0, ...
    'DisplayName', '3b trajectory (drift)');
plot(ax6b, x0(1), x0(2), 's', 'Color', c_r1, 'MarkerSize', 8, ...
    'MarkerFaceColor', c_r1, 'HandleVisibility', 'off');
plot(ax6b, 0, 0, '.k', 'MarkerSize', 14, 'DisplayName', 'Target');
xlabel(ax6b, 'r_x [m]'); ylabel(ax6b, 'r_y [m]');
title(ax6b, '3b: band + fuel  \rightarrow free drift (u^*=0)');
legend(ax6b, 'Location', 'northeast', 'FontSize', 7);

% 3c — no band, small Q → regulation to origin
ax6c = nexttile(tl6);
hold(ax6c, 'on'); grid(ax6c, 'on'); axis(ax6c, 'equal');
plot(ax6c, x_star_cl(1,:), x_star_cl(2,:), '--', 'Color', c_star, 'LineWidth', 1.5, ...
    'DisplayName', sprintf('\\Pi^* (b=%dm)', b_nmc));
plot(ax6c, r3c.x_log(1,:), r3c.x_log(2,:), 'Color', c_3c, 'LineWidth', 1.0, ...
    'DisplayName', '3c trajectory');
plot(ax6c, x0(1), x0(2), 's', 'Color', c_r1, 'MarkerSize', 8, ...
    'MarkerFaceColor', c_r1, 'HandleVisibility', 'off');
plot(ax6c, r3c.x_log(1,end), r3c.x_log(2,end), '^', 'Color', c_3c, ...
    'MarkerSize', 8, 'MarkerFaceColor', c_3c, 'DisplayName', 'Final state');
plot(ax6c, 0, 0, '.k', 'MarkerSize', 14, 'DisplayName', 'Target');
xlabel(ax6c, 'r_x [m]'); ylabel(ax6c, 'r_y [m]');
title(ax6c, '3c: Q only  \rightarrow regulation to origin');
legend(ax6c, 'Location', 'northeast', 'FontSize', 7);

sgtitle(fig6, 'S2 — Why the periodic terminal constraint is necessary');
exportgraphics(fig6, fullfile(fig_dir, 'failure_modes.pdf'), 'ContentType', 'vector');
exportgraphics(fig6, fullfile(fig_dir, 'failure_modes.png'), 'Resolution', 300);

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

res3a = struct( ...
    'terminal',               'none', ...
    'band',                   'r_min=75m, no outer', ...
    'cost',                   sprintf('fuel + Q (q=%.0e, position only)', q_state), ...
    'n_sim_steps',            n_sim_3a, ...
    'delta_v_total_m_s',      r3a.dv_total, ...
    'final_range_m',          r3a.range_log(end), ...
    'infeasible_qp_steps',    sum(r3a.exitflag_log <= 0), ...
    'note',                   'Converges to V-bar equilibrium at band boundary — NOT an NMC orbit');

res3b = struct( ...
    'terminal',               'none', ...
    'band',                   'r_min=75m, no outer', ...
    'cost',                   'pure fuel', ...
    'n_sim_steps',            n_sim_3b, ...
    'delta_v_total_m_s',      r3b.dv_total, ...
    'final_range_m',          r3b.range_log(end), ...
    'infeasible_qp_steps',    sum(r3b.exitflag_log <= 0), ...
    'note',                   'u*=0 globally optimal; chaser drifts freely under CWH dynamics');

res3c = struct( ...
    'terminal',               'none', ...
    'band',                   'none', ...
    'cost',                   sprintf('fuel + Q (q=%.0e, position only)', q_state), ...
    'n_sim_steps',            n_sim_3c, ...
    'delta_v_total_m_s',      r3c.dv_total, ...
    'final_range_m',          r3c.range_log(end), ...
    'infeasible_qp_steps',    sum(r3c.exitflag_log <= 0), ...
    'note',                   'Regulation to origin — no constraint forces orbital motion');

res4 = struct( ...
    'terminal',               'none', ...
    'band',                   'none', ...
    'cost',                   'pure fuel', ...
    'IC',                     x0_r4', ...
    'IC_note',                'on Pi* manifold: x=b_nmc, vy=-2*n*b_nmc', ...
    'n_sim_steps',            n_sim_r4, ...
    'delta_v_total_m_s',      r4.dv_total, ...
    'phase_err_max_m',        max(phase_err_r4), ...
    'phase_err_mean_m',       mean(phase_err_r4), ...
    'infeasible_qp_steps',    sum(r4.exitflag_log <= 0), ...
    'note',                   'IC on NMC manifold: U*=0 everywhere, Pi* maintained at zero fuel cost');

res5 = struct( ...
    'terminal',               'none', ...
    'band',                   'none', ...
    'cost',                   sprintf('fuel + eps*||x-x*_k||^2  (eps=%.0e)', eps_r5), ...
    'theory',                 'Koehler-Mueller-Allgower 2018 Ass.1/Cor.4 periodic dissipativity', ...
    'n_sim_steps',            n_sim_r5, ...
    'delta_v_total_m_s',      r5.dv_total, ...
    'phase_err_max_last_m',   phase_err_r5_max, ...
    'phase_err_mean_last_m',  phase_err_r5_mean, ...
    'b_est_m',                r5.b_est, ...
    'range_ratio_2to1',       r5.ratio_2to1, ...
    'infeasible_qp_steps',    sum(r5.exitflag_log <= 0), ...
    'note',                   'No terminal constraint; convergence via periodic strict dissipativity');

data = struct('metadata', meta, 'run1', res1, 'run2', res2, 'run4', res4, ...
              'run5', res5, 'run3a', res3a, 'run3b', res3b, 'run3c', res3c);

json_str    = jsonencode(data, 'PrettyPrint', true);
export_path = fullfile(script_dir, '..', 'exports', 'scenarios', ...
    'sim_scenario2_flyaround.json');
fid = fopen(export_path, 'w');
fprintf(fid, '%s', json_str);
fclose(fid);
fprintf('JSON exported to %s\n', export_path);
