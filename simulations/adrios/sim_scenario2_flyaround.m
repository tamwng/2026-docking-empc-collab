% sim_scenario2_flyaround.m
% Scenario 2: NMC Fly-Around (IC: V-bar hold [0;300;0;0;0;0])
%
% Run 1  - No terminal constraint: u*=0, chaser stays at V-bar hold (CWH equilibrium).
% Run 2  - Periodic terminal equality x(N)=Π*(k+N-1 mod P): converges to NMC orbit.
% Run 3a - Band ‖r‖≥75 m + small Q: converges to V-bar equilibrium, not an orbit.
% Run 3b - Band + pure fuel: u*=0, free CWH drift (no orbit selection).
% Run 3c - Small Q, no band: regulation toward origin (like std MPC, not orbit).
% Run 4  - IC on Π* manifold, pure fuel, no terminal constraint: u*=0, free orbit maintenance.
% Run 5  - Augmented stage cost ℓ=‖u‖²+ε‖x−x*_k‖²: strict P-periodic dissipativity,
%          no terminal constraint (Köhler-Müller-Allgöwer 2018).

clear; clc;
script_dir = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(script_dir, '..', '..', 'tools', 'matlab2tikz-master', 'src')));
constants;  

%%  Discretisation
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

%check
x_test   = [100; 0; 0; 0; -2*n*100; 0];
residual = norm(A_d^P * x_test - x_test);
fprintf('NMT closure residual: %.2e  (must be < 1e-4)\n', residual);
assert(residual < 1e-4, 'Discrete NMT does not close: check dt and P');

%% Prescribed NMC orbit Π* 
% 2:1 ellipse in Hill frame: x in [-b, b] radial, y in [-2b, 2b] along-track.
% IC on NMC manifold: vx=0, vy+2n*x=0 -> start at [b;0;0;0;-2nb;0].
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



%% QP dimensions and parameters
u_max = 1e-2;    % [m/s²]  symmetric per-axis thrust bound
N     = P;       % horizon = one full orbital period
n_sim = 6*P;     % total simulation length
n_x   = 6;
n_u   = 3;

%% Condensed prediction matrices Sx (6N×6) and Su (6N×3N)
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

% Position-row extracts
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

x0 = [0; 300; 0; 0; 0; 0];   % V-bar hold, CWH equilibrium (shared IC)

params_base.rho_min = 50;    
params_base.rho_max = 200;
params_base.u_max   = u_max;
params_base.N       = N;
params_base.P       = P;
params_base.n       = n;
params_base.dt      = dt;
params_base.n_sim   = n_sim;

%% RUN 1 - No terminal constraint

fprintf(' S2 Run 1 - No terminal constraint (V-bar hold, chaser stays put)\n');

params_r1     = params_base;
params_r1.x0  = x0;

run_cfg_r1.use_nmc_manifold      = false;
run_cfg_r1.use_periodic_terminal = false;
run_cfg_r1.x_star                = [];
run_cfg_r1.n_sim_override        = 2*P;
run_cfg_r1.use_band              = false;

r1 = run_flyaround(run_cfg_r1, matrices, params_r1);

%% RUN 2 - Periodic terminal to Π*
fprintf(' S2 Run 2 - Periodic terminal constraint -> stabilise Pi*\n');
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

% Phase error 
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

%%  RUN 4 - NMC manifold IC, pure fuel, no terminal constraint
% IC is placed exactly on Π* (satisfies vy = -2n·x, vx = 0).

fprintf(' S2 Run 4 - IC on NMC manifold, pure fuel, no terminal constraint\n');

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

% Phase error Run 4 
n_sim_r4      = r4.n_sim;
phase_err_r4  = zeros(1, n_sim_r4 + 1);
for k = 1 : n_sim_r4 + 1
    col              = mod(k-1, P) + 1;
    phase_err_r4(k)  = norm(r4.x_log(:, k) - x_star(:, col));
end
fprintf('Run 4 phase error:  max = %.2e m   mean = %.2e m  (expect machine-eps)\n', ...
    max(phase_err_r4), mean(phase_err_r4));

%% RUN 5 - Regularised dissipativity route (periodic stage cost)
% Stage cost: ℓ_aug(x,u,k) = ‖u‖² + ε‖x − x*_k‖²
% Phase-synchronisation prevents the "waiting" trap (Müller-Grüne 2016, Example 4):

fprintf(' S2 Run 5 - Periodic stage cost ell_aug=||u||^2+eps||x-x*_k||^2\n');
fprintf(' (Koehler-Mueller-Allgower 2018, Ass. 1/Cor. 4, no terminal cstr.)\n');


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

% Phase error Run 5
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

%% RUN T - Standard Tracking MPC baseline (quadratic cost, tracks Pi*)
fprintf(' S2 Run T - Standard Tracking MPC (tracks Pi*, quadratic cost)\n');

Q_track = eye(6);          % penalise full state error to the reference
R_track = eye(3);          % control effort (same scale as R_eco in Runs 2/5)
[~, P_track] = lqr_controller(A_d, B_d, Q_track, R_track);

% Short horizon on purpose: over a full period (N=P=92) Su'*Q*Su is severely
% ill-conditioned with Q=I, stalling quadprog. N_track=20 keeps the QP
% well-conditioned and matches Scenarios 1 & 4 (Runs 2/5 tolerate N=P only
% because their fuel-dominant Hessian stays H≈2I).
N_track = 20;

con_track.u_max = u_max;
con_track.u_ss  = zeros(3, 1);      % NMC is zero-input: suppress feedforward

x_track = zeros(6, n_sim + 1);
u_track = zeros(3, n_sim);
x_track(:, 1) = x0;
for k = 1:n_sim
    cols   = mod(k - 1 + (1:N_track), P) + 1;    % phase-synced N_track-step preview
    R_prev = x_star(:, cols);                     % 6 x N_track
    u_k    = mpc_tracking(x_track(:, k), R_prev, A_d, B_d, ...
                          Q_track, R_track, P_track, N_track, con_track);
    u_track(:, k)     = u_k;
    x_track(:, k + 1) = A_d * x_track(:, k) + B_d * u_k;
end

dv_track    = sum(vecnorm(u_track, 2, 1)) * dt;
phase_err_T = zeros(1, n_sim + 1);
for k = 1:n_sim + 1
    col            = mod(k - 1, P) + 1;
    phase_err_T(k) = norm(x_track(:, k) - x_star(:, col));
end

%% Capture-time comparison: Tracking vs Run 2 vs Run 5 (shared IC)
capture_tol = 5;     % [m] phase-error threshold defining "orbit captured"
cap_T  = capture_step(phase_err_T,  capture_tol);
cap_r2 = capture_step(phase_err,    capture_tol);
cap_r5 = capture_step(phase_err_r5, capture_tol);
dv_r2  = r2.dv_total;
dv_r5  = r5.dv_total;

fprintf('\n%s\n', repmat('-', 1, 70));
fprintf('Tracking vs Run 2 vs Run 5  (shared IC, capture tol = %d m)\n', capture_tol);
fprintf('%s\n', repmat('-', 1, 70));
fprintf('  %-26s  %9s  %9s  %13s\n', 'Controller', 'dv[m/s]', 'vs track', 'capture');
fprintf('  %-26s  %9.4f  %9s  %13s\n', 'Standard Tracking MPC', dv_track, '--', cap_str(cap_T, P));
fprintf('  %-26s  %9.4f  %8.1f%%  %13s\n', 'Run 2 (terminal eq.)', dv_r2, ...
    100*(dv_track - dv_r2)/dv_track, cap_str(cap_r2, P));
fprintf('  %-26s  %9.4f  %8.1f%%  %13s\n', 'Run 5 (periodic aug.)', dv_r5, ...
    100*(dv_track - dv_r5)/dv_track, cap_str(cap_r5, P));
fprintf('%s\n\n', repmat('-', 1, 70));

%% RUNS 3a / 3b / 3c
% Three configurations that do NOT acquire an NMC orbit, demonstrating that
% neither a state band nor a state cost alone can substitute for the periodic
% terminal constraint of Run 2.
%
%   3a: band (‖r(k)‖ ≥ 75 m) + small Q  -> cheapest feasible equilibrium is
%       the V-bar point [0; 75; 0; 0; 0; 0]; NOT an orbit.
%   3b: band + pure fuel               -> u*=0, chaser drifts freely under CWH.
%   3c: small Q, no band               -> regulation toward origin (like std MPC).

%% State cost matrices shared by 3a and 3c
q_state = 1e-6;
Q_state = q_state * blkdiag(eye(3), zeros(3));   % position-only penalty
Qbar    = kron(eye(N), Q_state);                  % block-diagonal over horizon
H_Q     = 2 * (Su' * Qbar * Su + eye(3*N));       % fuel + state
H_Q     = (H_Q + H_Q') / 2;
F_mat_Q = 2 * Su' * Qbar * Sx;                    % f_vec(k) = F_mat * x(k)

n_sim_3  = 3*P;   % 3 orbits, enough to see asymptotic behaviour

%% RUN 3a - band + small Q

fprintf(' S2 Run 3a - band (r_min=75 m) + small Q\n');
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

%% RUN 3b - band + pure fuel

fprintf(' S2 Run 3b - band (r_min=75 m) + pure fuel\n');


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

%% RUN 3c - small Q, no band

fprintf(' S2 Run 3c - small Q, no band\n');

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

n_sim_3a = r3a.n_sim;
n_sim_3b = r3b.n_sim;
n_sim_3c = r3c.n_sim;

%% FIGURES
fig_dir = fullfile(script_dir, '..', '..', 'results', 'figures', 's2_flyaround');
if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end

c_r2   = [0.00 0.45 0.70];
c_star = [0.93 0.69 0.13];
c_r1   = [0.85 0.33 0.10];
c_r5   = [0.80 0.20 0.00];   % dark orange-red

n_sim_r1  = r1.n_sim;
time_s_r2 = (0 : n_sim_r2) * dt;
x_star_cl = [x_star, x_star(:,1)];   % close the curve for plotting

% Fig 8: Standard Tracking MPC vs Run 2 vs Run 5: Hill-frame overlay + table
c_track = [0.55 0.10 0.75];   % purple, standard tracking baseline

% Presentation-facing run labels (thesis/slide terminology, not the code's
% internal run names), matching the Operation 1 (closing) convention.
lbl_MPC_S2 = 'Std.\ MPC';
lbl_R2_S2  = 'Per.\ term.\ eq.';
lbl_R5_S2  = 'Phase synchr.\ $\ell_{\mathrm{aug}}$';

% Shared axis-box aspect ratio for fig8 and fig6b (both datasets are
% already available here), so that scaling both .tikz files to the same
% \linewidth in LaTeX renders them at the same physical size. Target ratio
% is fig6b's own data proportions (its Pi* ellipse reads more naturally
% there); only the under-sized dimension of each figure is expanded to
% match it, so neither figure's data gets clipped.
all_x8_raw  = [x_star_cl(1,:), x_track(1,:), r2.x_log(1,:), r5.x_log(1,:), x0(1), 0];
all_y8_raw  = [x_star_cl(2,:), x_track(2,:), r2.x_log(2,:), r5.x_log(2,:), x0(2), 0];
all_x6b_raw = [x_star_cl(1,:), r3a.x_log(1,:), r3b.x_log(1,:), r3c.x_log(1,:), x0(1), 0];
all_y6b_raw = [x_star_cl(2,:), r3a.x_log(2,:), r3b.x_log(2,:), r3c.x_log(2,:), x0(2), 0];

pad_frac_s2 = 0.08;
rx8  = (max(all_x8_raw)  - min(all_x8_raw))  * (1 + 2*pad_frac_s2);
ry8  = (max(all_y8_raw)  - min(all_y8_raw))  * (1 + 2*pad_frac_s2);
rx6b = (max(all_x6b_raw) - min(all_x6b_raw)) * (1 + 2*pad_frac_s2);
ry6b = (max(all_y6b_raw) - min(all_y6b_raw)) * (1 + 2*pad_frac_s2);
target_ratio_s2 = ry6b / rx6b;   % height/width, from the failure-mode plot

if ry8 >= rx8 * target_ratio_s2
    rx8_f = ry8 / target_ratio_s2; ry8_f = ry8;
else
    rx8_f = rx8; ry8_f = rx8 * target_ratio_s2;
end
cx8 = (max(all_x8_raw) + min(all_x8_raw)) / 2;
cy8 = (max(all_y8_raw) + min(all_y8_raw)) / 2;
xlim8_final = [cx8 - rx8_f/2, cx8 + rx8_f/2];
ylim8_final = [cy8 - ry8_f/2, cy8 + ry8_f/2];

if ry6b >= rx6b * target_ratio_s2
    rx6b_f = ry6b / target_ratio_s2; ry6b_f = ry6b;
else
    rx6b_f = rx6b; ry6b_f = rx6b * target_ratio_s2;
end
cx6b = (max(all_x6b_raw) + min(all_x6b_raw)) / 2;
cy6b = (max(all_y6b_raw) + min(all_y6b_raw)) / 2;
xlim6b_final = [cx6b - rx6b_f/2, cx6b + rx6b_f/2];
ylim6b_final = [cy6b - ry6b_f/2, cy6b + ry6b_f/2];

fig8 = figure('Name', 'Op 2 - Tracking vs Run 2 vs Run 5');
ax8  = axes(fig8); hold(ax8, 'on'); grid(ax8, 'on'); axis(ax8, 'equal');
plot(ax8, x_star_cl(1,:), x_star_cl(2,:), '--', 'Color', c_star, 'LineWidth', 2.0, ...
    'DisplayName', sprintf('$\\Pi^*$ ($b=%d$\\,m)', b_nmc));
% Distinct linestyles (not just colors) for these three: they can sit close
% enough to each other that color alone isn't enough to tell them apart.
plot(ax8, x_track(1,:), x_track(2,:), '-', 'Color', c_track, 'LineWidth', 1.4, ...
    'DisplayName', lbl_MPC_S2);
plot(ax8, r2.x_log(1,:), r2.x_log(2,:), '--', 'Color', c_r2, 'LineWidth', 1.4, ...
    'DisplayName', lbl_R2_S2);
plot(ax8, r5.x_log(1,:), r5.x_log(2,:), '-.', 'Color', c_r5, 'LineWidth', 1.4, ...
    'DisplayName', lbl_R5_S2);
plot(ax8, x0(1), x0(2), 's', 'Color', [0.4 0.4 0.4], 'MarkerFaceColor', [0.4 0.4 0.4], ...
    'MarkerSize', 8, 'DisplayName', 'IC (shared)');
plot(ax8, 0, 0, '.k', 'MarkerSize', 14, 'DisplayName', 'Target');
% NOTE: this plot's own convention is radial=horizontal, along-track=
% vertical (x_log(1,:) plotted first) -- the OPPOSITE of S1/S4 and the
% validity-chapter hillframe plots (along-track horizontal, radial
% vertical). Only the label text was changed here, not the data order.
xlabel(ax8, 'Radial $x$ [m]', 'Interpreter', 'latex');
ylabel(ax8, 'Along-track $y$ [m]', 'Interpreter', 'latex');
% Inside the axis, upper-left: verified data-driven (zero trajectory points
% fall within that corner's footprint, vs. hundreds in both lower corners).
legend(ax8, 'Location', 'northwest', 'FontSize', 8, 'Interpreter', 'latex');
% Widen further (left side, where the now-inset legend sits) so the legend
% box fits inside the axis without overflowing -- doesn't need to be as
% wide as the pre-aspect-matching version, just enough for the text.
xlim8_final(1) = xlim8_final(1) - 90;
xlim(ax8, xlim8_final);
ylim(ax8, ylim8_final);

exportgraphics(fig8, fullfile(fig_dir, 'tracking_vs_run2_run5.pdf'), 'ContentType', 'vector');
exportgraphics(fig8, fullfile(fig_dir, 'tracking_vs_run2_run5.png'), 'Resolution', 300);
if exist('matlab2tikz', 'file')
    tikz_path8 = fullfile(fig_dir, 'tracking_vs_run2_run5.tikz');
    matlab2tikz(tikz_path8, 'figurehandle', fig8, 'showInfo', false, 'parseStrings', false);
    % legend(...,'FontSize',...) is silently dropped by this matlab2tikz
    % version -- patch directly for "a little bigger".
    txt8 = fileread(tikz_path8);
    txt8 = regexprep(txt8, 'legend style=\{', 'legend style={font=\\tiny, ', 'once');
    fid8 = fopen(tikz_path8, 'w'); fprintf(fid8, '%s', txt8); fclose(fid8);
end

% LaTeX comparison table (Tracking vs Run 2 vs Run 5)
fid_ct = fopen(fullfile(fig_dir, 'tracking_comparison_table.tex'), 'w');
fprintf(fid_ct, '%% Auto-generated by sim_scenario2_flyaround.m\n');
fprintf(fid_ct, '\\begin{tabular}{l r r r}\n\\toprule\n');
fprintf(fid_ct, ['Controller & $\\Delta v$ [m/s] & saving vs track.\\ [\\%%] & ', ...
    'capture [orbits] \\\\\n\\midrule\n']);
fprintf(fid_ct, 'Standard Tracking MPC & %.4f & -- & %s \\\\\n', dv_track, cap_tex(cap_T, P));
fprintf(fid_ct, 'Run 2 (terminal eq.) & %.4f & %.1f & %s \\\\\n', ...
    dv_r2, 100*(dv_track - dv_r2)/dv_track, cap_tex(cap_r2, P));
fprintf(fid_ct, 'Run 5 (periodic $\\ell_{\\mathrm{aug}}$) & %.4f & %.1f & %s \\\\\n', ...
    dv_r5, 100*(dv_track - dv_r5)/dv_track, cap_tex(cap_r5, P));
fprintf(fid_ct, '\\bottomrule\n\\end{tabular}\n');
fclose(fid_ct);

% Fig 6b: Failure modes overlaid in a single Hill-frame (3a/3b/3c)
c_3a = [0.13 0.63 0.37];
c_3b = [0.49 0.18 0.56];
c_3c = [0.93 0.53 0.18];
fig6b  = figure('Name', 'Op 2 - Failure modes overlay (3a/3b/3c)');
ax6bb  = axes(fig6b); hold(ax6bb, 'on'); grid(ax6bb, 'on'); axis(ax6bb, 'equal');
plot(ax6bb, x_star_cl(1,:), x_star_cl(2,:), '--', 'Color', c_star, 'LineWidth', 1.8, ...
    'DisplayName', sprintf('$\\Pi^*$ ($b=%d$\\,m)', b_nmc));
plot(ax6bb, r3a.x_log(1,:), r3a.x_log(2,:), '-', 'Color', c_3a, 'LineWidth', 1.4, ...
    'DisplayName', 'band $+\,Q$');
plot(ax6bb, r3b.x_log(1,:), r3b.x_log(2,:), '-', 'Color', c_3b, 'LineWidth', 1.4, ...
    'DisplayName', 'band $+$ fuel');
plot(ax6bb, r3c.x_log(1,:), r3c.x_log(2,:), '-', 'Color', c_3c, 'LineWidth', 1.4, ...
    'DisplayName', '$Q$ only');
plot(ax6bb, x0(1), x0(2), 's', 'Color', [0.4 0.4 0.4], 'MarkerFaceColor', [0.4 0.4 0.4], ...
    'MarkerSize', 8, 'DisplayName', 'IC (shared)');
plot(ax6bb, 0, 0, '.k', 'MarkerSize', 14, 'DisplayName', 'Target');
% Same axis-order note as tracking_vs_run2_run5.tikz above: radial is the
% horizontal (plotted first) here, opposite of S1/S4's convention.
xlabel(ax6bb, 'Radial $x$ [m]', 'Interpreter', 'latex');
ylabel(ax6bb, 'Along-track $y$ [m]', 'Interpreter', 'latex');
legend(ax6bb, 'Location', 'northwest', 'FontSize', 8, 'Interpreter', 'latex');
% Widen left side (where the legend now sits) so it fits without
% overflowing the axis, matching the same treatment as fig8 above.
xlim6b_final(1) = xlim6b_final(1) - 70;
xlim(ax6bb, xlim6b_final);
ylim(ax6bb, ylim6b_final);

exportgraphics(fig6b, fullfile(fig_dir, 'failure_modes_overlay.pdf'), 'ContentType', 'vector');
exportgraphics(fig6b, fullfile(fig_dir, 'failure_modes_overlay.png'), 'Resolution', 300);
if exist('matlab2tikz', 'file')
    tikz_path6b = fullfile(fig_dir, 'failure_modes_overlay.tikz');
    matlab2tikz(tikz_path6b, 'figurehandle', fig6b, 'showInfo', false, 'parseStrings', false);
    txt6b = fileread(tikz_path6b);
    txt6b = regexprep(txt6b, 'legend style=\{', 'legend style={font=\\tiny, ', 'once');
    fid6b = fopen(tikz_path6b, 'w'); fprintf(fid6b, '%s', txt6b); fclose(fid6b);
end

%% JSON EXPORT
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
    'note',                   'CWH equilibrium; U*=0 globally optimal, chaser stays at IC');

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
    'note',                   'Converges to V-bar equilibrium at band boundary, NOT an NMC orbit');

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
    'note',                   'Regulation to origin, no constraint forces orbital motion');

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

resT = struct( ...
    'controller',             'standard_tracking_mpc', ...
    'cost',                   'e''Qe + u''Ru  (Q=I6, R=I3, P=DARE)', ...
    'IC',                     x0', ...
    'n_sim_steps',            n_sim, ...
    'delta_v_total_m_s',      dv_track, ...
    'dv_saving_run2_pct',     100*(dv_track - dv_r2)/dv_track, ...
    'dv_saving_run5_pct',     100*(dv_track - dv_r5)/dv_track, ...
    'capture_orbits',         (cap_T  - 1)/P, ...
    'capture_orbits_run2',    (cap_r2 - 1)/P, ...
    'capture_orbits_run5',    (cap_r5 - 1)/P, ...
    'capture_tol_m',          capture_tol, ...
    'phase_err_last_m',       phase_err_T(end), ...
    'note',                   'Standard tracking MPC baseline for Runs 2 & 5; shared IC x0');

data = struct('metadata', meta, 'run1', res1, 'run2', res2, 'run4', res4, ...
              'run5', res5, 'runT', resT, 'run3a', res3a, 'run3b', res3b, 'run3c', res3c);

json_str    = jsonencode(data, 'PrettyPrint', true);
export_path = fullfile(script_dir, '..', '..', 'exports', 'scenarios', ...
    'sim_scenario2_flyaround.json');
fid = fopen(export_path, 'w');
fprintf(fid, '%s', json_str);
fclose(fid);
fprintf('JSON exported to %s\n', export_path);

%% Local functions

function ks = capture_step(perr, tol)
% First step index from which the phase error stays at/below tol for the
% remainder of the run (i.e. the orbit is captured and held). NaN if never.
    below = perr <= tol;
    ks    = NaN;
    for k = 1:numel(below)
        if all(below(k:end)); ks = k; break; end
    end
end

function s = cap_str(k, P)
% Console string for a capture step, expressed in orbital periods.
    if isnan(k); s = 'not captured'; else; s = sprintf('%.2f orb', (k-1)/P); end
end

function s = cap_tex(k, P)
% LaTeX-table string for a capture step, in orbital periods.
    if isnan(k); s = 'n/a'; else; s = sprintf('%.2f', (k-1)/P); end
end
