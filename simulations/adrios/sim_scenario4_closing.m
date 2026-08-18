% sim_scenario4_closing.m
% Scenario 4 — Three EMPC stabilisation mechanisms on the closing phase.
%
% Run A1     — pure fuel (Q=0, P=0), terminal ball ON: converges.
% Run A2     — same, ball OFF: u≡0, chaser stuck
% Run B-noball — fuel + small state penalty ℓ=‖u‖²+ρc‖e‖², no ball: converges
%               via strict dissipativity
% Run B-ball — same cost, ball ON: confirms ball is inactive.
% Run C      — pure fuel, CLF terminal cost Vf=e'P_clf e, no ball: asymptotically
%               stable without strict dissipativity (Amrit et al. 2011).
% Run A1-ps  — A1 config + passive safety SCA: verifies ps_frac ≈ 100%.
% Run MPC    — standard regulation MPC (Q,R≠0, LQR terminal cost): fuel-hungry baseline.
% Sweeps     — horizon/ρc tradeoff (B-noball), ball radius (A1),
%              passive safety vs N and vs r_KOS (Run C).

clear; clc;

%% Path setup
script_dir = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(script_dir, '..', '..', 'tools', 'matlab2tikz-master', 'src')));

%% Constants
constants;

%% Shared scenario parameters
dt        = 120;                    
N         = 20;                       
x0        = [0; 10000; 0; 0; 0; 0]; 
x_switch  = [0;   300; 0; 0; 0; 0];  % [m] handover point (V-bar, 300 m)
u_max     = 1e-2;                     % [m/s²] thrust bound
r_KOS     = 50;                       % [m] keep-out sphere radius
P_avail   = 1e-5;                     % [m²/s⁴] avg-power budget
rho_c     = 1e-11;                     % state regularisation weight 
ps_margin = 150;                       % [m] passive-safety activation margin
t_final   = 5 * T;                    % [s] ~7.7 h (5 orbital periods)

%% Cost structs
cost_A = struct('Q', zeros(6),    'R', eye(3), 'P', zeros(6));
cost_B = struct('Q', rho_c*eye(6),'R', eye(3), 'P', zeros(6));

%% Precompute CLF terminal cost (Run C)
Ac_cwh = [0      0     0    1     0    0  ;
          0      0     0    0     1    0  ;
          0      0     0    0     0    1  ;
          3*n^2  0     0    0     2*n  0  ;
          0      0     0   -2*n   0    0  ;
          0      0    -n^2  0     0    0 ];
Bc_cwh = [zeros(3,3); eye(3)];
[Ad_clf, Bd_clf] = discretize(Ac_cwh, Bc_cwh, dt);

Q_dare = diag([1e-2, 1e-2, 1e-2, 1e0, 1e0, 1e0]);
R_dare = diag([1e4,  1e4,  1e4]);
[K_clf, P_clf] = lqr_controller(Ad_clf, Bd_clf, Q_dare, R_dare);
R_eco_clf      = eye(3);

% Verify Amrit et al. (2011) Assumption 6 at dt=120s.
% M1 = P - Acl'*P*Acl - K'*R_eco*K = Q_dare + K'*(R_dare - R_eco)*K >= 0.
Acl_clf = Ad_clf - Bd_clf * K_clf;
M1      = Q_dare + K_clf' * (R_dare - R_eco_clf) * K_clf;
M1      = (M1 + M1') / 2;
lam_M1  = min(eig(M1));
rho_Acl = max(abs(eig(Acl_clf)));
if lam_M1 > 1e-10 && rho_Acl < 1
    fprintf('  => P_clf satisfies Assumption 6. Terminal cost is a valid CLF.  ✓\n\n');
else
    warning('CLF condition not satisfied — check DARE weights.');
end

cost_C = struct('Q', zeros(6), 'R', eye(3), 'P', P_clf);

%% Shared params base
params_base.x0       = x0;
params_base.x_switch = x_switch;
params_base.dt       = dt;
params_base.N        = N;
params_base.u_max    = u_max;
params_base.t_final  = t_final;
params_base.r_KOS    = r_KOS;

%% Run A1: pure energy, terminal ball ON

con_A1.u_max                  = u_max;
con_A1.y_min_active           = false;
con_A1.los_cone.active        = false;
con_A1.use_avg_power          = false;
con_A1.avg_power.P_avail      = P_avail;
con_A1.use_terminal_ball      = true;
con_A1.terminal_ball.r_switch = 200;   % [m] L-inf position ball
con_A1.use_passive_safety     = false;

params_A1 = params_base;

fprintf(' Run A1 — pure energy ℓ=‖u‖², ball ON\n');
resA1 = run_closing(cost_A, con_A1, params_A1);

%% Run A2: pure energy, terminal ball OFF

con_A2.u_max              = u_max;
con_A2.y_min_active       = false;
con_A2.los_cone.active    = false;
con_A2.use_avg_power      = false;
con_A2.avg_power.P_avail  = P_avail;
con_A2.use_terminal_ball  = false;
con_A2.use_passive_safety = false;

params_A2                  = params_base;
params_A2.n_steps_override = 60;

fprintf(' Run A2 — pure energy ℓ=‖u‖², ball OFF\n');
resA2 = run_closing(cost_A, con_A2, params_A2);

max_disp_A2 = max(vecnorm(resA2.X(1:3,:) - x0(1:3), 2, 1));
max_u_A2    = max([resA2.u_norm_seq, 0]);
fprintf('A2 check — max position displacement: %.4f m   (expect ≈ 0)\n', max_disp_A2);
fprintf('         — max ||u||:                 %.2e m/s² (expect ≈ 0)\n\n', max_u_A2);

%% Run B-noball: strict dissipativity, no terminal ball
% Adding ρc‖e‖² restores strict dissipativity (storage fn λ≡0, α_ℓ(r)=ρc·r²).

N_B_noball = N;

con_B_noball.u_max              = u_max;
con_B_noball.y_min_active       = false;
con_B_noball.los_cone.active    = false;
con_B_noball.use_avg_power      = false;
con_B_noball.avg_power.P_avail  = P_avail;
con_B_noball.use_terminal_ball  = false;
con_B_noball.use_passive_safety = false;

params_B_noball   = params_base;
params_B_noball.N = N_B_noball;

fprintf(' Run B-noball — ρc=%.0e, N=%d, ball OFF\n', rho_c, N_B_noball);
resB_noball = run_closing(cost_B, con_B_noball, params_B_noball);

if isnan(resB_noball.conv_step)
    N_B_noball        = 40;
    params_B_noball.N = N_B_noball;
    fprintf('\n[Retry] N=20 did not converge; retrying with N=%d...\n\n', N_B_noball);
    fprintf('══════════════════════════════════════════════\n');
    fprintf(' Run B-noball (retry) — N=%d, ball OFF\n', N_B_noball);
    fprintf('══════════════════════════════════════════════\n');
    resB_noball = run_closing(cost_B, con_B_noball, params_B_noball);
end
fprintf('B-noball: used N=%d, converged=%s\n', N_B_noball, ...
    mat2str(~isnan(resB_noball.conv_step)));

%% Run B-ball: strict dissipativity + terminal ball ON 
% terminal constraint never binds.  Active-fraction ≈ 0 confirms this.

r_switch_ball_B = 200;   

con_B_ball = con_B_noball;
con_B_ball.use_terminal_ball      = true;
con_B_ball.terminal_ball.r_switch = r_switch_ball_B;

params_B_ball = params_base;   % N=20 

fprintf(' Run B-ball — ρc=%.0e, N=20, ball ON r=%d m\n', rho_c, r_switch_ball_B);
resB_ball = run_closing(cost_B, con_B_ball, params_B_ball);

%% Run C: pure energy + CLF terminal cost, no ball 
% Terminal cost Vf = e'P_clf e (P_clf from DARE) is a valid local CLF for
% the error dynamics.

con_C.u_max              = u_max;
con_C.y_min_active       = false;
con_C.los_cone.active    = false;
con_C.use_avg_power      = true;
con_C.avg_power.P_avail  = P_avail;
con_C.use_terminal_ball  = false;
con_C.use_passive_safety = false;

params_C = params_base;


fprintf(' Run C — energy ℓ=‖u‖², CLF terminal cost V_f=e''P_clf e\n');
resC = run_closing(cost_C, con_C, params_C);

%% Run MPC: Standard regulation MPC (comparison baseline)

cost_std = struct('Q', Q_dare, 'R', R_dare, 'P', P_clf);

con_MPC.u_max              = u_max;
con_MPC.y_min_active       = false;
con_MPC.los_cone.active    = false;
con_MPC.use_avg_power      = true;
con_MPC.avg_power.P_avail  = P_avail;
con_MPC.use_terminal_ball  = false;
con_MPC.use_passive_safety = false;

params_MPC = params_base;

fprintf(' Run MPC — standard regulation ℓ=e''Qe+u''Ru  (Q,R = DARE weights, P = LQR cost-to-go)\n');
resMPC = run_closing(cost_std, con_MPC, params_MPC);


%% Postprocessing

% Ball-role evidence (trajectory comparison) 

n_B_cmp  = min(resB_noball.n_ctrl, resB_ball.n_ctrl);
dX_B_max = max(abs(resB_noball.X(1:3, 1:n_B_cmp+1) - resB_ball.X(1:3, 1:n_B_cmp+1)), [], 'all');

fprintf('\nBall-role evidence:\n');
fprintf('  A2 total Δv (no ball, Q=P=0):     %.2e m/s  → ball is load-bearing in A1\n', resA2.dv);
fprintf('  ||pos(B-noball) - pos(B-ball)||:  %.2e m   → ball is slack in B-ball\n\n', dX_B_max);

% MPC comparison — reuse the first-class Run MPC (resMPC) computed above
X_mpc  = resMPC.X;
U_mpc  = resMPC.U;
t_mpc  = resMPC.t_vec;
dv_mpc = resMPC.dv;

fprintf('\n── Δv summary ───────────────────────────────────────────────\n');
fprintf('  %-32s  %10s\n', 'Run', 'Δv [m/s]');
fprintf('  %-32s  %10.4f\n', 'A1 — energy, ball ON',               resA1.dv);
fprintf('  %-32s  %10.4f\n', sprintf('B-noball (N=%d)',N_B_noball),    resB_noball.dv);
fprintf('  %-32s  %10.4f\n', 'B-ball  — strict dissip.',              resB_ball.dv);
fprintf('  %-32s  %10.4f\n', 'C  — energy + CLF terminal cost',       resC.dv);
fprintf('  %-32s  %10.4f\n', 'MPC — standard regulation',             dv_mpc);
fprintf('  A1 vs MPC savings: %.4f m/s  (%.1f%%)\n', ...
    dv_mpc - resA1.dv, 100*(dv_mpc - resA1.dv)/dv_mpc);
fprintf('  C  vs MPC savings: %.4f m/s  (%.1f%%)\n', ...
    dv_mpc - resC.dv,  100*(dv_mpc - resC.dv) /dv_mpc);
fprintf('─────────────────────────────────────────────────────────────\n\n');

% Run A1-ps: energy + ball ON + passive safety SCA 
% SCA per-step constraint forces free-drift safety at every predicted step.

con_A1ps                             = con_A1;
con_A1ps.use_passive_safety          = true;
con_A1ps.passive_safety.r_KOS        = r_KOS;
con_A1ps.passive_safety.safety_margin = ps_margin;
con_A1ps.passive_safety.x_switch     = x_switch;

fprintf(' Run A1-ps — energy, ball ON, passive safety SCA ON\n');
resA1_ps = run_closing(cost_A, con_A1ps, params_A1);
fprintf('A1-ps: passive_safe_fraction = %.1f%%\n\n', 100*resA1_ps.ps_frac);

% Horizon / ρc sweep (B-noball, 60 steps)
% Visualises σ(N) ≈ √(δ/(ρc·N)) tradeoff: smaller ρc needs larger N to
% achieve ball-free convergence within 60 steps.

N_sweep     = [10, 20, 40];
rho_c_sweep = [1e-12, 1e-11, 1e-10];
n_sw        = numel(N_sweep) * numel(rho_c_sweep);
sw_N        = zeros(1, n_sw);
sw_rho_c    = zeros(1, n_sw);
sw_conv     = false(1, n_sw);
sw_steps    = NaN(1, n_sw);
sw_dv       = zeros(1, n_sw);
sw_td       = zeros(1, n_sw);

fprintf('Horizon/ρc sweep  (B-noball config, n_steps_override=60)\n\n');

idx = 0;
for nn = N_sweep
    for rc = rho_c_sweep
        idx         = idx + 1;
        cost_sw     = struct('Q', rc*eye(6), 'R', eye(3), 'P', zeros(6));
        params_sw   = params_base;
        params_sw.N                = nn;
        params_sw.n_steps_override = 60;
        params_sw.verbose          = false;

        fprintf('── sweep %d/%d  N=%2d  ρc=%.0e  ...', idx, n_sw, nn, rc);
        r_sw = run_closing(cost_sw, con_B_noball, params_sw);

        sw_N(idx)    = nn;
        sw_rho_c(idx) = rc;
        sw_conv(idx)  = ~isnan(r_sw.conv_step);
        sw_steps(idx) = r_sw.conv_step;
        sw_dv(idx)    = r_sw.dv;
        sw_td(idx)    = r_sw.term_dist;
        fprintf('  converged=%s  term_dist=%.0f m\n', mat2str(sw_conv(idx)), sw_td(idx));
    end
end

fprintf('\n── Horizon/ρc sweep results ─────────────────────────────────────────────\n');
fprintf('  %4s  %8s  %9s  %6s  %10s  %12s\n', ...
    'N', 'rho_c', 'converged', 'steps', 'dv [m/s]', 'term_dist [m]');
for ii = 1:n_sw
    cstr = 'no';  sstr = '—';
    if sw_conv(ii), cstr = 'yes'; end
    if ~isnan(sw_steps(ii)), sstr = num2str(sw_steps(ii)); end
    fprintf('  %4d  %8.0e  %9s  %6s  %10.4f  %12.2f\n', ...
        sw_N(ii), sw_rho_c(ii), cstr, sstr, sw_dv(ii), sw_td(ii));
end
fprintf('──────────────────────────────────────────────────────────────────────────\n\n');

% Precompute free-drift STMs for passive-safety sweeps
% Same CWH system at dt=120s, shared across all sweep sections below.
N_safe_sw  = 554;
Ad_pow_sw  = zeros(6, 6, N_safe_sw);
Atmp_sw    = Ad_clf;
for jj_sw  = 1:N_safe_sw
    Ad_pow_sw(:,:,jj_sw) = Atmp_sw;
    Atmp_sw = Atmp_sw * Ad_clf;
end

% Sweep: passive safety vs r_KOS (existing trajectories)

r_KOS_vec  = [10, 25, 50, 75, 100, 150, 200, 300];
ps_runs    = {resA1, resC, resB_noball};
ps_labels  = {'A1', 'C', 'B-noball'};
n_ps_runs  = numel(ps_runs);
ps_vs_rkos = zeros(n_ps_runs, numel(r_KOS_vec));

for ri = 1:n_ps_runs
    X_ri = ps_runs{ri}.X;
    n_ri = size(X_ri, 2);
    for ki = 1:numel(r_KOS_vec)
        is_s = check_passive_safety(X_ri, n_ri, Ad_pow_sw, r_KOS_vec(ki), N_safe_sw);
        ps_vs_rkos(ri, ki) = mean(is_s);
    end
end

% Sweep: passive safety vs N for Run C (re-runs at r_KOS=50m)
% Show how prediction horizon affects the trajectory's inherent passive safety.
N_ps_vec = [5, 10, 15, 20, 30, 40];
ps_vs_N  = zeros(1, numel(N_ps_vec));
fprintf('\nPassive safety vs N sweep (Run C config, r_KOS=%d m)...\n', r_KOS);
for ni = 1:numel(N_ps_vec)
    params_psN         = params_base;
    params_psN.N       = N_ps_vec(ni);
    params_psN.verbose = false;
    r_psN              = run_closing(cost_C, con_C, params_psN);
    ps_vs_N(ni)        = r_psN.ps_frac;
    fprintf('  N=%2d: ps_frac=%.1f%%  (converged=%s)\n', ...
        N_ps_vec(ni), 100*r_psN.ps_frac, mat2str(~isnan(r_psN.conv_step)));
end

% Sweep: shrinking terminal ball → turnpike behavior

r_ball_vec = [500, 300, 200, 150, 100, 50];
n_ball     = numel(r_ball_vec);
ball_res   = cell(1, n_ball);

fprintf('\nShrinking ball sweep (A1 config, N=20)...\n');
for bi = 1:n_ball
    con_bi                            = con_A1;
    con_bi.terminal_ball.r_switch     = r_ball_vec(bi);
    params_bi                         = params_base;
    params_bi.verbose                 = false;
    fprintf('  r_ball=%3d m  ...  ', r_ball_vec(bi));
    r_bi = run_closing(cost_A, con_bi, params_bi);
    ball_res{bi} = r_bi;
    conv_str = 'DNF';
    if ~isnan(r_bi.conv_step), conv_str = sprintf('step %d', r_bi.conv_step); end
    fprintf('dv=%.3f m/s  conv=%s  tp_frac=%.2f\n', ...
        r_bi.dv, conv_str, r_bi.turnpike_fraction);
end

%% Figures
c_A1   = [0.85 0.33 0.10];  
c_A2   = [0.55 0.55 0.55];   
c_Bnb  = [0.20 0.63 0.17];   
c_Bbal = [0.58 0.40 0.74];  
c_C    = [0.00 0.68 0.68];   
c_mpc  = [0.12 0.47 0.71];   

t_A1_h      = resA1.t_vec       / 3600;
t_Bnb_h     = resB_noball.t_vec / 3600;
t_Bbal_h    = resB_ball.t_vec   / 3600;
t_C_h       = resC.t_vec        / 3600;
n_ctrl_A1   = resA1.n_ctrl;
n_ctrl_Bnb  = resB_noball.n_ctrl;
n_ctrl_Bbal = resB_ball.n_ctrl;
n_ctrl_C    = resC.n_ctrl;

% Figure 1 — Hill-frame: A1 vs A2
fig1 = figure('Name', 'S4 — Hill Frame: A1 (ball ON) vs A2 (ball OFF)');
ax1  = axes(fig1);
hold(ax1,'on'); grid(ax1,'on'); axis(ax1,'equal');
plot(ax1, resA1.X(2,:), resA1.X(1,:), '-', 'Color', c_A1, ...
    'LineWidth', 1.4, 'DisplayName', 'A1 — ball ON (converges)');
plot(ax1, resA2.X(2,1), resA2.X(1,1), 'x', 'Color', c_A2, ...
    'MarkerSize', 14, 'LineWidth', 2.5, ...
    'DisplayName', 'A2 — ball OFF  ($u\equiv 0$, stuck at IC)');
plot(ax1, x0(2), x0(1), 'o', 'Color', [0.4 0.4 0.4], ...
    'MarkerFaceColor', [0.4 0.4 0.4], 'MarkerSize', 6, 'HandleVisibility', 'off');
plot(ax1, x_switch(2), x_switch(1), 's', 'Color', [0.2 0.6 0.2], ...
    'MarkerFaceColor', [0.2 0.6 0.2], 'MarkerSize', 7, 'DisplayName', 'Switch point');
plot(ax1, 0, 0, '+k', 'MarkerSize', 9, 'LineWidth', 1.5, 'DisplayName', 'Target');
xlabel(ax1, 'Along-track $y$ [m]', 'Interpreter', 'latex');
ylabel(ax1, 'Radial $x$ [m]',      'Interpreter', 'latex');
title(ax1, 'Hill-Frame Trajectory — A1 (ball ON) vs A2 (ball OFF)', 'Interpreter', 'latex');
legend(ax1, 'Location', 'northeast', 'Interpreter', 'latex');
if exist('matlab2tikz', 'file')
    matlab2tikz('results/figures/scenario4/AB_hill_frame.tikz', ...
        'figurehandle', fig1, 'showInfo', false);
end

% Figure 2 — Cumulative Δv: A1 vs C vs MPC
fig2 = figure('Name', 'S4 — Cumulative \Deltav: A1 vs C vs MPC');
ax2  = axes(fig2);
hold(ax2,'on'); grid(ax2,'on');
plot(ax2, t_A1_h(1:n_ctrl_A1), ...
    cumsum(vecnorm(resA1.U(:,1:n_ctrl_A1), 2, 1))*dt, ...
    'Color', c_A1,  'LineWidth', 1.4, 'DisplayName', 'A1 — EMPC, ball ON');
plot(ax2, t_C_h(1:n_ctrl_C), ...
    cumsum(vecnorm(resC.U(:,1:n_ctrl_C), 2, 1))*dt, ...
    'Color', c_C,   'LineWidth', 1.4, 'DisplayName', 'C — EMPC, CLF $V_f$');
plot(ax2, t_mpc(1:end-1)/3600, ...
    cumsum(vecnorm(U_mpc(:,1:end-1), 2, 1))*dt, ...
    'Color', c_mpc, 'LineWidth', 1.4, 'DisplayName', 'MPC — standard regulation');
xlabel(ax2, 'Time [h]');
ylabel(ax2, 'Cumulative $\Delta v$ [m/s]', 'Interpreter', 'latex');
title(ax2, 'Cumulative $\Delta v$ — A1 vs C vs Standard MPC', 'Interpreter', 'latex');
legend(ax2, 'Location', 'northwest', 'Interpreter', 'latex');
if exist('matlab2tikz', 'file')
    matlab2tikz('results/figures/scenario4/AB_cumulative_dv.tikz', ...
        'figurehandle', fig2, 'showInfo', false);
end

% Figure 3 — Thrust profile A1 (log scale)
% Plateau for A1 is an energy-cost cruise: near-constant thrust driven by the
% terminal ball with no state-cost attraction to any point (Q=0).
tp_label_A1 = 'Energy-cost cruise';

fig3 = figure('Name', 'S4 — A1 Thrust Profile');
ax3  = axes(fig3);
hold(ax3,'on'); grid(ax3,'on');
set(ax3, 'YScale', 'log');
k_vec  = 1:n_ctrl_A1;
u_plot = max(resA1.u_norm_seq, 1e-20);
stem(ax3, k_vec, u_plot, 'Color', c_A1, 'LineWidth', 0.8, ...
    'MarkerSize', 3, 'DisplayName', '$\|u_k\|$');
if ~isempty(resA1.turnpike_steps)
    yl = ylim(ax3);
    if yl(1) <= 0 || ~isfinite(log10(yl(1)))
        pv = u_plot(u_plot>1e-20); yl(1)=min(pv)*0.01; ylim(ax3,yl); yl=ylim(ax3);
    end
    tp_lo = resA1.turnpike_steps(1)-0.5;  tp_hi = resA1.turnpike_steps(end)+0.5;
    h_p = patch(ax3,[tp_lo tp_hi tp_hi tp_lo],[yl(1) yl(1) yl(2) yl(2)],[0.8 0.8 0.8], ...
        'FaceAlpha',0.35,'EdgeColor','none','HandleVisibility','off');
    uistack(h_p,'bottom');
    ym = 10^(0.5*(log10(yl(1))+log10(yl(2))));
    if resA1.turnpike_steps(1) > 2
        text(ax3, mean([0.5,tp_lo]), ym, 'Entry', ...
            'HorizontalAlignment','center','FontSize',8,'Color',[0.2 0.2 0.2]);
    end
    text(ax3, mean([tp_lo,tp_hi]), ym, tp_label_A1, ...
        'HorizontalAlignment','center','FontSize',8,'Color',[0.2 0.2 0.2]);
    if resA1.turnpike_steps(end) < n_ctrl_A1-1
        text(ax3, mean([tp_hi,n_ctrl_A1+0.5]), ym, 'Exit', ...
            'HorizontalAlignment','center','FontSize',8,'Color',[0.2 0.2 0.2]);
    end
end
xlabel(ax3, 'Step $k$',             'Interpreter', 'latex');
ylabel(ax3, '$\|u_k\|$ [m/s$^2$]', 'Interpreter', 'latex');
title(ax3, 'A1 — $V_f=0$ EMPC, pure energy, $N=20$, ball ON', 'Interpreter', 'latex');
legend(ax3, 'Interpreter', 'latex', 'Location', 'northeast');
if exist('matlab2tikz', 'file')
    matlab2tikz('results/figures/scenario4/A1_thrust_profile.tikz', ...
        'figurehandle', fig3, 'showInfo', false);
end

% Figure 4 — Hill-frame: A1, B-noball, C, MPC  (full view + switch-point zoom)
conv_tol_s4 = 50;         % [m] run_closing convergence threshold
zoom_hw     = 150;        % [m] half-width of the switch-point zoom window

fig4 = figure('Name', 'S4 — Hill Frame: A1, B-noball, C, MPC (full + zoom)');
fig4.Position(3:4) = [760, 300];

% Panel A: full trajectory (IC → switch)
ax4a = subplot(1, 2, 1);
hold(ax4a,'on'); grid(ax4a,'on'); axis(ax4a,'equal');
plot(ax4a, resA1.X(2,:),       resA1.X(1,:),       '-', 'Color', c_A1, 'LineWidth', 1.2, ...
    'DisplayName', 'A1 — energy, ball ON');
plot(ax4a, resB_noball.X(2,:), resB_noball.X(1,:), '-', 'Color', c_Bnb, 'LineWidth', 1.2, ...
    'DisplayName', sprintf('B-noball ($\\rho_c$=%.0e, N=%d)', rho_c, N_B_noball));
plot(ax4a, resC.X(2,:),        resC.X(1,:),        '-', 'Color', c_C,  'LineWidth', 1.2, ...
    'DisplayName', 'C — energy + CLF $V_f$');
plot(ax4a, resMPC.X(2,:),      resMPC.X(1,:),      '-', 'Color', c_mpc, 'LineWidth', 1.2, ...
    'DisplayName', 'MPC — standard regulation');
plot(ax4a, x0(2), x0(1), 'o', 'Color', [0.4 0.4 0.4], ...
    'MarkerFaceColor', [0.4 0.4 0.4], 'MarkerSize', 6, 'DisplayName', 'IC (10 km)');
plot(ax4a, x_switch(2), x_switch(1), 's', 'Color', [0.2 0.6 0.2], ...
    'MarkerFaceColor', [0.2 0.6 0.2], 'MarkerSize', 7, 'DisplayName', 'Switch point');
plot(ax4a, 0, 0, '+k', 'MarkerSize', 9, 'LineWidth', 1.5, 'DisplayName', 'Target');
% Dashed rectangle showing the zoom window
rectangle(ax4a, 'Position', [x_switch(2)-zoom_hw, x_switch(1)-zoom_hw, 2*zoom_hw, 2*zoom_hw], ...
    'EdgeColor', [0.4 0.4 0.4], 'LineStyle', '--', 'LineWidth', 0.8);
xlabel(ax4a, 'Along-track $y$ [m]', 'Interpreter', 'latex');
ylabel(ax4a, 'Radial $x$ [m]',      'Interpreter', 'latex');
title(ax4a, 'Full closing trajectory', 'Interpreter', 'latex');
legend(ax4a, 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 6);

% Panel B: zoom on the switch point
ax4b = subplot(1, 2, 2);
hold(ax4b,'on'); grid(ax4b,'on'); axis(ax4b,'equal');
plot(ax4b, resA1.X(2,:),       resA1.X(1,:),       '-', 'Color', c_A1, 'LineWidth', 1.2, ...
    'DisplayName', 'A1 — energy, ball ON');
plot(ax4b, resB_noball.X(2,:), resB_noball.X(1,:), '-', 'Color', c_Bnb, 'LineWidth', 1.2, ...
    'DisplayName', sprintf('B-noball ($\\rho_c$=%.0e, N=%d)', rho_c, N_B_noball));
plot(ax4b, resC.X(2,:),        resC.X(1,:),        '-', 'Color', c_C,  'LineWidth', 1.2, ...
    'DisplayName', 'C — energy + CLF $V_f$');
plot(ax4b, resMPC.X(2,:),      resMPC.X(1,:),      '-', 'Color', c_mpc, 'LineWidth', 1.2, ...
    'DisplayName', 'MPC — standard regulation');
theta_c = linspace(0, 2*pi, 200);
plot(ax4b, x_switch(2)+conv_tol_s4*cos(theta_c), x_switch(1)+conv_tol_s4*sin(theta_c), ...
    ':k', 'LineWidth', 0.9, 'DisplayName', sprintf('conv. tol %d m', conv_tol_s4));
plot(ax4b, x_switch(2), x_switch(1), 's', 'Color', [0.2 0.6 0.2], ...
    'MarkerFaceColor', [0.2 0.6 0.2], 'MarkerSize', 8, 'DisplayName', 'Switch point');
xlim(ax4b, x_switch(2) + [-zoom_hw, zoom_hw]);
ylim(ax4b, x_switch(1) + [-zoom_hw, zoom_hw]);
xlabel(ax4b, 'Along-track $y$ [m]', 'Interpreter', 'latex');
ylabel(ax4b, 'Radial $x$ [m]',      'Interpreter', 'latex');
title(ax4b, 'Zoom: switch-point arrival', 'Interpreter', 'latex');
legend(ax4b, 'Location', 'southeast', 'Interpreter', 'latex', 'FontSize', 6);

sgtitle(fig4, 'Hill-Frame Trajectory — A1, B-noball, C, MPC', 'Interpreter', 'latex');
if exist('matlab2tikz', 'file')
    matlab2tikz('results/figures/scenario4/B_hill_frame.tikz', ...
        'figurehandle', fig4, 'showInfo', false);
end

% Figure 5 — Horizon / ρc sweep: convergence steps + Δv (two panels)
N_vals  = unique(sw_N,     'sorted');
rc_vals = unique(sw_rho_c, 'sorted');
n_N  = length(N_vals);
n_rc = length(rc_vals);

% Build n_N × n_rc matrices from flat sweep vectors
step_mat = NaN(n_N, n_rc);
dv_mat   = NaN(n_N, n_rc);
for ii = 1:n_sw
    ni = find(N_vals  == sw_N(ii));
    ri = find(rc_vals == sw_rho_c(ii));
    if sw_conv(ii), step_mat(ni, ri) = sw_steps(ii); end
    dv_mat(ni, ri) = sw_dv(ii);   % show Δv even for DNF (reveals wasted fuel)
end

dnf_cap  = 65;
disp_mat = step_mat;
disp_mat(isnan(disp_mat)) = dnf_cap;

c_rc3     = [0.85 0.33 0.10; 0.20 0.63 0.17; 0.12 0.47 0.71];
bar_grp_w = 0.8;   % MATLAB default total group width
xtick_lbl = arrayfun(@(x) sprintf('N=%d', x), N_vals, 'UniformOutput', false);

fig5 = figure('Name', 'S4 — B-noball Convergence vs Horizon/rho_c');
fig5.Position(3:4) = [900, 380];

% Panel A: steps to convergence
ax5a = subplot(1, 2, 1);
hold(ax5a, 'on'); grid(ax5a, 'on');
hb_a = bar(ax5a, disp_mat, 'grouped');
for ri = 1:n_rc
    hb_a(ri).FaceColor   = c_rc3(ri, :);
    hb_a(ri).DisplayName = sprintf('$\\rho_c = %.0e$', rc_vals(ri));
end
for ni = 1:n_N
    for ri = 1:n_rc
        if isnan(step_mat(ni, ri))
            xc = ni + (ri - (n_rc + 1) / 2) * (bar_grp_w / n_rc);
            text(ax5a, xc, dnf_cap + 1.5, 'DNF', ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
                'FontSize', 7, 'Color', [0.6 0.1 0.1]);
        end
    end
end
yline(ax5a, 60, '--k', 'LineWidth', 0.8, 'HandleVisibility', 'off');
set(ax5a, 'XTick', 1:n_N, 'XTickLabel', xtick_lbl);
xlabel(ax5a, 'Prediction horizon $N$', 'Interpreter', 'latex');
ylabel(ax5a, 'Steps to convergence  (budget: 60)');
title(ax5a, 'Convergence steps', 'Interpreter', 'latex');
legend(ax5a, 'Interpreter', 'latex', 'Location', 'northeast');
ylim(ax5a, [0, dnf_cap + 8]);

% Panel B: total Δv
ax5b = subplot(1, 2, 2);
hold(ax5b, 'on'); grid(ax5b, 'on');
hb_b = bar(ax5b, dv_mat, 'grouped');
for ri = 1:n_rc
    hb_b(ri).FaceColor        = c_rc3(ri, :);
    hb_b(ri).HandleVisibility = 'off';   % legend already in Panel A
end
for ni = 1:n_N
    for ri = 1:n_rc
        if isnan(step_mat(ni, ri))
            xc = ni + (ri - (n_rc + 1) / 2) * (bar_grp_w / n_rc);
            text(ax5b, xc, dv_mat(ni, ri) + 0.5, '(DNF)', ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
                'FontSize', 6, 'Color', [0.6 0.1 0.1]);
        end
    end
end
set(ax5b, 'XTick', 1:n_N, 'XTickLabel', xtick_lbl);
xlabel(ax5b, 'Prediction horizon $N$', 'Interpreter', 'latex');
ylabel(ax5b, 'Total $\Delta v$ [m/s]', 'Interpreter', 'latex');
title(ax5b, 'Total $\Delta v$', 'Interpreter', 'latex');

sgtitle(fig5, 'B-noball: convergence vs $(N,\,\rho_c)$ — strict dissipativity, no ball', ...
    'Interpreter', 'latex');

if exist('matlab2tikz', 'file')
    matlab2tikz('results/figures/scenario4/horizon_sweep.tikz', ...
        'figurehandle', fig5, 'showInfo', false);
end

if exist('matlab2tikz', 'file')
    fprintf('TikZ figures exported to results/figures/scenario4/\n\n');
else
    fprintf('TikZ export skipped (matlab2tikz not on path — run setup.m first).\n\n');
end

% Figure 6 — A1 vs A1-ps: Hill-frame safety coloring + thrust profiles
n_ctrl_A1ps  = resA1_ps.n_ctrl;
c_A1ps       = [0.12 0.47 0.71];   % blue  — A1-ps
c_safe_dot   = [0.18 0.55 0.34];   % green — passively safe states
c_unsafe_dot = [0.85 0.20 0.20];   % red   — passively unsafe states

fig6 = figure('Name', 'S4 — A1 vs A1-ps');
fig6.Position(3:4) = [900, 380];

% Panel A: Hill frame, A1 states safety-coloured + A1-ps trajectory
ax6a = subplot(1, 2, 1);
hold(ax6a, 'on'); grid(ax6a, 'on'); axis(ax6a, 'equal');

safe_msk   = resA1.is_safe;
unsafe_msk = ~resA1.is_safe;
scatter(ax6a, resA1.X(2,  safe_msk),  resA1.X(1,  safe_msk),  14, ...
    c_safe_dot,   'filled', 'DisplayName', 'A1 — passively safe');
scatter(ax6a, resA1.X(2, unsafe_msk), resA1.X(1, unsafe_msk), 40, ...
    c_unsafe_dot, 'filled', '^', 'DisplayName', 'A1 — passively unsafe');
plot(ax6a, resA1_ps.X(2,:), resA1_ps.X(1,:), '-', 'Color', c_A1ps, ...
    'LineWidth', 1.4, 'DisplayName', 'A1-ps — SCA ON (all safe)');
plot(ax6a, x_switch(2), x_switch(1), 's', 'Color', [0.2 0.6 0.2], ...
    'MarkerFaceColor', [0.2 0.6 0.2], 'MarkerSize', 7, 'DisplayName', 'Switch point');
plot(ax6a, 0, 0, '+k', 'MarkerSize', 9, 'LineWidth', 1.5, 'DisplayName', 'Target');
plot(ax6a, x0(2), x0(1), 'o', 'Color', [0.4 0.4 0.4], ...
    'MarkerFaceColor', [0.4 0.4 0.4], 'MarkerSize', 5, 'HandleVisibility', 'off');
xlabel(ax6a, 'Along-track $y$ [m]', 'Interpreter', 'latex');
ylabel(ax6a, 'Radial $x$ [m]',      'Interpreter', 'latex');
title(ax6a, 'Hill-frame trajectory', 'Interpreter', 'latex');
legend(ax6a, 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 8);

% Panel B: Thrust profiles (log scale)
ax6b = subplot(1, 2, 2);
hold(ax6b, 'on'); grid(ax6b, 'on');
set(ax6b, 'YScale', 'log');

stem(ax6b, 1:n_ctrl_A1, max(resA1.u_norm_seq, 1e-20), ...
    'Color', c_A1, 'LineWidth', 0.7, 'MarkerSize', 3, ...
    'DisplayName', sprintf('A1  (step %d,  %.2f m/s,  ps=%.0f%%)', ...
        resA1.conv_step, resA1.dv, 100*resA1.ps_frac));
stem(ax6b, 1:n_ctrl_A1ps, max(resA1_ps.u_norm_seq, 1e-20), ...
    'Color', c_A1ps, 'LineWidth', 0.7, 'MarkerSize', 3, ...
    'DisplayName', sprintf('A1-ps  (step %d,  %.2f m/s,  ps=%.0f%%)', ...
        resA1_ps.conv_step, resA1_ps.dv, 100*resA1_ps.ps_frac));
xlabel(ax6b, 'Step $k$',             'Interpreter', 'latex');
ylabel(ax6b, '$\|u_k\|$ [m/s$^2$]', 'Interpreter', 'latex');
title(ax6b, 'Thrust profile (log scale)', 'Interpreter', 'latex');
legend(ax6b, 'Interpreter', 'latex', 'Location', 'northeast', 'FontSize', 8);

sgtitle(fig6, 'A1 vs A1-ps — effect of SCA passive-safety constraints', ...
    'Interpreter', 'latex');
if exist('matlab2tikz', 'file')
    matlab2tikz('results/figures/scenario4/A1_vs_A1ps.tikz', ...
        'figurehandle', fig6, 'showInfo', false);
end

% Figure 8 — Passive safety sensitivity (2 panels)
% Panel A: ps_frac vs r_KOS for A1/C/B-noball at N=20.
% Panel B: ps_frac vs N for Run C at r_KOS=50m.

c_ps = {c_A1, c_C, c_Bnb};   % colours match run colours

fig8 = figure('Name', 'S4 — Passive Safety Sensitivity');
fig8.Position(3:4) = [900, 380];

ax8a = subplot(1, 2, 1);
hold(ax8a, 'on'); grid(ax8a, 'on');
for ri = 1:n_ps_runs
    plot(ax8a, r_KOS_vec, 100*ps_vs_rkos(ri,:), '-o', ...
        'Color', c_ps{ri}, 'LineWidth', 1.4, 'MarkerSize', 5, ...
        'DisplayName', ps_labels{ri});
end
xline(ax8a, r_KOS, '--k', 'LineWidth', 0.8, 'HandleVisibility', 'off');
xlabel(ax8a, 'KOS radius $r_\mathrm{KOS}$ [m]', 'Interpreter', 'latex');
ylabel(ax8a, 'Passive-safe fraction [\%]');
title(ax8a, 'ps\_frac vs $r_\mathrm{KOS}$  ($N=20$)', 'Interpreter', 'latex');
legend(ax8a, 'Location', 'southwest');
ylim(ax8a, [-5, 105]);

ax8b = subplot(1, 2, 2);
hold(ax8b, 'on'); grid(ax8b, 'on');
plot(ax8b, N_ps_vec, 100*ps_vs_N, '-s', 'Color', c_C, 'LineWidth', 1.4, 'MarkerSize', 6);
xline(ax8b, N, '--k', 'LineWidth', 0.8);
xlabel(ax8b, 'Prediction horizon $N$', 'Interpreter', 'latex');
ylabel(ax8b, 'Passive-safe fraction [\%]');
title(ax8b, 'Run C: ps\_frac vs $N$  ($r_\mathrm{KOS}=50$ m)', 'Interpreter', 'latex');
ylim(ax8b, [-5, 105]);

sgtitle(fig8, 'Passive safety sensitivity — $r_\mathrm{KOS}$ and horizon $N$', ...
    'Interpreter', 'latex');
if exist('matlab2tikz', 'file')
    matlab2tikz('results/figures/scenario4/ps_sensitivity.tikz', ...
        'figurehandle', fig8, 'showInfo', false);
end

% Figure 9 — Shrinking ball: thrust profiles + summary metrics (2 panels)

% Colour ramp from orange (large ball) to dark red (small ball)
n_b    = n_ball;
c_ball = zeros(n_b, 3);
for bi = 1:n_b
    t_bi       = (bi - 1) / max(n_b - 1, 1);  % 0 → 1 (large → small)
    c_ball(bi,:) = (1 - t_bi)*[0.95 0.60 0.20] + t_bi*[0.55 0.05 0.05];
end

fig9 = figure('Name', 'S4 — Shrinking Ball: Thrust & Metrics');
fig9.Position(3:4) = [900, 400];

ax9a = subplot(1, 2, 1);
hold(ax9a, 'on'); grid(ax9a, 'on');
set(ax9a, 'YScale', 'log');
for bi = 1:n_b
    r_bi   = ball_res{bi};
    nc_bi  = r_bi.n_ctrl;
    if nc_bi < 1, continue; end
    u_bi   = max(r_bi.u_norm_seq, 1e-20);
    plot(ax9a, 1:nc_bi, u_bi, '-', 'Color', c_ball(bi,:), 'LineWidth', 1.0, ...
        'DisplayName', sprintf('$r_\\mathrm{ball}=%d$ m', r_ball_vec(bi)));
end
xlabel(ax9a, 'Step $k$',             'Interpreter', 'latex');
ylabel(ax9a, '$\|u_k\|$ [m/s$^2$]', 'Interpreter', 'latex');
title(ax9a, 'Thrust profiles (log scale)', 'Interpreter', 'latex');
legend(ax9a, 'Interpreter', 'latex', 'Location', 'northeast', 'FontSize', 7);

ax9b = subplot(1, 2, 2);
hold(ax9b, 'on'); grid(ax9b, 'on');
ball_dv        = cellfun(@(r) r.dv,             ball_res);
ball_conv      = cellfun(@(r) r.conv_step,      ball_res);
ball_tp_frac   = cellfun(@(r) r.turnpike_fraction, ball_res);
ball_conv(isnan(ball_conv)) = NaN;

yyaxis(ax9b, 'left');
plot(ax9b, r_ball_vec, ball_dv, '-o', 'Color', c_A1, 'LineWidth', 1.4, ...
    'MarkerSize', 6, 'DisplayName', '$\Delta v$ [m/s]');
ylabel(ax9b, 'Total $\Delta v$ [m/s]', 'Interpreter', 'latex');

yyaxis(ax9b, 'right');
plot(ax9b, r_ball_vec, ball_tp_frac, '-s', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.4, ...
    'MarkerSize', 6, 'DisplayName', 'Turnpike fraction');
ylabel(ax9b, 'Turnpike fraction [-]');

xlabel(ax9b, 'Ball radius $r_\mathrm{ball}$ [m]', 'Interpreter', 'latex');
title(ax9b, '$\Delta v$ and turnpike fraction vs $r_\mathrm{ball}$', 'Interpreter', 'latex');
legend(ax9b, 'Interpreter', 'latex', 'Location', 'best', 'FontSize', 8);

sgtitle(fig9, 'Shrinking terminal ball — effect on thrust profile and turnpike', ...
    'Interpreter', 'latex');
if exist('matlab2tikz', 'file')
    matlab2tikz('results/figures/scenario4/shrinking_ball.tikz', ...
        'figurehandle', fig9, 'showInfo', false);
end

% Figure 10 — Terminal ball as an actuator-authority knob (table + focused plot)
% Turnpike framing dropped: characterise how shrinking the ball trades
% convergence time against peak control effort, with a feasibility floor
% (peak thrust rises toward u_max as the ball tightens, then the QP goes
% infeasible / DNF).
out_dir_s4 = 'results/figures/scenario4';
[~, ~]     = mkdir(out_dir_s4);

ball_dv      = zeros(1, n_ball);
ball_peaku   = zeros(1, n_ball);
ball_satfrac = zeros(1, n_ball);
ball_tconv_h = nan(1, n_ball);
ball_ok      = false(1, n_ball);

fprintf('\n── Terminal-ball sweep: control demand vs radius ──────────────\n');
fprintf('  %8s  %6s  %9s  %9s  %10s  %6s\n', ...
    'r_ball', 'conv', 't_conv[h]', 'dv[m/s]', 'peak|u|', 'sat%');
for bi = 1:n_ball
    r = ball_res{bi};
    ball_dv(bi)    = r.dv;
    ball_peaku(bi) = max([r.u_norm_seq, 0]);
    if r.n_ctrl >= 1
        sat = any(abs(r.U(:, 1:r.n_ctrl)) >= 0.999*u_max, 1);   % per-axis saturation
        ball_satfrac(bi) = mean(sat);
    end
    ball_ok(bi) = ~isnan(r.conv_step);
    if ball_ok(bi), ball_tconv_h(bi) = r.conv_step * dt / 3600; end
    conv_str = 'DNF'; if ball_ok(bi), conv_str = sprintf('%d', r.conv_step); end
    fprintf('  %6d m  %6s  %9.2f  %9.4f  %10.2e  %5.0f%%\n', ...
        r_ball_vec(bi), conv_str, ball_tconv_h(bi), ball_dv(bi), ...
        ball_peaku(bi), 100*ball_satfrac(bi));
end
fprintf('───────────────────────────────────────────────────────────────\n');

% LaTeX table for the report (plain yes/no — no pifont dependency)
fid_bt = fopen(fullfile(out_dir_s4, 'ball_sweep_table.tex'), 'w');
fprintf(fid_bt, '%% Auto-generated by sim_scenario4_closing.m\n');
fprintf(fid_bt, '\\begin{tabular}{r c r r r r}\n\\toprule\n');
fprintf(fid_bt, ['$r_{\\mathrm{ball}}$ [m] & converged & $t_{\\mathrm{conv}}$ [h] & ', ...
    '$\\Delta v$ [m/s] & peak $\\|u\\|$ [m/s$^2$] & sat.\\ [\\%%] \\\\\n\\midrule\n']);
for bi = 1:n_ball
    conv_tex  = 'no';  if ball_ok(bi), conv_tex  = 'yes'; end
    tconv_tex = '--';  if ball_ok(bi), tconv_tex = sprintf('%.2f', ball_tconv_h(bi)); end
    fprintf(fid_bt, '%d & %s & %s & %.4f & %.2e & %.0f \\\\\n', ...
        r_ball_vec(bi), conv_tex, tconv_tex, ball_dv(bi), ball_peaku(bi), 100*ball_satfrac(bi));
end
fprintf(fid_bt, '\\bottomrule\n\\end{tabular}\n');
fclose(fid_bt);

% Focused figure: convergence time vs control demand vs ball radius
fig10 = figure('Name', 'S4 — Terminal ball: control demand vs radius');
ax10  = axes(fig10); hold(ax10, 'on'); grid(ax10, 'on');
set(ax10, 'XDir', 'reverse');   % shrinking ball reads left→right

yyaxis(ax10, 'left');
plot(ax10, r_ball_vec(ball_ok), ball_tconv_h(ball_ok), '-o', ...
    'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'Convergence time');
ylabel(ax10, 'Convergence time [h]', 'Interpreter', 'latex');

yyaxis(ax10, 'right');
plot(ax10, r_ball_vec, ball_peaku, '-s', ...
    'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'Peak $\|u\|$');
yline(ax10, u_max, ':', '$u_{\max}$ (per axis)', ...
    'Interpreter', 'latex', 'HandleVisibility', 'off');
ylabel(ax10, 'Peak $\|u\|$ [m/s$^2$]', 'Interpreter', 'latex');

if any(~ball_ok)
    yl10 = ylim(ax10);
    plot(ax10, r_ball_vec(~ball_ok), yl10(1)*ones(1, sum(~ball_ok)), 'x', ...
        'Color', [0.8 0 0], 'MarkerSize', 10, 'LineWidth', 1.5, ...
        'DisplayName', 'infeasible (DNF)');
end

xlabel(ax10, 'Terminal ball radius $r_{\mathrm{ball}}$ [m]', 'Interpreter', 'latex');
title(ax10, 'Terminal ball: convergence speed vs control demand', 'Interpreter', 'latex');
legend(ax10, 'Interpreter', 'latex', 'Location', 'north');

if exist('matlab2tikz', 'file')
    matlab2tikz(fullfile(out_dir_s4, 'ball_control_demand.tikz'), ...
        'figurehandle', fig10, 'showInfo', false);
end

% Figure 7 — Run C thrust profile (log scale)
tp_label_C = 'CLF-guided cruise';

fig7 = figure('Name', 'S4 — C Thrust Profile');
ax7  = axes(fig7);
hold(ax7,'on'); grid(ax7,'on');
set(ax7, 'YScale', 'log');
k_vec_C  = 1:n_ctrl_C;
u_plot_C = max(resC.u_norm_seq, 1e-20);
stem(ax7, k_vec_C, u_plot_C, 'Color', c_C, 'LineWidth', 0.8, ...
    'MarkerSize', 3, 'DisplayName', '$\|u_k\|$');
if ~isempty(resC.turnpike_steps)
    yl7 = ylim(ax7);
    if yl7(1) <= 0 || ~isfinite(log10(yl7(1)))
        pv7 = u_plot_C(u_plot_C>1e-20);
        if ~isempty(pv7), yl7(1)=min(pv7)*0.01; ylim(ax7,yl7); yl7=ylim(ax7); end
    end
    tp_lo7 = resC.turnpike_steps(1)-0.5;  tp_hi7 = resC.turnpike_steps(end)+0.5;
    h_p7 = patch(ax7,[tp_lo7 tp_hi7 tp_hi7 tp_lo7],[yl7(1) yl7(1) yl7(2) yl7(2)], ...
        [0.8 0.8 0.8],'FaceAlpha',0.35,'EdgeColor','none','HandleVisibility','off');
    uistack(h_p7,'bottom');
    ym7 = 10^(0.5*(log10(yl7(1))+log10(yl7(2))));
    text(ax7, mean([tp_lo7,tp_hi7]), ym7, tp_label_C, ...
        'HorizontalAlignment','center','FontSize',8,'Color',[0.2 0.2 0.2]);
end
xlabel(ax7, 'Step $k$',             'Interpreter', 'latex');
ylabel(ax7, '$\|u_k\|$ [m/s$^2$]', 'Interpreter', 'latex');
title(ax7, 'C — $V_f=e^\top P_{\mathrm{clf}}e$ EMPC, pure energy, $N=20$, no ball', ...
    'Interpreter', 'latex');
legend(ax7, 'Interpreter', 'latex', 'Location', 'northeast');
if exist('matlab2tikz', 'file')
    matlab2tikz('results/figures/scenario4/C_thrust_profile.tikz', ...
        'figurehandle', fig7, 'showInfo', false);
end

%% PNG export
[~, ~] = mkdir('results/figures/scenario4');
print(fig1, 'results/figures/scenario4/AB_hill_frame',    '-dpng', '-r150');
print(fig2, 'results/figures/scenario4/AB_cumulative_dv', '-dpng', '-r150');
print(fig3, 'results/figures/scenario4/A1_thrust_profile','-dpng', '-r150');
print(fig4, 'results/figures/scenario4/B_hill_frame',     '-dpng', '-r150');
print(fig5, 'results/figures/scenario4/horizon_sweep',    '-dpng', '-r150');
print(fig6, 'results/figures/scenario4/A1_vs_A1ps',       '-dpng', '-r150');
print(fig7, 'results/figures/scenario4/C_thrust_profile', '-dpng', '-r150');
print(fig8, 'results/figures/scenario4/ps_sensitivity',   '-dpng', '-r150');
print(fig9, 'results/figures/scenario4/shrinking_ball',   '-dpng', '-r150');
print(fig10,'results/figures/scenario4/ball_control_demand','-dpng', '-r150');
fprintf('PNG figures saved to results/figures/scenario4/\n\n');

%% JSON export 
json_data.runA1 = make_run_entry(resA1, 'A1', 'pure_energy', true, 0, N, ...
    dt, x0, x_switch, u_max, con_A1.terminal_ball.r_switch, false, resA2.dv, true);

json_data.runA2 = make_run_entry(resA2, 'A2', 'pure_energy', false, 0, N, ...
    dt, x0, x_switch, u_max, NaN, false, NaN, true);

json_data.runB_noball = make_run_entry(resB_noball, 'B_noball', ...
    sprintf('energy_rho_c_%.0e', rho_c), false, rho_c, N_B_noball, ...
    dt, x0, x_switch, u_max, NaN, false, NaN, true);

json_data.runB_ball = make_run_entry(resB_ball, 'B_ball', ...
    sprintf('energy_rho_c_%.0e', rho_c), true, rho_c, N, ...
    dt, x0, x_switch, u_max, r_switch_ball_B, false, dX_B_max, true);

json_data.runC = make_run_entry(resC, 'C', 'energy_clf_terminal', false, 0, N, ...
    dt, x0, x_switch, u_max, NaN, false, NaN, false);

json_data.runMPC = make_run_entry(resMPC, 'MPC', 'standard_regulation_Qdare_Rdare', ...
    false, 0, N, dt, x0, x_switch, u_max, NaN, false, NaN, false);

json_data.sweep_horizon = struct( ...
    'N',           sw_N, ...
    'rho_c',       sw_rho_c, ...
    'converged',   sw_conv, ...
    'steps',       sw_steps, ...
    'dv_ms',       sw_dv, ...
    'term_dist_m', sw_td, ...
    'note', 'n_steps_override=60 for all; sigma(N)~sqrt(delta/(rho_c*N))');

json_str = jsonencode(json_data, 'PrettyPrint', true);
fid = fopen('exports/scenarios/sim_scenario4_closing.json', 'w');
fprintf(fid, '%s', json_str);
fclose(fid);
fprintf('JSON exported to exports/scenarios/sim_scenario4_closing.json\n');

%%  Local functions

function s = make_run_entry(res, run_label, cost_label, ball_on, rho_c_val, N_used, ...
                             dt_val, x0_val, x_switch_val, u_max_val, r_ball_val, ...
                             use_ps, ball_evidence, v_f_zero)
    if nargin < 14, v_f_zero = true; end
    s.metadata = struct( ...
        'run',                run_label, ...
        'cost',               cost_label, ...
        'ball_on',            ball_on, ...
        'rho_c',              rho_c_val, ...
        'N',                  N_used, ...
        'dt',                 dt_val, ...
        'x0',                 x0_val', ...
        'x_switch',           x_switch_val', ...
        'u_max',              u_max_val, ...
        'r_ball',             r_ball_val, ...
        'use_passive_safety', use_ps, ...
        'V_f_zero',           v_f_zero);
    s.results = struct( ...
        'total_dv_ms',            res.dv, ...
        'n_steps_to_convergence', res.conv_step, ...
        'passive_safe_fraction',  res.ps_frac, ...
        'terminal_distance_m',    res.term_dist, ...
        'ball_evidence',          ball_evidence, ...
        'n_ps_rows_avg',          mean(res.n_ps_h), ...
        'n_ps_rows_max',          max([res.n_ps_h, 0]));
    s.turnpike = struct( ...
        'u_thresh',          1e-4 * u_max_val, ...
        'turnpike_steps',    res.turnpike_steps, ...
        'turnpike_fraction', res.turnpike_fraction, ...
        'tp_label',          res.tp_label, ...
        'u_norm_sequence',   res.u_norm_seq);
    s.t       = res.t_vec;
    s.X       = res.X';
    s.U       = res.U';
    s.ef      = res.ef;
    s.is_safe = res.is_safe;
end
