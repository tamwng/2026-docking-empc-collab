% sim_scenario4_closing.m
% Scenario 4 — terminal-ball study (Run A) and strict-dissipativity study (Run B).
%
% Run A — pure energy cost ℓ=‖u‖², Q=0, P=0.
%   No dissipativity certificate: marginally stable CWH zero-fuel modes tie at
%   the optimum, so the optimizer coasts and u≡0 is globally optimal without
%   a terminal constraint.
%   A1 (ball ON):  converges — ball is the sole forcing mechanism.
%                  Ref: Grüne (2013) Thm 6.4; Angeli, Rawlings & Amrit (2012) Thm 1.
%   A2 (ball OFF): u≡0 at every step, chaser stuck at IC — proves ball load-bearing.
%
% Run B — energy + small state regularisation ℓ=‖u‖²+ρc‖e‖², λ≡0, α_ℓ=ρc·r².
%   Strict dissipativity is restored trivially; turnpike drives chaser toward
%   x_switch autonomously.  Terminal ball becomes slack.
%   B-noball: converges without ball — Grüne (2013) Thm 5.6/4.2, Remark 7.7.
%   B-ball:   ball ON but inactive (active-frac ≈ 0) — confirms slackness.

clear; clc;

%% Path setup
script_dir = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(script_dir, '..', 'tools', 'matlab2tikz-master', 'src')));

%% Constants
constants;

%% Shared scenario parameters
dt        = 120;                      % [s] sampling time
N         = 20;                       % default prediction horizon
x0        = [0; 10000; 0; 0; 0; 0];  % [m, m/s] initial condition
x_switch  = [0;   300; 0; 0; 0; 0];  % [m] handover point (V-bar, 300 m)
u_max     = 1e-2;                     % [m/s²] thrust bound
r_KOS     = 50;                       % [m] keep-out sphere radius
P_avail   = 1e-5;                     % [m²/s⁴] avg-power budget (non-binding)
rho_c     = 1e-6;                     % state regularisation weight (Run B)
ps_margin = 150;                       % [m] passive-safety activation margin
t_final   = 5 * T;                    % [s] ~7.7 h (5 orbital periods)

%% Cost structs
cost_A = struct('Q', zeros(6),    'R', eye(3), 'P', zeros(6));
cost_B = struct('Q', rho_c*eye(6),'R', eye(3), 'P', zeros(6));

%% Shared params base
params_base.x0       = x0;
params_base.x_switch = x_switch;
params_base.dt       = dt;
params_base.N        = N;
params_base.u_max    = u_max;
params_base.t_final  = t_final;
params_base.r_KOS    = r_KOS;

%% ── Run A1: pure energy, terminal ball ON ────────────────────────────────────
% Ball is the ONLY mechanism forcing convergence.
% Ref: Grüne Thm 6.4 (ball closes avg-performance gap);
%      Angeli, Rawlings & Amrit (2012) Thm 1 (avg-performance lower bound).

con_A1.u_max                  = u_max;
con_A1.y_min_active           = false;
con_A1.los_cone.active        = false;
con_A1.use_avg_power          = true;
con_A1.avg_power.P_avail      = P_avail;
con_A1.use_terminal_ball      = true;
con_A1.terminal_ball.r_switch = 200;   % [m] L-inf position ball
con_A1.use_passive_safety     = false;

params_A1 = params_base;

fprintf('\n══════════════════════════════════════════════\n');
fprintf(' Run A1 — pure energy ℓ=‖u‖², ball ON\n');
fprintf('══════════════════════════════════════════════\n');
resA1 = run_closing(cost_A, con_A1, params_A1);

%% ── Run A2: pure energy, terminal ball OFF ───────────────────────────────────
% Without ball, min ‖u‖² → U=0.  Chaser stays at IC.

con_A2.u_max              = u_max;
con_A2.y_min_active       = false;
con_A2.los_cone.active    = false;
con_A2.use_avg_power      = true;
con_A2.avg_power.P_avail  = P_avail;
con_A2.use_terminal_ball  = false;
con_A2.use_passive_safety = false;

params_A2                  = params_base;
params_A2.n_steps_override = 60;

fprintf('\n══════════════════════════════════════════════\n');
fprintf(' Run A2 — pure energy ℓ=‖u‖², ball OFF\n');
fprintf(' Expect: u≡0, chaser stuck at IC.\n');
fprintf('══════════════════════════════════════════════\n');
resA2 = run_closing(cost_A, con_A2, params_A2);

max_disp_A2 = max(vecnorm(resA2.X(1:3,:) - x0(1:3), 2, 1));
max_u_A2    = max([resA2.u_norm_seq, 0]);
fprintf('A2 check — max position displacement: %.4f m   (expect ≈ 0)\n', max_disp_A2);
fprintf('         — max ||u||:                 %.2e m/s² (expect ≈ 0)\n\n', max_u_A2);

%% ── Run B-noball: strict dissipativity, no terminal ball ─────────────────────
% Adding ρc‖e‖² restores strict dissipativity (storage fn λ≡0, α_ℓ(r)=ρc·r²).
% Turnpike drives chaser to x_switch autonomously — no terminal constraint needed.
% Ref: Grüne (2013) Thm 5.6 / Thm 4.2; Grüne & Pannek (2011) Remark 7.7.

N_B_noball = N;

con_B_noball.u_max              = u_max;
con_B_noball.y_min_active       = false;
con_B_noball.los_cone.active    = false;
con_B_noball.use_avg_power      = true;
con_B_noball.avg_power.P_avail  = P_avail;
con_B_noball.use_terminal_ball  = false;
con_B_noball.use_passive_safety = false;

params_B_noball   = params_base;
params_B_noball.N = N_B_noball;

fprintf('\n══════════════════════════════════════════════\n');
fprintf(' Run B-noball — ρc=%.0e, N=%d, ball OFF\n', rho_c, N_B_noball);
fprintf(' Expect: converges (strict dissipativity, turnpike)\n');
fprintf('══════════════════════════════════════════════\n');
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

%% ── Run B-ball: strict dissipativity + terminal ball ON ─────────────────────
% Ball should be SLACK: turnpike already drives convergence, so the terminal
% constraint never binds.  Active-fraction ≈ 0 confirms this.

r_switch_ball_B = 200;   % [m]

con_B_ball = con_B_noball;
con_B_ball.use_terminal_ball      = true;
con_B_ball.terminal_ball.r_switch = r_switch_ball_B;

params_B_ball = params_base;   % N=20 (ball + strict dissip. guarantee convergence)

fprintf('\n══════════════════════════════════════════════\n');
fprintf(' Run B-ball — ρc=%.0e, N=20, ball ON r=%d m\n', rho_c, r_switch_ball_B);
fprintf(' Expect: converges; ball inactive (active-frac ≈ 0%%)\n');
fprintf('══════════════════════════════════════════════\n');
resB_ball = run_closing(cost_B, con_B_ball, params_B_ball);

%% ── Ball-role evidence (trajectory comparison) ────────────────────────────────
% Ball is load-bearing in A1 iff removing it stalls the chaser.
%   Proof: Run A2 (same cost, no ball) → U*≡0, Δv≈0.
% Ball is slack in B-ball iff the trajectory is unchanged without it.
%   Proof: ||pos(B-noball) − pos(B-ball)|| at machine precision.
n_B_cmp  = min(resB_noball.n_ctrl, resB_ball.n_ctrl);
dX_B_max = max(abs(resB_noball.X(1:3, 1:n_B_cmp+1) - resB_ball.X(1:3, 1:n_B_cmp+1)), [], 'all');

fprintf('\nBall-role evidence:\n');
fprintf('  A2 total Δv (no ball, Q=P=0):     %.2e m/s  → ball is load-bearing in A1\n', resA2.dv);
fprintf('  ||pos(B-noball) - pos(B-ball)||:  %.2e m   → ball is slack in B-ball\n\n', dX_B_max);

%% ── MPC comparison ───────────────────────────────────────────────────────────
fprintf('Running standard regulation MPC (comparison baseline)...\n');
p_mpc = struct('x0', x0, 'x_switch', x_switch, 'dt', dt, ...
               'N', N, 'u_max', u_max, 't_final', t_final);
[X_mpc, U_mpc, t_mpc, dv_mpc] = sim_mpc_regulation(p_mpc);

fprintf('\n── Δv summary ───────────────────────────────────────────────\n');
fprintf('  %-28s  %10s\n', 'Run', 'Δv [m/s]');
fprintf('  %-28s  %10.4f\n', 'A1 — energy, ball ON',            resA1.dv);
fprintf('  %-28s  %10.4f\n', sprintf('B-noball (N=%d)',N_B_noball), resB_noball.dv);
fprintf('  %-28s  %10.4f\n', 'B-ball  — strict dissip.',         resB_ball.dv);
fprintf('  %-28s  %10.4f\n', 'MPC — standard regulation',        dv_mpc);
fprintf('  A1 vs MPC savings: %.4f m/s  (%.1f%%)\n', ...
    dv_mpc - resA1.dv, 100*(dv_mpc - resA1.dv)/dv_mpc);
fprintf('─────────────────────────────────────────────────────────────\n\n');

%% ── Run A1-ps: energy + ball ON + passive safety SCA ────────────────────────
% SCA per-step constraint forces free-drift safety at every predicted step.
% Should recover ps_frac ≈ 100%.

con_A1ps                             = con_A1;
con_A1ps.use_passive_safety          = true;
con_A1ps.passive_safety.r_KOS        = r_KOS;
con_A1ps.passive_safety.safety_margin = ps_margin;
con_A1ps.passive_safety.x_switch     = x_switch;
% Phi_array lazily precomputed inside run_closing

fprintf('\n══════════════════════════════════════════════\n');
fprintf(' Run A1-ps — energy, ball ON, passive safety SCA ON\n');
fprintf(' Expect: ps_frac ≈ 100%%\n');
fprintf('══════════════════════════════════════════════\n');
resA1_ps = run_closing(cost_A, con_A1ps, params_A1);
fprintf('A1-ps: passive_safe_fraction = %.1f%%\n\n', 100*resA1_ps.ps_frac);

%% ── Horizon / ρc sweep (B-noball, 60 steps) ──────────────────────────────────
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

%% ── Figures ───────────────────────────────────────────────────────────────────
c_A1   = [0.85 0.33 0.10];   % burnt orange  — A1
c_A2   = [0.55 0.55 0.55];   % grey          — A2 (stuck)
c_Bnb  = [0.20 0.63 0.17];   % green         — B-noball
c_Bbal = [0.58 0.40 0.74];   % purple        — B-ball
c_mpc  = [0.12 0.47 0.71];   % blue          — MPC comparison

t_A1_h     = resA1.t_vec      / 3600;
t_Bnb_h    = resB_noball.t_vec / 3600;
t_Bbal_h   = resB_ball.t_vec   / 3600;
n_ctrl_A1  = resA1.n_ctrl;
n_ctrl_Bnb = resB_noball.n_ctrl;
n_ctrl_Bbal = resB_ball.n_ctrl;

%% Figure 1 — Hill-frame: A1 vs A2
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

%% Figure 2 — Cumulative Δv: A1 vs MPC
fig2 = figure('Name', 'S4 — Cumulative \Deltav: A1 vs MPC');
ax2  = axes(fig2);
hold(ax2,'on'); grid(ax2,'on');
plot(ax2, t_A1_h(1:n_ctrl_A1), ...
    cumsum(vecnorm(resA1.U(:,1:n_ctrl_A1), 2, 1))*dt, ...
    'Color', c_A1,  'LineWidth', 1.4, 'DisplayName', 'A1 — EMPC, ball ON');
plot(ax2, t_mpc(1:end-1)/3600, ...
    cumsum(vecnorm(U_mpc(:,1:end-1), 2, 1))*dt, ...
    'Color', c_mpc, 'LineWidth', 1.4, 'DisplayName', 'MPC — standard regulation');
xlabel(ax2, 'Time [h]');
ylabel(ax2, 'Cumulative $\Delta v$ [m/s]', 'Interpreter', 'latex');
title(ax2, 'Cumulative $\Delta v$ — A1 vs Standard MPC', 'Interpreter', 'latex');
legend(ax2, 'Location', 'northwest');
if exist('matlab2tikz', 'file')
    matlab2tikz('results/figures/scenario4/AB_cumulative_dv.tikz', ...
        'figurehandle', fig2, 'showInfo', false);
end

%% Figure 3 — Thrust profile A1 (log scale)
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

%% Figure 4 — Hill-frame: A1 + B-noball + B-ball
fig4 = figure('Name', 'S4 — Hill Frame: A1, B-noball, B-ball');
ax4  = axes(fig4);
hold(ax4,'on'); grid(ax4,'on'); axis(ax4,'equal');
plot(ax4, resA1.X(2,:),       resA1.X(1,:),       '-',  'Color', c_A1,  'LineWidth', 1.2, ...
    'DisplayName', 'A1 — energy, ball ON');
plot(ax4, resB_noball.X(2,:), resB_noball.X(1,:),  '-',  'Color', c_Bnb,  'LineWidth', 1.2, ...
    'DisplayName', sprintf('B-noball (\\rho_c=%.0e, N=%d)', rho_c, N_B_noball));
plot(ax4, resB_ball.X(2,:),   resB_ball.X(1,:),    '--', 'Color', c_Bbal, 'LineWidth', 1.2, ...
    'DisplayName', sprintf('B-ball (\\rho_c=%.0e, ball slack)', rho_c));
plot(ax4, x0(2), x0(1), 'o', 'Color', [0.4 0.4 0.4], ...
    'MarkerFaceColor', [0.4 0.4 0.4], 'MarkerSize', 6, 'HandleVisibility', 'off');
plot(ax4, x_switch(2), x_switch(1), 's', 'Color', [0.2 0.6 0.2], ...
    'MarkerFaceColor', [0.2 0.6 0.2], 'MarkerSize', 7, 'DisplayName', 'Switch point');
plot(ax4, 0, 0, '+k', 'MarkerSize', 9, 'LineWidth', 1.5, 'DisplayName', 'Target');
xlabel(ax4, 'Along-track $y$ [m]', 'Interpreter', 'latex');
ylabel(ax4, 'Radial $x$ [m]',      'Interpreter', 'latex');
title(ax4, 'Hill-Frame Trajectory — A1, B-noball, B-ball', 'Interpreter', 'latex');
legend(ax4, 'Location', 'northeast');
if exist('matlab2tikz', 'file')
    matlab2tikz('results/figures/scenario4/B_hill_frame.tikz', ...
        'figurehandle', fig4, 'showInfo', false);
end

%% Figure 5 — Horizon / ρc sweep: convergence steps + Δv (two panels)
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

%% Panel A: steps to convergence
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

%% Panel B: total Δv
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

%% Figure 6 — A1 vs A1-ps: Hill-frame safety coloring + thrust profiles
n_ctrl_A1ps  = resA1_ps.n_ctrl;
c_A1ps       = [0.12 0.47 0.71];   % blue  — A1-ps
c_safe_dot   = [0.18 0.55 0.34];   % green — passively safe states
c_unsafe_dot = [0.85 0.20 0.20];   % red   — passively unsafe states

fig6 = figure('Name', 'S4 — A1 vs A1-ps');
fig6.Position(3:4) = [900, 380];

%% Panel A: Hill frame, A1 states safety-coloured + A1-ps trajectory
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

%% Panel B: Thrust profiles (log scale)
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

%% PNG export
[~, ~] = mkdir('results/figures/scenario4');
print(fig1, 'results/figures/scenario4/AB_hill_frame',    '-dpng', '-r150');
print(fig2, 'results/figures/scenario4/AB_cumulative_dv', '-dpng', '-r150');
print(fig3, 'results/figures/scenario4/A1_thrust_profile','-dpng', '-r150');
print(fig4, 'results/figures/scenario4/B_hill_frame',     '-dpng', '-r150');
print(fig5, 'results/figures/scenario4/horizon_sweep',    '-dpng', '-r150');
print(fig6, 'results/figures/scenario4/A1_vs_A1ps',       '-dpng', '-r150');
fprintf('PNG figures saved to results/figures/scenario4/\n\n');

%% ── JSON export ───────────────────────────────────────────────────────────────
% One top-level object: runA1, runA2, runB_noball, runB_ball, sweep_horizon.
% Local function make_run_entry is defined at the end of the file (MATLAB
% requires local functions after all executable script code).

% ball_evidence: A1/A2 pass A2.dv (load-bearing proof); B_ball passes dX_B_max (slack proof).
json_data.runA1 = make_run_entry(resA1, 'A1', 'pure_energy', true, 0, N, ...
    dt, x0, x_switch, u_max, con_A1.terminal_ball.r_switch, false, resA2.dv);

json_data.runA2 = make_run_entry(resA2, 'A2', 'pure_energy', false, 0, N, ...
    dt, x0, x_switch, u_max, NaN, false, NaN);

json_data.runB_noball = make_run_entry(resB_noball, 'B_noball', ...
    sprintf('energy_rho_c_%.0e', rho_c), false, rho_c, N_B_noball, ...
    dt, x0, x_switch, u_max, NaN, false, NaN);

json_data.runB_ball = make_run_entry(resB_ball, 'B_ball', ...
    sprintf('energy_rho_c_%.0e', rho_c), true, rho_c, N, ...
    dt, x0, x_switch, u_max, r_switch_ball_B, false, dX_B_max);

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

%% ── Local functions (must follow all executable code) ────────────────────────

function s = make_run_entry(res, run_label, cost_label, ball_on, rho_c_val, N_used, ...
                             dt_val, x0_val, x_switch_val, u_max_val, r_ball_val, ...
                             use_ps, ball_evidence)
% ball_evidence meaning depends on run:
%   A1 : resA2.dv  — A2 Δv≈0 proves ball is load-bearing in A1
%   B_ball : dX_B_max — max pos deviation vs B_noball proves ball is slack
%   others : NaN
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
        'V_f_zero',           true);
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
