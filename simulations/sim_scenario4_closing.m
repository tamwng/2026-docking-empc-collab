% sim_scenario4_closing.m
% Scenario 4 — Passively Safe Closing Phase (EMPC, V_f=0, averaged power constraint)
%
% Long-range closing transfer from IC=[0;10000;0;0;0;0] m to switch point
% x_sw=[0;300;0;0;0;0] m.  EMPC stage cost ℓ(u)=||u||² (energy minimisation).
%
% V_f=0 + averaged power constraint achieves near-optimal performance via the
% turnpike property (Grüne 2013, Angeli 2012).
%
% Run A: baseline EMPC (no passive safety in QP)
% Run B: passive safety added as linearised KOS half-space per prediction step
%        (SCA linearisation, constraint active when predicted free-drift < 150 m).

clear; clc;

%% Ensure matlab2tikz is on path
script_dir = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(script_dir, '..', 'tools', 'matlab2tikz-master', 'src')));

%% Constants
constants;

%% Simulation parameters
dt      = 120;       % [s]
t_final = 5 * T;    % [s]  ~7.7 hours
n_steps = round(t_final / dt);

%% MPC parameters
N = 20;   % prediction horizon [-]

%% Scenario 4 specific parameters
P_avail  = 1e-5;   % [m²/s⁴] averaged power budget — non-binding
r_switch = 200;    % [m] terminal ball radius (L-inf box approximation of L2 ball)

%% Constraints (base struct)
con.u_max            = 1e-2;
con.y_min_active     = false;
con.los_cone.active  = false;
con.los_cone.half_angle = deg2rad(20);
con.los_cone.n_faces    = 10;

%% CWH state space
Ac = [0      0     0    1     0    0  ;
      0      0     0    0     1    0  ;
      0      0     0    0     0    1  ;
      3*n^2  0     0    0     2*n  0  ;
      0      0     0   -2*n   0    0  ;
      0      0    -n^2  0     0    0 ];
Bc = [zeros(3,3); eye(3)];
[Ad, Bd] = discretize(Ac, Bc, dt);

%% Initial condition and switch point
x0       = [0; 10000; 0; 0; 0; 0];   % [m]
x_switch = [0;   300; 0; 0; 0; 0];   % S2 handover — 300 m on V-bar

%% DARE terminal cost (retained for std-MPC path; NOT used by EMPC)
Q_dare = diag([1e-2, 1e-2, 1e-2, 1, 1, 1]);
R_dare = diag([1e4,  1e4,  1e4]);
[~, P] = lqr_controller(Ad, Bd, Q_dare, R_dare);   %#ok<ASGLU>

%% Passive safety parameters
r_KOS  = 50;     % [m]  keep-out sphere radius
N_safe = 554;    % free-drift steps (~9.2 h, ~6 orbital periods)

Ad_pow_safe = zeros(6, 6, N_safe);
Ad_cur = Ad;
for j = 1:N_safe
    Ad_pow_safe(:,:,j) = Ad_cur;
    Ad_cur = Ad_cur * Ad;
end

%% Precompute CW STMs for passive safety SCA constraint
fprintf('Precomputing CW STMs (%d matrices)... ', N_safe);
C_pos   = [eye(3), zeros(3,3)];
Phi_ps  = zeros(6, 6, N_safe);
for i = 1:N_safe
    Phi_ps(:,:,i) = cw_stm(n, i * dt);
end
fprintf('done.\n');
safety_margin = 2 * r_KOS;   % 100 m — activation band beyond r_KOS

%% EMPC cost matrices (shared by both runs)
Q_eco = zeros(6);
R_eco = eye(3);
P_eco = zeros(6);   % V_f=0

%% Base EMPC constraint struct
con_eco                        = con;
con_eco.avg_power.P_avail      = P_avail;
con_eco.avg_power.U_warm       = zeros(3*N, 1);   % reset before each run
con_eco.terminal_ball.active   = true;
con_eco.terminal_ball.r_switch = r_switch;

%% Pack shared parameters for the run function
p.x0          = x0;
p.x_switch    = x_switch;
p.Ad          = Ad;
p.Bd          = Bd;
p.Q           = Q_eco;
p.R_ctrl      = R_eco;
p.P           = P_eco;
p.N           = N;
p.n_steps     = n_steps;
p.dt          = dt;
p.r_KOS       = r_KOS;
p.N_safe      = N_safe;
p.Ad_pow_safe = Ad_pow_safe;

%% RUN A: BASELINE (no passive safety in QP)
con_A = con_eco;
con_A.passive_safety.active = false;
p.con_eco = con_A;
p.label   = 'RUN A — BASELINE (no passive safety constraint)';

fprintf('\n%s\n', repmat('=',1,60));
fprintf('  P_avail = %.2e  |  r_switch = %d m  |  V_f = 0 (turnpike)\n', ...
    P_avail, r_switch);
fprintf('%s\n\n', repmat('=',1,60));
RA = run_empc_s4(p);

%% Verify: minimum free-drift distance at step 14 using CW STM
fprintf('\n[Step 14 free-drift check]\n');
if size(RA.X, 2) >= 14
    x14  = RA.X(:, 14);
    d14  = inf;
    for ii = 1:N_safe
        p14 = C_pos * cw_stm(n, ii*dt) * x14;
        if norm(p14) < d14; d14 = norm(p14); end
    end
    fprintf('  Actual state at step 14: x=%.1f m, y=%.1f m  (||pos||=%.1f m)\n', ...
        x14(1), x14(2), norm(x14(1:3)));
    fprintf('  Min free-drift distance: %.2f m  (r_KOS=%d m, margin=%d m)\n', ...
        d14, r_KOS, safety_margin);
    if d14 < r_KOS + safety_margin
        fprintf('  → Within activation band: constraint would fire.\n');
    else
        fprintf('  → Outside activation band: constraint silent at this step.\n');
    end
end

%% RUN B: WITH PASSIVE SAFETY CONSTRAINT
con_B = con_eco;
con_B.passive_safety.active        = true;
con_B.passive_safety.Phi_array     = Phi_ps;
con_B.passive_safety.C_pos         = C_pos;
con_B.passive_safety.r_KOS         = r_KOS;
con_B.passive_safety.safety_margin = safety_margin;
con_B.passive_safety.x_switch      = x_switch;
con_B.passive_safety.U_warm        = zeros(3*N, 1);

p.con_eco = con_B;
p.label   = 'RUN B — WITH PASSIVE SAFETY CONSTRAINT';
RB = run_empc_s4(p);

%% Check: if any Run B QP was infeasible, re-run with tighter margin
if ~isempty(RB.bad_steps)
    fprintf('\n[Feasibility recovery] Re-running Run B with safety_margin = 0.5*r_KOS = %d m\n', ...
        round(0.5 * r_KOS));
    con_B2 = con_B;
    con_B2.passive_safety.safety_margin = 0.5 * r_KOS;
    p.con_eco = con_B2;
    p.label   = 'RUN B2 — PASSIVE SAFETY (tight margin = 0.5*r_KOS)';
    RB2 = run_empc_s4(p);
    if isempty(RB2.bad_steps)
        fprintf('  Feasibility recovered with tighter margin.\n');
        RB = RB2;   % promote tight-margin run as the "B" result
    else
        fprintf('  Still infeasible at steps: %s\n', mat2str(RB2.bad_steps));
    end
end

%% ── ALIAS RUN A → EXISTING VARIABLE NAMES (for sanity-checks / figures) ──
X_eco          = RA.X;
U_eco          = RA.U;
ef_eco         = RA.ef;
tf_eco         = RA.tf;
flags_eco      = RA.flags;
avg_power_used = RA.apow;
u_norm_seq     = RA.u_norm_seq;
all_X_pred     = RA.all_X_pred;
is_safe_eco    = RA.is_safe;
ps_frac_eco    = RA.ps_frac;
first_viol_eco = RA.first_viol;
conv_eco       = RA.conv;
n_eco          = RA.n_run;
n_ctrl         = RA.n_ctrl;
dv_eco         = RA.dv;
term_dist      = RA.term_dist;
t_eco          = RA.t;
avg_power_active = any(avg_power_used(1:end-1) >= 0.999 * P_avail);


%% SANITY CHECKS (Run A)
fprintf('\n══════════════════ SANITY CHECKS (RUN A) ══════════════════\n');

fprintf('\n[1] Convergence\n');
if ~isnan(conv_eco)
    fprintf('  EMPC:  converged at step %d  (%.2f h)\n', conv_eco, t_eco(conv_eco)/3600);
else
    fprintf('  EMPC:  DID NOT CONVERGE.\n');
end
fprintf('  Terminal distance ||x_final - x_switch|| = %.2f m\n', term_dist);

fprintf('\n[2] Fuel\n');
fprintf('  EMPC total Δv: %.4f m/s\n', dv_eco);

fprintf('\n[3] Averaged power constraint  (P_avail = %.2e)\n', P_avail);
fprintf('  Max avg power used: %.2e m^2/s^4\n', max(avg_power_used(1:end-1)));
if avg_power_active
    fprintf('  Constraint was ACTIVE at >= 1 step.\n');
else
    fprintf('  Constraint was NOT active (non-binding at this P_avail).\n');
end

fprintf('\n[4] Passive safety\n');
fprintf('  EMPC: %.1f%% of steps passively safe\n', ps_frac_eco * 100);
if ~isempty(first_viol_eco)
    fprintf('  WARNING: first violation at step %d (%.2f h)\n', ...
        first_viol_eco, t_eco(first_viol_eco)/3600);
else
    fprintf('  All steps passively safe.\n');
end

fprintf('\n[5] Solver health\n');
bad_eco = find(flags_eco ~= 1 & flags_eco ~= 0);
if isempty(bad_eco)
    fprintf('  EMPC: all solves optimal.\n');
else
    fprintf('  EMPC: non-optimal exits at steps: %s\n', mat2str(bad_eco));
end
fprintf('════════════════════════════════════════════════════════════\n\n');


%% TURNPIKE ANALYSIS (Run A)
u_thresh         = 1e-4 * con.u_max;
tp_mask          = u_norm_seq < u_thresh;
turnpike_steps   = find(tp_mask);
turnpike_fraction = length(turnpike_steps) / n_ctrl;
tp_label         = 'Turnpike (near-zero thrust)';

fprintf('[6] Turnpike (u_thresh = %.2e m/s^2)\n', u_thresh);
fprintf('  Turnpike steps: %d / %d  (%.1f%%)\n', ...
    length(turnpike_steps), n_ctrl, turnpike_fraction * 100);
if ~isempty(turnpike_steps)
    fprintf('  Turnpike region: step %d to step %d\n', ...
        turnpike_steps(1), turnpike_steps(end));
else
    fprintf('  No near-zero coast at this threshold.\n');
end

if isempty(turnpike_steps) && n_ctrl >= 5
    log_u   = log(max(u_norm_seq, 1e-20));
    dlog_u  = abs(diff(log_u));
    pl_mask = [dlog_u < 0.15, false];
    pl_pad  = [false, pl_mask, false];
    d_pl    = diff(pl_pad);
    pl_st   = find(d_pl == 1);
    pl_en   = find(d_pl == -1) - 1;
    if ~isempty(pl_st)
        [~, li]           = max(pl_en - pl_st);
        turnpike_steps    = pl_st(li):pl_en(li);
        turnpike_fraction = length(turnpike_steps) / n_ctrl;
        tp_label          = 'Turnpike (steady cruise)';
        fprintf('  Plateau detected (log-rate < 15%%): step %d to %d (%.1f%%).\n', ...
            turnpike_steps(1), turnpike_steps(end), turnpike_fraction * 100);
    else
        fprintf('  No turnpike structure detected.\n');
    end
end

%% Select 4 steps for prediction visualisation
if ~isempty(turnpike_steps) && length(turnpike_steps) >= 2
    i_tp_start = turnpike_steps(1);
    i_tp_mid   = turnpike_steps(round(end/2));
    i_tp_end   = turnpike_steps(end);
    i_entry    = max(1, i_tp_start - 1);
    i_exit     = min(n_ctrl, i_tp_end + 1);
    vis_steps  = unique([i_entry, i_tp_start, i_tp_mid, i_exit]);
elseif ~isempty(turnpike_steps)
    vis_steps = unique([1, turnpike_steps(1), ...
        round(n_ctrl * 0.6), min(n_ctrl, turnpike_steps(1) + 1)]);
else
    vis_steps = unique(round(linspace(1, n_ctrl, 4)));
end
while length(vis_steps) < 4 && length(vis_steps) < n_ctrl
    candidates = setdiff(1:n_ctrl, vis_steps);
    if isempty(candidates); break; end
    vis_steps  = sort([vis_steps, candidates(round(length(candidates) / 2))]);
end
vis_steps = vis_steps(1:min(4, end));
fprintf('  Visualisation steps: %s\n\n', mat2str(vis_steps));


%% ── PASSIVE SAFETY COMPARISON REPORT ─────────────────────────────────────
fprintf('══════════════ PASSIVE SAFETY COMPARISON ══════════════\n\n');

runs    = {RA, RB};
labels  = {'A (baseline)', 'B (with PS constraint)'};
for ri = 1:2
    R   = runs{ri};
    lbl = labels{ri};
    fprintf('  Run %s\n', lbl);
    fprintf('    Passive safe fraction : %.1f%%\n',  R.ps_frac * 100);
    fprintf('    Total Δv              : %.4f m/s\n', R.dv);
    fprintf('    Steps to convergence  : %d\n',       R.conv);
    fprintf('    Avg solve time        : %.4f s\n',   mean(R.tf(1:R.n_ctrl)));
    if ri == 2   % passive safety stats only for Run B
        if R.n_ctrl > 0
            fprintf('    PS rows — avg: %.1f / max: %d per step\n', ...
                mean(R.n_ps_h), max(R.n_ps_h));
            act = find(R.n_ps_h > 0);
            if ~isempty(act)
                fprintf('    Steps with PS rows: %s\n', mat2str(act));
            else
                fprintf('    No PS rows added at any step.\n');
            end
        end
    end
    if ~isempty(R.bad_steps)
        fprintf('    *** Non-optimal QP exits at steps: %s\n', mat2str(R.bad_steps));
        for bi = 1:length(R.bad_steps)
            bk = R.bad_steps(bi);
            if bk <= size(R.X, 2)
                xbk = R.X(:, bk);
                d_bk = inf;
                for ii = 1:N_safe
                    p_bk = C_pos * cw_stm(n, ii*dt) * xbk;
                    if norm(p_bk) < d_bk; d_bk = norm(p_bk); end
                end
                fprintf('      Step %d: min_dist = %.2f m\n', bk, d_bk);
            end
        end
    else
        fprintf('    All QP solves optimal.\n');
    end
    fprintf('\n');
end
fprintf('════════════════════════════════════════════════════════\n\n');


%% FIGURES (Run A data throughout)
c_eco   = [0.85 0.33 0.10];
t_eco_h = t_eco / 3600;

%% Figure 1 — Distance to switch point
fig1 = figure('Name', 'S4 — Distance to Switch Point');
ax = axes(fig1);
hold(ax, 'on'); grid(ax, 'on');
semilogy(ax, t_eco_h, ef_eco, 'Color', c_eco, 'LineWidth', 1.4, 'DisplayName', 'EMPC');
yline(ax, 50, '--k', 'LineWidth', 0.9, 'DisplayName', '50 m threshold');
xlabel(ax, 'Time [h]');
ylabel(ax, '$\|e_k\|$ [m]', 'Interpreter', 'latex');
title(ax, 'Distance to Switch Point');
legend(ax, 'Location', 'southwest');

if exist('matlab2tikz', 'file')
    matlab2tikz('results/figures/scenario4/dist_convergence.tikz', ...
        'figurehandle', fig1, 'showInfo', false);
end

%% Figure 2 — Thrust norm + averaged power
fig2 = figure('Name', 'S4 — Control Effort and Averaged Power');

ax2a = subplot(2,1,1);
hold(ax2a, 'on'); grid(ax2a, 'on');
plot(ax2a, t_eco_h(1:end-1), vecnorm(U_eco(:,1:end-1), 2, 1), ...
    'Color', c_eco, 'LineWidth', 1.2);
yline(ax2a, con.u_max, '--k', 'LineWidth', 0.8, 'HandleVisibility', 'off');
xlabel(ax2a, 'Time [h]');
ylabel(ax2a, '$\|u_k\|$ [m/s$^2$]', 'Interpreter', 'latex');
title(ax2a, 'Per-step Thrust Norm (bang–coast–bang structure expected)');

ax2b = subplot(2,1,2);
hold(ax2b, 'on'); grid(ax2b, 'on');
plot(ax2b, t_eco_h(1:end-1), avg_power_used(1:end-1), ...
    'Color', c_eco, 'LineWidth', 1.2);
yline(ax2b, P_avail, '--k', 'LineWidth', 0.8, 'DisplayName', ...
    sprintf('P_{avail} = %.0e', P_avail));
xlabel(ax2b, 'Time [h]');
ylabel(ax2b, '$\frac{1}{N}\|U\|^2$ [m$^2$/s$^4$]', 'Interpreter', 'latex');
title(ax2b, 'Averaged Power per Step (SCA constraint)');
legend(ax2b, 'Location', 'northeast');

if exist('matlab2tikz', 'file')
    matlab2tikz('results/figures/scenario4/control_effort.tikz', ...
        'figurehandle', fig2, 'showInfo', false);
end

%% Figure 3 — Cumulative delta-v
fig3 = figure('Name', 'S4 — Cumulative Delta-v');
ax3 = axes(fig3);
hold(ax3, 'on'); grid(ax3, 'on');
plot(ax3, t_eco_h(1:end-1), cumsum(vecnorm(U_eco(:,1:end-1), 2, 1)) * dt, ...
    'Color', c_eco, 'LineWidth', 1.4, 'DisplayName', 'EMPC');
xlabel(ax3, 'Time [h]');
ylabel(ax3, 'Cumulative $\Delta v$ [m/s]', 'Interpreter', 'latex');
title(ax3, 'Cumulative $\Delta v$', 'Interpreter', 'latex');
legend(ax3, 'Location', 'northwest');

if exist('matlab2tikz', 'file')
    matlab2tikz('results/figures/scenario4/cumulative_dv.tikz', ...
        'figurehandle', fig3, 'showInfo', false);
end

%% Figure 4 — Passive safety status
fig4 = figure('Name', 'S4 — Passive Safety Status');
ax4 = axes(fig4);
hold(ax4, 'on'); grid(ax4, 'on');
stem(ax4, t_eco_h, double(is_safe_eco), 'Color', c_eco, ...
    'LineWidth', 0.8, 'Marker', 'none');
ylim(ax4, [-0.1 1.5]);
yticks(ax4, [0 1]);  yticklabels(ax4, {'Unsafe', 'Safe'});
xlabel(ax4, 'Time [h]');
title(ax4, 'Passive Safety Status — EMPC');

if exist('matlab2tikz', 'file')
    matlab2tikz('results/figures/scenario4/passive_safety_status.tikz', ...
        'figurehandle', fig4, 'showInfo', false);
end

%% Figure 5 — Hill frame trajectory
fig5 = figure('Name', 'S4 — Hill Frame Trajectory');

ax5a = subplot(1,2,1);
hold(ax5a, 'on'); grid(ax5a, 'on'); axis(ax5a, 'equal');
plot(ax5a, X_eco(2,:), X_eco(1,:), 'Color', c_eco, 'LineWidth', 1.4, 'DisplayName', 'EMPC');
plot(ax5a, x0(2),       x0(1),       'o', 'Color', [0.4 0.4 0.4], ...
    'MarkerFaceColor', [0.4 0.4 0.4], 'MarkerSize', 6, 'HandleVisibility', 'off');
plot(ax5a, x_switch(2), x_switch(1), 's', 'Color', [0.2 0.6 0.2], ...
    'MarkerFaceColor', [0.2 0.6 0.2], 'MarkerSize', 7, 'DisplayName', 'Switch point');
plot(ax5a, 0, 0, '+k', 'MarkerSize', 9, 'LineWidth', 1.5, 'DisplayName', 'Target');
xlabel(ax5a, 'Along-track $y$ [m]', 'Interpreter', 'latex');
ylabel(ax5a, 'Radial $x$ [m]',      'Interpreter', 'latex');
title(ax5a, 'Full Approach (10 km $\rightarrow$ switch)', 'Interpreter', 'latex');
legend(ax5a, 'Location', 'northeast');

ax5b = subplot(1,2,2);
hold(ax5b, 'on'); grid(ax5b, 'on'); axis(ax5b, 'equal');
theta_kos = linspace(0, 2*pi, 200);
plot(ax5b, r_KOS * cos(theta_kos), r_KOS * sin(theta_kos), ...
    '--', 'Color', [0.8 0.0 0.0], 'LineWidth', 1.0, ...
    'DisplayName', sprintf('KOS r=%d m', r_KOS));
y_zoom = 1500;
mask_eco = X_eco(2,:) <= y_zoom;
plot(ax5b, X_eco(2, mask_eco), X_eco(1, mask_eco), 'Color', c_eco, ...
    'LineWidth', 1.4, 'DisplayName', 'EMPC');
plot(ax5b, x_switch(2), x_switch(1), 's', 'Color', [0.2 0.6 0.2], ...
    'MarkerFaceColor', [0.2 0.6 0.2], 'MarkerSize', 7, 'DisplayName', 'Switch point');
plot(ax5b, 0, 0, '+k', 'MarkerSize', 9, 'LineWidth', 1.5, 'DisplayName', 'Target');
xlabel(ax5b, 'Along-track $y$ [m]', 'Interpreter', 'latex');
ylabel(ax5b, 'Radial $x$ [m]',      'Interpreter', 'latex');
title(ax5b, 'Final Approach -- Switch Region', 'Interpreter', 'latex');
legend(ax5b, 'Location', 'northeast');
xlim(ax5b, [-50 y_zoom]);

if exist('matlab2tikz', 'file')
    matlab2tikz('results/figures/scenario4/hill_frame.tikz', ...
        'figurehandle', fig5, 'showInfo', false);
end

%% Figure 6 — Turnpike: thrust profile
fig6 = figure('Name', 'S4 — Turnpike: Thrust Profile');
ax6 = axes(fig6);
hold(ax6, 'on'); grid(ax6, 'on');
set(ax6, 'YScale', 'log');

k_vec  = 1:n_ctrl;
u_plot = max(u_norm_seq, 1e-20);
stem(ax6, k_vec, u_plot, 'Color', c_eco, 'LineWidth', 0.8, ...
    'MarkerSize', 3, 'DisplayName', '$\|u_k\|$');

if ~isempty(turnpike_steps)
    yl = ylim(ax6);
    if yl(1) <= 0 || ~isfinite(log10(yl(1)))
        pos_vals = u_plot(u_plot > 1e-20);
        yl(1) = min(pos_vals) * 0.01;
        ylim(ax6, yl); yl = ylim(ax6);
    end
    tp_lo   = turnpike_steps(1) - 0.5;
    tp_hi   = turnpike_steps(end) + 0.5;
    h_patch = patch(ax6, [tp_lo tp_hi tp_hi tp_lo], ...
        [yl(1) yl(1) yl(2) yl(2)], [0.8 0.8 0.8], ...
        'FaceAlpha', 0.35, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    uistack(h_patch, 'bottom');
    y_label = 10^(0.5 * (log10(yl(1)) + log10(yl(2))));
    if turnpike_steps(1) > 2
        text(ax6, mean([0.5, tp_lo]), y_label, 'Entry', ...
            'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', [0.2 0.2 0.2]);
    end
    text(ax6, mean([tp_lo, tp_hi]), y_label, tp_label, ...
        'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', [0.2 0.2 0.2]);
    if turnpike_steps(end) < n_ctrl - 1
        text(ax6, mean([tp_hi, n_ctrl + 0.5]), y_label, 'Exit', ...
            'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', [0.2 0.2 0.2]);
    end
end

xlabel(ax6, 'Step $k$', 'Interpreter', 'latex');
ylabel(ax6, '$\|u_k\|$ [m/s$^2$]', 'Interpreter', 'latex');
title(ax6, '$V_f=0$ EMPC, $N=20$, $P_{\rm avail}$ non-binding', 'Interpreter', 'latex');
legend(ax6, 'Interpreter', 'latex', 'Location', 'northeast');

if exist('matlab2tikz', 'file')
    matlab2tikz('results/figures/scenario4/turnpike_thrust_profile.tikz', ...
        'figurehandle', fig6, 'showInfo', false);
end

%% Figure 7 — Open-loop predictions vs closed-loop trajectory
colors_pred = [0.12 0.47 0.71;
               0.20 0.63 0.17;
               0.89 0.10 0.11;
               0.60 0.31 0.64];

fig7 = figure('Name', 'S4 — Turnpike: Open-Loop Predictions');
ax7 = axes(fig7);
hold(ax7, 'on'); grid(ax7, 'on'); axis(ax7, 'equal');
plot(ax7, X_eco(2, :), X_eco(1, :), '-k', 'LineWidth', 2.0, ...
    'DisplayName', 'CL trajectory');

for ii = 1:length(vis_steps)
    kk = vis_steps(ii);
    if kk <= length(all_X_pred) && ~isempty(all_X_pred{kk})
        XP = all_X_pred{kk};
        plot(ax7, XP(2,:), XP(1,:), '-', 'Color', colors_pred(ii,:), ...
            'LineWidth', 0.9, 'DisplayName', sprintf('Prediction $k=%d$', kk));
        plot(ax7, XP(2,1), XP(1,1), 'o', 'Color', colors_pred(ii,:), ...
            'MarkerFaceColor', colors_pred(ii,:), 'MarkerSize', 5, ...
            'HandleVisibility', 'off');
    end
end
plot(ax7, x_switch(2), x_switch(1), 's', 'Color', [0.2 0.6 0.2], ...
    'MarkerFaceColor', [0.2 0.6 0.2], 'MarkerSize', 8, 'DisplayName', 'Switch point');

% Annotate turnpike arc from first visible prediction
kk_ann = vis_steps(1);
if ~isempty(all_X_pred{kk_ann})
    XP_ann  = all_X_pred{kk_ann};
    ann_col = round(size(XP_ann, 2) * 0.4);
    text(ax7, XP_ann(2, ann_col), XP_ann(1, ann_col), '  Turnpike arc', ...
        'FontSize', 8, 'Color', [0.3 0.3 0.3], 'FontAngle', 'italic');
end

xlabel(ax7, 'Along-track $y$ [m]', 'Interpreter', 'latex');
ylabel(ax7, 'Radial $x$ [m]',      'Interpreter', 'latex');
title(ax7, 'Open-Loop Predictions vs Closed-Loop Trajectory', 'Interpreter', 'latex');
legend(ax7, 'Location', 'northeast', 'Interpreter', 'latex');

if exist('matlab2tikz', 'file')
    matlab2tikz('results/figures/scenario4/turnpike_predictions.tikz', ...
        'figurehandle', fig7, 'showInfo', false);
end


%% JSON EXPORT (both runs)
metadata = struct( ...
    'scenario',       'S4_closing', ...
    'dynamics_model', 'CWH', ...
    'dt',             dt, ...
    'N',              N, ...
    'u_max',          con.u_max, ...
    'x0',             x0', ...
    'x_switch',       x_switch', ...
    'r_KOS',          r_KOS, ...
    'N_safe',         N_safe, ...
    'P_avail',        P_avail, ...
    'V_f_zero',       true, ...
    'r_switch',       r_switch, ...
    'safety_margin',  safety_margin);

empc_results = struct( ...
    'total_dv_ms',            dv_eco, ...
    'n_steps_to_convergence', conv_eco, ...
    'passive_safe_fraction',  ps_frac_eco, ...
    'avg_power_active',       avg_power_active, ...
    'terminal_distance_m',    term_dist, ...
    'Q_stage',                'zeros(6)', ...
    'R_empc',                 'eye(3)', ...
    'P_terminal',             'zeros(6)');

empc_ps_results = struct( ...
    'total_dv_ms',            RB.dv, ...
    'n_steps_to_convergence', RB.conv, ...
    'passive_safe_fraction',  RB.ps_frac, ...
    'terminal_distance_m',    RB.term_dist, ...
    'n_ps_rows_avg',          mean(RB.n_ps_h), ...
    'n_ps_rows_max',          max(RB.n_ps_h), ...
    'ps_active_steps',        RB.ps_active_steps, ...
    'bad_steps',              RB.bad_steps);

turnpike_data = struct( ...
    'u_thresh',          u_thresh, ...
    'turnpike_steps',    turnpike_steps, ...
    'turnpike_fraction', turnpike_fraction, ...
    'u_norm_sequence',   u_norm_seq);

data = struct( ...
    'metadata',           metadata, ...
    'empc_baseline',      empc_results, ...
    'empc_passive_safety',empc_ps_results, ...
    'turnpike',           turnpike_data, ...
    't_eco',              t_eco, ...
    'X_eco',              X_eco', ...
    'U_eco',              U_eco', ...
    'ef_eco',             ef_eco, ...
    'avg_power_used',     avg_power_used, ...
    'is_safe_eco',        is_safe_eco);

json_str = jsonencode(data, 'PrettyPrint', true);
fid = fopen('exports/scenarios/sim_scenario4_closing.json', 'w');
fprintf(fid, '%s', json_str);
fclose(fid);

fprintf('JSON exported to exports/scenarios/sim_scenario4_closing.json\n');
if exist('matlab2tikz', 'file')
    fprintf('TikZ figures exported to results/figures/scenario4/\n');
else
    fprintf('TikZ export skipped (matlab2tikz not on path — run setup.m first).\n');
end


%% ── LOCAL FUNCTIONS ──────────────────────────────────────────────────────

function R = run_empc_s4(p)
% Runs the EMPC loop for Scenario 4 and returns all results in struct R.
% p fields: x0, x_switch, Ad, Bd, Q, R_ctrl, P, N, n_steps, dt,
%           con_eco, r_KOS, N_safe, Ad_pow_safe, label.

fprintf('\n%s\n\n', p.label);

U_warm = zeros(3 * p.N, 1);
X      = zeros(6, p.n_steps + 1);
U      = zeros(3, p.n_steps + 1);
ef     = zeros(1, p.n_steps + 1);
tf     = zeros(1, p.n_steps + 1);
flags  = zeros(1, p.n_steps + 1);
apow   = zeros(1, p.n_steps + 1);
n_ps_h = zeros(1, p.n_steps);
u_norm = zeros(1, p.n_steps);
xpred  = cell(1, p.n_steps);

X(:,1) = p.x0;
e_k    = p.x0 - p.x_switch;
conv   = NaN;
n_run  = 0;
con    = p.con_eco;

fprintf('  %5s  %10s  %12s  %8s  %4s  %6s  %8s\n', ...
    'step', '||u||', 'avg_power', 'active', 'n_ps', 'flag', 'time[s]');

for k = 1:p.n_steps
    ef(k) = norm(e_k);

    if ef(k) < 50
        conv  = k;
        n_run = k;
        break;
    end

    con.avg_power.U_warm = U_warm;
    if isfield(con, 'passive_safety') && con.passive_safety.active
        con.passive_safety.U_warm = U_warm;
    end

    t0 = tic;
    [u_opt, U_opt, X_pred_k, exitflag, n_ps_step] = mpc_regulation( ...
        e_k, p.Ad, p.Bd, p.Q, p.R_ctrl, p.P, p.N, con);
    tf(k)     = toc(t0);
    flags(k)  = exitflag;
    n_ps_h(k) = n_ps_step;
    u_norm(k) = norm(u_opt, 2);
    xpred{k}  = X_pred_k + p.x_switch;   % error → Hill frame

    apow(k) = (1 / p.N) * norm(U_opt)^2;
    active  = apow(k) >= 0.999 * p.con_eco.avg_power.P_avail;

    fprintf('  %5d  %10.3e  %12.3e  %8s  %4d  %6d  %8.4f\n', ...
        k, norm(u_opt), apow(k), mat2str(active), n_ps_step, exitflag, tf(k));

    U_warm = [U_opt(4:end); U_opt(end-2:end)];

    U(:,k)    = u_opt;
    X(:,k+1)  = p.Ad * X(:,k) + p.Bd * u_opt;
    e_k       = X(:,k+1) - p.x_switch;
    n_run     = k + 1;
end

X     = X(:,    1:n_run);
U     = U(:,    1:n_run);
ef    = ef(     1:n_run);
tf    = tf(     1:n_run);
flags = flags(  1:n_run);
apow  = apow(   1:n_run);
n_ctrl = n_run - 1;
n_ps_h = n_ps_h(1:n_ctrl);
u_norm = u_norm(1:n_ctrl);
xpred  = xpred( 1:n_ctrl);

is_safe   = check_passive_safety(X, n_run, p.Ad_pow_safe, p.r_KOS, p.N_safe);
ps_frac   = mean(is_safe);
first_viol = find(~is_safe, 1, 'first');

fprintf('\n  done: %d steps,  max solve time = %.3f s\n', n_run, max(tf));

R.X           = X;
R.U           = U;
R.ef          = ef;
R.tf          = tf;
R.flags       = flags;
R.apow        = apow;
R.u_norm_seq  = u_norm;
R.all_X_pred  = xpred;
R.is_safe     = is_safe;
R.ps_frac     = ps_frac;
R.first_viol  = first_viol;
R.conv        = conv;
R.n_run       = n_run;
R.n_ctrl      = n_ctrl;
R.n_ps_h      = n_ps_h;
R.dv          = sum(vecnorm(U(:,1:end-1), 2, 1)) * p.dt;
R.term_dist   = norm(X(:,end) - p.x_switch);
R.ps_active_steps = find(n_ps_h > 0);
R.bad_steps   = find(flags ~= 1 & flags ~= 0);
R.t           = (0:n_run-1) * p.dt;
end


function is_safe = check_passive_safety(X_log, n_logged, Ad_pow_safe, r_KOS, N_safe)
    is_safe = false(1, n_logged);
    for k = 1:n_logged
        xk   = X_log(:, k);
        safe = true;
        for j = 1:N_safe
            x_free = Ad_pow_safe(:,:,j) * xk;
            if norm(x_free(1:3)) < r_KOS
                safe = false;
                break;
            end
        end
        is_safe(k) = safe;
    end
end
