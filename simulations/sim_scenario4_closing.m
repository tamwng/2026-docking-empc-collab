% sim_scenario4_closing.m
% Scenario 4 — Passively Safe Closing Phase
%
% EMPC (V_f=0) long-range closing from x0=[0;10000;0;0;0;0] m to the
% switch point x_switch=[0;300;0;0;0;0] m.  Stage cost ℓ(u)=‖u‖² (energy).
% Turnpike property (Grüne 2013) ensures near-optimal performance without
% a stabilising terminal cost.  Three optional hard constraints: averaged
% power cap, terminal ball, and SCA-linearised passive safety (KOS
% avoidance under free drift).

clear; clc;

%% Path setup
script_dir = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(script_dir, '..', 'tools', 'matlab2tikz-master', 'src')));

%% Constants
constants;

%% Simulation and Scenario parameters
dt       = 120;                      % [s] sampling time
N        = 20;                       % prediction horizon
x0       = [0; 10000; 0; 0; 0; 0];  % [m, m/s] initial condition
x_switch = [0;   300; 0; 0; 0; 0];  % [m] handover point (V-bar, 300 m)
u_max    = 1e-2;                     % [m/s²] thrust bound
r_KOS    = 50;                       % [m] keep-out sphere radius

%% Derived parameters
t_final = 5 * T;              % [s] ~7.7 h (5 orbital periods)
n_steps = round(t_final / dt);
N_safe  = 554;                % free-drift steps (~9.2 h, ~6 orbital periods)

%% CWH state space
Ac = [0      0     0    1     0    0  ;
      0      0     0    0     1    0  ;
      0      0     0    0     0    1  ;
      3*n^2  0     0    0     2*n  0  ;
      0      0     0   -2*n   0    0  ;
      0      0    -n^2  0     0    0 ];
Bc = [zeros(3,3); eye(3)];
[Ad, Bd] = discretize(Ac, Bc, dt);

%% Precompute Ad powers for post-hoc passive safety check
Ad_pow_safe = zeros(6, 6, N_safe);
Ad_cur = Ad;
for j = 1:N_safe
    Ad_pow_safe(:,:,j) = Ad_cur;
    Ad_cur = Ad_cur * Ad;
end

%% Precompute CW STMs for passive safety SCA constraint
fprintf('Precomputing CW STMs (%d matrices)... ', N_safe);
C_pos  = [eye(3), zeros(3,3)];
Phi_ps = zeros(6, 6, N_safe);
for i = 1:N_safe
    Phi_ps(:,:,i) = cw_stm(n, i * dt);
end
fprintf('done.\n');

%% EMPC cost matrices  (V_f = 0, energy stage cost)
Q_eco = zeros(6);
R_eco = eye(3);
P_eco = zeros(6);

%% Constraint toggles (for variations)
use_avg_power      = true;
P_avail            = 1e-5;  % [m²/s⁴] averaged power budget
use_terminal_ball  = true;
r_switch_ball      = 200;   % [m] terminal ball radius (L-inf)
use_passive_safety = false;
ps_margin          = 50;    % [m] passive safety activation margin

%% Constraint struct
con.u_max           = u_max;
con.y_min_active    = false;
con.los_cone.active = false;

if use_avg_power
    con.avg_power.P_avail = P_avail;
    con.avg_power.U_warm  = zeros(3*N, 1);  % warm start, updated each step
end

if use_terminal_ball
    con.terminal_ball.active   = true;
    con.terminal_ball.r_switch = r_switch_ball;
end

if use_passive_safety
    con.passive_safety.active        = true;
    con.passive_safety.Phi_array     = Phi_ps;
    con.passive_safety.C_pos         = C_pos;
    con.passive_safety.r_KOS         = r_KOS;
    con.passive_safety.safety_margin = ps_margin;
    con.passive_safety.x_switch      = x_switch;
    con.passive_safety.U_warm        = zeros(3*N, 1);  % warm start, updated each step
end

%% Allocate
U_warm     = zeros(3*N, 1);
X          = zeros(6, n_steps + 1);
U          = zeros(3, n_steps + 1);
ef         = zeros(1, n_steps + 1);
tf         = zeros(1, n_steps + 1);
flags      = zeros(1, n_steps + 1);
avg_power  = zeros(1, n_steps + 1);
n_ps_h     = zeros(1, n_steps);
u_norm_seq = zeros(1, n_steps);
all_X_pred = cell(1, n_steps);

X(:,1)    = x0;
e_k       = x0 - x_switch;
conv_step = NaN;
n_run     = 0;

%% EMPC loop
fprintf('\n  %5s  %10s  %12s  %8s  %4s  %6s  %8s\n', ...
    'step', '||u||', 'avg_power', 'active', 'n_ps', 'flag', 'time[s]');

for k = 1:n_steps
    ef(k) = norm(e_k);

    if ef(k) < 50
        conv_step = k;
        n_run     = k;
        break;
    end

    if use_avg_power
        con.avg_power.U_warm = U_warm;
    end
    if use_passive_safety
        con.passive_safety.U_warm = U_warm;
    end

    t0 = tic;
    [u_opt, U_opt, X_pred_k, exitflag, n_ps_step] = mpc_regulation( ...
        e_k, Ad, Bd, Q_eco, R_eco, P_eco, N, con);
    tf(k)         = toc(t0);
    flags(k)      = exitflag;
    n_ps_h(k)     = n_ps_step;
    u_norm_seq(k) = norm(u_opt);
    all_X_pred{k} = X_pred_k + x_switch;   % error → Hill frame

    avg_power(k) = (1/N) * norm(U_opt)^2;
    active = use_avg_power && (avg_power(k) >= 0.999 * P_avail);

    fprintf('  %5d  %10.3e  %12.3e  %8s  %4d  %6d  %8.4f\n', ...
        k, u_norm_seq(k), avg_power(k), mat2str(active), n_ps_step, exitflag, tf(k));

    U_warm   = [U_opt(4:end); U_opt(end-2:end)];
    U(:,k)   = u_opt;
    X(:,k+1) = Ad * X(:,k) + Bd * u_opt;
    e_k      = X(:,k+1) - x_switch;
    n_run    = k + 1;
end

%% Trim
X          = X(:,       1:n_run);
U          = U(:,       1:n_run);
ef         = ef(        1:n_run);
tf         = tf(        1:n_run);
flags      = flags(     1:n_run);
avg_power  = avg_power( 1:n_run);
n_ctrl     = n_run - 1;
n_ps_h     = n_ps_h(    1:n_ctrl);
u_norm_seq = u_norm_seq(1:n_ctrl);
all_X_pred = all_X_pred(1:n_ctrl);

%% Post-hoc passive safety check
is_safe   = check_passive_safety(X, n_run, Ad_pow_safe, r_KOS, N_safe);
ps_frac   = mean(is_safe);
dv        = sum(vecnorm(U(:,1:end-1), 2, 1)) * dt;
term_dist = norm(X(:,end) - x_switch);
t_vec     = (0:n_run-1) * dt;
t_h       = t_vec / 3600;

%% Summary
fprintf('\n── Summary ──────────────────────────────────────────────────\n');
if ~isnan(conv_step)
    fprintf('  Converged at step %d  (%.2f h)\n', conv_step, t_vec(conv_step)/3600);
else
    fprintf('  DID NOT CONVERGE in %d steps.\n', n_steps);
end
fprintf('  Total Δv:           %.4f m/s\n', dv);
fprintf('  Terminal distance:  %.2f m\n',   term_dist);
fprintf('  Passive safe:       %.1f%%\n',   ps_frac * 100);
if use_passive_safety && any(n_ps_h > 0)
    fprintf('  PS rows — avg: %.1f / max: %d\n', mean(n_ps_h), max(n_ps_h));
end
bad = find(flags(1:n_ctrl) ~= 1 & flags(1:n_ctrl) ~= 0);
if isempty(bad)
    fprintf('  All QP solves:      optimal.\n');
else
    fprintf('  Non-optimal QPs at: %s\n', mat2str(bad));
end
fprintf('─────────────────────────────────────────────────────────────\n\n');

%% Turnpike detection (used for figures 6 and 7)
u_thresh          = 1e-4 * u_max;
tp_mask           = u_norm_seq < u_thresh;
turnpike_steps    = find(tp_mask);
turnpike_fraction = length(turnpike_steps) / n_ctrl;
tp_label          = 'Turnpike (near-zero thrust)';

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

%% Figures
c = [0.85 0.33 0.10];   % burnt orange

% Figure 1 — Distance to switch point
fig1 = figure('Name', 'S4 — Distance to Switch Point');
ax = axes(fig1);
hold(ax, 'on'); grid(ax, 'on');
semilogy(ax, t_h, ef, 'Color', c, 'LineWidth', 1.4, 'DisplayName', 'EMPC');
yline(ax, 50, '--k', 'LineWidth', 0.9, 'DisplayName', '50 m threshold');
xlabel(ax, 'Time [h]');
ylabel(ax, '$\|e_k\|$ [m]', 'Interpreter', 'latex');
title(ax, 'Distance to Switch Point');
legend(ax, 'Location', 'southwest');

if exist('matlab2tikz', 'file')
    matlab2tikz('results/figures/scenario4/dist_convergence.tikz', ...
        'figurehandle', fig1, 'showInfo', false);
end

% Figure 2 — Thrust norm + averaged power
fig2 = figure('Name', 'S4 — Control Effort and Averaged Power');

ax2a = subplot(2,1,1);
hold(ax2a, 'on'); grid(ax2a, 'on');
plot(ax2a, t_h(1:end-1), vecnorm(U(:,1:end-1), 2, 1), ...
    'Color', c, 'LineWidth', 1.2);
yline(ax2a, u_max, '--k', 'LineWidth', 0.8, 'HandleVisibility', 'off');
xlabel(ax2a, 'Time [h]');
ylabel(ax2a, '$\|u_k\|$ [m/s$^2$]', 'Interpreter', 'latex');
title(ax2a, 'Per-step Thrust Norm');

ax2b = subplot(2,1,2);
hold(ax2b, 'on'); grid(ax2b, 'on');
plot(ax2b, t_h(1:end-1), avg_power(1:end-1), ...
    'Color', c, 'LineWidth', 1.2);
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

% Figure 3 — Cumulative delta-v
fig3 = figure('Name', 'S4 — Cumulative Delta-v');
ax3 = axes(fig3);
hold(ax3, 'on'); grid(ax3, 'on');
plot(ax3, t_h(1:end-1), cumsum(vecnorm(U(:,1:end-1), 2, 1)) * dt, ...
    'Color', c, 'LineWidth', 1.4, 'DisplayName', 'EMPC');
xlabel(ax3, 'Time [h]');
ylabel(ax3, 'Cumulative $\Delta v$ [m/s]', 'Interpreter', 'latex');
title(ax3, 'Cumulative $\Delta v$', 'Interpreter', 'latex');
legend(ax3, 'Location', 'northwest');

if exist('matlab2tikz', 'file')
    matlab2tikz('results/figures/scenario4/cumulative_dv.tikz', ...
        'figurehandle', fig3, 'showInfo', false);
end

% Figure 4 — Passive safety status
fig4 = figure('Name', 'S4 — Passive Safety Status');
ax4 = axes(fig4);
hold(ax4, 'on'); grid(ax4, 'on');
stem(ax4, t_h, double(is_safe), 'Color', c, ...
    'LineWidth', 0.8, 'Marker', 'none');
ylim(ax4, [-0.1 1.5]);
yticks(ax4, [0 1]);  yticklabels(ax4, {'Unsafe', 'Safe'});
xlabel(ax4, 'Time [h]');
title(ax4, 'Passive Safety Status — EMPC');

if exist('matlab2tikz', 'file')
    matlab2tikz('results/figures/scenario4/passive_safety_status.tikz', ...
        'figurehandle', fig4, 'showInfo', false);
end

% Figure 5 — Hill frame trajectory
fig5 = figure('Name', 'S4 — Hill Frame Trajectory');

ax5a = subplot(1,2,1);
hold(ax5a, 'on'); grid(ax5a, 'on'); axis(ax5a, 'equal');
plot(ax5a, X(2,:), X(1,:), 'Color', c, 'LineWidth', 1.4, 'DisplayName', 'EMPC');
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
mask_zoom = X(2,:) <= y_zoom;
plot(ax5b, X(2, mask_zoom), X(1, mask_zoom), 'Color', c, ...
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

% Figure 6 — Turnpike: thrust profile (log scale)
fig6 = figure('Name', 'S4 — Turnpike: Thrust Profile');
ax6 = axes(fig6);
hold(ax6, 'on'); grid(ax6, 'on');
set(ax6, 'YScale', 'log');

k_vec  = 1:n_ctrl;
u_plot = max(u_norm_seq, 1e-20);
stem(ax6, k_vec, u_plot, 'Color', c, 'LineWidth', 0.8, ...
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

% Figure 7 — Open-loop predictions vs closed-loop trajectory
colors_pred = [0.12 0.47 0.71;
               0.20 0.63 0.17;
               0.89 0.10 0.11;
               0.60 0.31 0.64];

fig7 = figure('Name', 'S4 — Turnpike: Open-Loop Predictions');
ax7 = axes(fig7);
hold(ax7, 'on'); grid(ax7, 'on'); axis(ax7, 'equal');
plot(ax7, X(2,:), X(1,:), '-k', 'LineWidth', 2.0, 'DisplayName', 'CL trajectory');

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

%% JSON export
metadata = struct( ...
    'scenario',           'S4_closing', ...
    'dynamics_model',     'CWH', ...
    'dt',                 dt, ...
    'N',                  N, ...
    'u_max',              u_max, ...
    'x0',                 x0', ...
    'x_switch',           x_switch', ...
    'r_KOS',              r_KOS, ...
    'N_safe',             N_safe, ...
    'use_avg_power',      use_avg_power, ...
    'P_avail',            P_avail, ...
    'use_terminal_ball',  use_terminal_ball, ...
    'r_switch_ball',      r_switch_ball, ...
    'use_passive_safety', use_passive_safety, ...
    'ps_margin',          ps_margin, ...
    'V_f_zero',           true);

results = struct( ...
    'total_dv_ms',            dv, ...
    'n_steps_to_convergence', conv_step, ...
    'passive_safe_fraction',  ps_frac, ...
    'terminal_distance_m',    term_dist, ...
    'n_ps_rows_avg',          mean(n_ps_h), ...
    'n_ps_rows_max',          max(n_ps_h));

turnpike_data = struct( ...
    'u_thresh',          u_thresh, ...
    'turnpike_steps',    turnpike_steps, ...
    'turnpike_fraction', turnpike_fraction, ...
    'u_norm_sequence',   u_norm_seq);

data = struct( ...
    'metadata',  metadata, ...
    'results',   results, ...
    'turnpike',  turnpike_data, ...
    't',         t_vec, ...
    'X',         X', ...
    'U',         U', ...
    'ef',        ef, ...
    'avg_power', avg_power, ...
    'is_safe',   is_safe);

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


%% ── Local functions ──────────────────────────────────────────────────────────

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
