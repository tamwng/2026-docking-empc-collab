% sim_scenario4_closing.m
% Scenario 4 — Passively Safe Closing Phase.
%
% Compares Standard MPC vs Economic MPC for a long-range closing manoeuvre
% from 10 km to the 300 m V-bar hold point (S2 handover). Both controllers
% operate on the error state e_k = x_k - x_switch, driving e -> 0.
% Passive safety is evaluated post-hoc: for each logged state, the free-drift
% trajectory is propagated for N_safe steps and checked against a keep-out
% sphere of radius r_KOS around the target.
%
% Cost formulation:
%   Std MPC:  min  sum Q*e + u'R u  +  e_N'P e_N
%   EMPC:     min  sum u'R_eco u    +  e_N'P e_N
%   Both use the same P from DARE(Ad,Bd,Q_dare,R_dare).

clear; clc;

%% Ensure matlab2tikz is on path
script_dir = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(script_dir, '..', 'tools', 'matlab2tikz-master', 'src')));

%% Constants
constants;

%% Simulation parameters
dt      = 120;           % [s]
t_final = 5 * T;        % [s]  ~7.7 hours
n_steps = round(t_final / dt);

%% MPC parameters
N     = 20;             % prediction horizon [-]

%% Constraints
con.u_max        = 1e-2;    % [m/s^2]  symmetric thrust bound
con.y_min_active = false;
con.los_cone.active     = false;
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
x0       = [0; 10000; 0; 0; 0; 0];   % [m, m, m, m/s, m/s, m/s]
x_switch = [0;   300; 0; 0; 0; 0];   % S2 handover — 300 m on V-bar

%% Terminal cost  (DARE with Std MPC weights)
Q_dare = diag([1e-2, 1e-2, 1e-2, 1, 1, 1]);
R_dare = diag([1e4,  1e4,  1e4]);
[~, P] = lqr_controller(Ad, Bd, Q_dare, R_dare);

%% Passive safety parameters
r_KOS  = 50;    % [m]  keep-out sphere radius (50 meters instead of 200 should be sufficent)
N_safe = 554;    % free-drift steps (~9.2 h, ~6 orbital periods)

% Pre-compute Ad powers for free-drift propagation  [6 x 6 x N_safe]
Ad_pow_safe = zeros(6, 6, N_safe);
Ad_cur = Ad;
for j = 1:N_safe
    Ad_pow_safe(:,:,j) = Ad_cur;
    Ad_cur = Ad_cur * Ad;
end

%% STANDARD MPC
fprintf('Running Standard MPC...\n');

Q_std = Q_dare;
R_std = R_dare;

X_std    = zeros(6, n_steps);
U_std    = zeros(3, n_steps);
ef_std   = zeros(1, n_steps);   % ||e_k||
tf_std   = zeros(1, n_steps);   % solve time
flags_std = zeros(1, n_steps);  % quadprog exit flags

X_std(:,1) = x0;
e_k        = x0 - x_switch;
conv_std   = NaN;
n_std      = 0;

for k = 1:n_steps
    ef_std(k) = norm(e_k);

    if ef_std(k) < 50
        conv_std = k;
        n_std    = k;
        break;
    end

    t0 = tic;
    [u_opt, ~, ~, exitflag] = mpc_regulation( ...
        e_k, Ad, Bd, Q_std, R_std, P, N, con);
    tf_std(k) = toc(t0);
    flags_std(k) = exitflag;

    U_std(:,k)    = u_opt;
    X_std(:,k+1)  = Ad * X_std(:,k) + Bd * u_opt;
    e_k           = X_std(:,k+1) - x_switch;
    n_std         = k + 1;
end

X_std  = X_std(:, 1:n_std);
U_std  = U_std(:, 1:n_std);
ef_std = ef_std(1:n_std);
tf_std = tf_std(1:n_std);
flags_std = flags_std(1:n_std);
t_std  = (0:n_std-1) * dt;

fprintf('  done (%d steps).\n', n_std);


%% ECONOMIC MPC
fprintf('Running Economic MPC...\n');

Q_eco = zeros(6);
R_eco = eye(3);

X_eco    = zeros(6, n_steps);
U_eco    = zeros(3, n_steps);
ef_eco   = zeros(1, n_steps);
tf_eco   = zeros(1, n_steps);
flags_eco = zeros(1, n_steps);

X_eco(:,1) = x0;
e_k        = x0 - x_switch;
conv_eco   = NaN;
n_eco      = 0;

for k = 1:n_steps
    ef_eco(k) = norm(e_k);

    if ef_eco(k) < 50
        conv_eco = k;
        n_eco    = k;
        break;
    end

    t0 = tic;
    [u_opt, ~, ~, exitflag] = mpc_regulation( ...
        e_k, Ad, Bd, Q_eco, R_eco, P, N, con);
    tf_eco(k) = toc(t0);
    flags_eco(k) = exitflag;

    U_eco(:,k)    = u_opt;
    X_eco(:,k+1)  = Ad * X_eco(:,k) + Bd * u_opt;
    e_k           = X_eco(:,k+1) - x_switch;
    n_eco         = k + 1;
end

X_eco  = X_eco(:, 1:n_eco);
U_eco  = U_eco(:, 1:n_eco);
ef_eco = ef_eco(1:n_eco);
tf_eco = tf_eco(1:n_eco);
flags_eco = flags_eco(1:n_eco);
t_eco  = (0:n_eco-1) * dt;

fprintf('  done (%d steps).\n', n_eco);


%% PASSIVE SAFETY DIAGNOSTIC
fprintf('Running passive safety diagnostic...\n');

is_safe_std = check_passive_safety(X_std, n_std, Ad_pow_safe, r_KOS, N_safe);
is_safe_eco = check_passive_safety(X_eco, n_eco, Ad_pow_safe, r_KOS, N_safe);

ps_frac_std = mean(is_safe_std);
ps_frac_eco = mean(is_safe_eco);

first_viol_std = find(~is_safe_std, 1, 'first');
first_viol_eco = find(~is_safe_eco, 1, 'first');


%% SANITY CHECKS 
dv_std = sum(vecnorm(U_std, 2, 1)) * dt;
dv_eco = sum(vecnorm(U_eco, 2, 1)) * dt;
dv_saving = (dv_std - dv_eco) / dv_std * 100;

fprintf('\n══════════════ SANITY CHECKS ══════════════\n');

fprintf('\n[1] Convergence\n');
if ~isnan(conv_std)
    fprintf('  Std MPC:  converged at step %d  (%.2f h)\n', conv_std, t_std(conv_std)/3600);
else
    fprintf('  Std MPC:  DID NOT CONVERGE within t_final — consider increasing t_final or retuning weights.\n');
end
if ~isnan(conv_eco)
    fprintf('  EMPC:     converged at step %d  (%.2f h)\n', conv_eco, t_eco(conv_eco)/3600);
else
    fprintf('  EMPC:     DID NOT CONVERGE within t_final — consider increasing t_final or retuning weights.\n');
end

fprintf('\n[2] Fuel comparison\n');
fprintf('  Std MPC  total Δv: %.4f m/s\n', dv_std);
fprintf('  EMPC     total Δv: %.4f m/s\n', dv_eco);
if dv_saving >= 0
    fprintf('  EMPC fuel saving: %.2f%%\n', dv_saving);
else
    fprintf('  WARNING: EMPC used MORE fuel than Std MPC (%.2f%%). Check Q/R weights.\n', -dv_saving);
end

fprintf('\n[3] Passive safety\n');
fprintf('  Std MPC:  %.1f%% of steps passively safe\n', ps_frac_std * 100);
fprintf('  EMPC:     %.1f%% of steps passively safe\n', ps_frac_eco * 100);
if ps_frac_std < 0.05 && ps_frac_eco < 0.05
    fprintf('  NOTE: Both near zero — expected at 10 km without a hard safety constraint.\n');
    fprintf('        This motivates adding a hard KOS constraint as a future extension.\n');
end
if ~isempty(first_viol_std)
    fprintf('  Std MPC:  first passive safety violation at step %d (%.2f h)\n', ...
        first_viol_std, t_std(first_viol_std)/3600);
end
if ~isempty(first_viol_eco)
    fprintf('  EMPC:     first passive safety violation at step %d (%.2f h)\n', ...
        first_viol_eco, t_eco(first_viol_eco)/3600);
end

fprintf('\n[4] Solver health\n');
bad_std = find(flags_std ~= 1 & flags_std ~= 0);
bad_eco = find(flags_eco ~= 1 & flags_eco ~= 0);
if isempty(bad_std)
    fprintf('  Std MPC:  all solves optimal.\n');
else
    fprintf('  Std MPC:  non-optimal exits at steps: %s\n', mat2str(bad_std));
end
if isempty(bad_eco)
    fprintf('  EMPC:     all solves optimal.\n');
else
    fprintf('  EMPC:     non-optimal exits at steps: %s\n', mat2str(bad_eco));
end
fprintf('═══════════════════════════════════════════\n\n');


%% FIGURES
c_std = [0.00 0.45 0.70];   % blue
c_eco = [0.85 0.33 0.10];   % red/orange

t_std_h = t_std / 3600;
t_eco_h = t_eco / 3600;

%% Figure 1 — Distance to switch point
fig1 = figure('Name', 'S4 — Distance to Switch Point');
ax = axes(fig1);
hold(ax, 'on'); grid(ax, 'on');
semilogy(ax, t_std_h, ef_std, 'Color', c_std, 'LineWidth', 1.4, 'DisplayName', 'Std MPC');
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

%% Figure 2 — Control effort
fig2 = figure('Name', 'S4 — Control Effort');

ax2a = subplot(2,1,1);
hold(ax2a, 'on'); grid(ax2a, 'on');
plot(ax2a, t_std_h(1:end-1), vecnorm(U_std(:,1:end-1), 2, 1), ...
    'Color', c_std, 'LineWidth', 1.2);
yline(ax2a, con.u_max, '--k', 'LineWidth', 0.8);
xlabel(ax2a, 'Time [h]');
ylabel(ax2a, '$\|u_k\|$ [m/s$^2$]', 'Interpreter', 'latex');
title(ax2a, 'Control Effort — Std MPC');

ax2b = subplot(2,1,2);
hold(ax2b, 'on'); grid(ax2b, 'on');
plot(ax2b, t_eco_h(1:end-1), vecnorm(U_eco(:,1:end-1), 2, 1), ...
    'Color', c_eco, 'LineWidth', 1.2);
yline(ax2b, con.u_max, '--k', 'LineWidth', 0.8);
xlabel(ax2b, 'Time [h]');
ylabel(ax2b, '$\|u_k\|$ [m/s$^2$]', 'Interpreter', 'latex');
title(ax2b, 'Control Effort — EMPC');

if exist('matlab2tikz', 'file')
    matlab2tikz('results/figures/scenario4/control_effort.tikz', ...
        'figurehandle', fig2, 'showInfo', false);
end

%% Figure 3 — Cumulative delta-v
fig3 = figure('Name', 'S4 — Cumulative Delta-v');
ax3 = axes(fig3);
hold(ax3, 'on'); grid(ax3, 'on');
plot(ax3, t_std_h(1:end-1), cumsum(vecnorm(U_std(:,1:end-1), 2, 1)) * dt, ...
    'Color', c_std, 'LineWidth', 1.4, 'DisplayName', 'Std MPC');
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

ax4a = subplot(2,1,1);
hold(ax4a, 'on'); grid(ax4a, 'on');
stem(ax4a, t_std_h, double(is_safe_std), 'Color', c_std, ...
    'LineWidth', 0.8, 'Marker', 'none');
ylim(ax4a, [-0.1 1.5]);
yticks(ax4a, [0 1]); yticklabels(ax4a, {'Unsafe', 'Safe'});
xlabel(ax4a, 'Time [h]');
title(ax4a, 'Passive Safety — Std MPC');

ax4b = subplot(2,1,2);
hold(ax4b, 'on'); grid(ax4b, 'on');
stem(ax4b, t_eco_h, double(is_safe_eco), 'Color', c_eco, ...
    'LineWidth', 0.8, 'Marker', 'none');
ylim(ax4b, [-0.1 1.5]);
yticks(ax4b, [0 1]); yticklabels(ax4b, {'Unsafe', 'Safe'});
xlabel(ax4b, 'Time [h]');
title(ax4b, 'Passive Safety — EMPC');

if exist('matlab2tikz', 'file')
    matlab2tikz('results/figures/scenario4/passive_safety_status.tikz', ...
        'figurehandle', fig4, 'showInfo', false);
end

%% Figure 5 — Hill frame trajectory
fig5 = figure('Name', 'S4 — Hill Frame Trajectory');

%% Subplot 1: full approach (10 km → switch point)
ax5a = subplot(1,2,1);
hold(ax5a, 'on'); grid(ax5a, 'on'); axis(ax5a, 'equal');
plot(ax5a, X_std(2,:), X_std(1,:), 'Color', c_std, 'LineWidth', 1.4, ...
    'DisplayName', 'Std MPC');
plot(ax5a, X_eco(2,:), X_eco(1,:), 'Color', c_eco, 'LineWidth', 1.4, ...
    'DisplayName', 'EMPC');
plot(ax5a, x0(2),       x0(1),       'o', 'Color', [0.4 0.4 0.4], ...
    'MarkerFaceColor', [0.4 0.4 0.4], 'MarkerSize', 6, 'HandleVisibility', 'off');
plot(ax5a, x_switch(2), x_switch(1), 's', 'Color', [0.2 0.6 0.2], ...
    'MarkerFaceColor', [0.2 0.6 0.2], 'MarkerSize', 7, 'DisplayName', 'Switch point');
plot(ax5a, 0, 0, '+k', 'MarkerSize', 9, 'LineWidth', 1.5, 'DisplayName', 'Target');
xlabel(ax5a, 'Along-track $y$ [m]', 'Interpreter', 'latex');
ylabel(ax5a, 'Radial $x$ [m]',      'Interpreter', 'latex');
title(ax5a, 'Full Approach (10 km $\rightarrow$ switch)', 'Interpreter', 'latex');
legend(ax5a, 'Location', 'northeast');

%% Subplot 2: zoomed final approach — switch region + KOS
ax5b = subplot(1,2,2);
hold(ax5b, 'on'); grid(ax5b, 'on'); axis(ax5b, 'equal');

% KOS circle
theta_kos = linspace(0, 2*pi, 200);
plot(ax5b, r_KOS * cos(theta_kos), r_KOS * sin(theta_kos), ...
    '--', 'Color', [0.8 0.0 0.0], 'LineWidth', 1.0, 'DisplayName', ...
    sprintf('KOS r=%d m', r_KOS));

% Trajectories — only the portion within the zoom window
y_zoom = 1500;   % [m] show up to 1500 m along-track
mask_std = X_std(2,:) <= y_zoom;
mask_eco = X_eco(2,:) <= y_zoom;
plot(ax5b, X_std(2, mask_std), X_std(1, mask_std), 'Color', c_std, ...
    'LineWidth', 1.4, 'DisplayName', 'Std MPC');
plot(ax5b, X_eco(2, mask_eco), X_eco(1, mask_eco), 'Color', c_eco, ...
    'LineWidth', 1.4, 'DisplayName', 'EMPC');

% Switch point and target markers
plot(ax5b, x_switch(2), x_switch(1), 's', 'Color', [0.2 0.6 0.2], ...
    'MarkerFaceColor', [0.2 0.6 0.2], 'MarkerSize', 7, 'DisplayName', 'Switch point');
plot(ax5b, 0, 0, '+k', 'MarkerSize', 9, 'LineWidth', 1.5, 'DisplayName', 'Target');

xlabel(ax5b, 'Along-track $y$ [m]', 'Interpreter', 'latex');
ylabel(ax5b, 'Radial $x$ [m]',      'Interpreter', 'latex');
title(ax5b, 'Final Approach — Switch Region', 'Interpreter', 'latex');
legend(ax5b, 'Location', 'northeast');
xlim(ax5b, [-50 y_zoom]);

if exist('matlab2tikz', 'file')
    matlab2tikz('results/figures/scenario4/hill_frame.tikz', ...
        'figurehandle', fig5, 'showInfo', false);
end


%% JSON EXPOR
metadata = struct( ...
    'scenario',       'S4_closing', ...
    'dynamics_model', 'CWH', ...
    'dt',             dt, ...
    'N',              N, ...
    'u_max',          con.u_max, ...
    'x0',             x0', ...
    'x_switch',       x_switch', ...
    'r_KOS',          r_KOS, ...
    'N_safe',         N_safe);

std_mpc = struct( ...
    'total_dv_ms',             dv_std, ...
    'n_steps_to_convergence',  conv_std, ...
    'passive_safety_fraction', ps_frac_std, ...
    'Q',                       'diag([1e-2,1e-2,1e-2,1,1,1])', ...
    'R',                       'diag([1e4,1e4,1e4])');

empc = struct( ...
    'total_dv_ms',             dv_eco, ...
    'n_steps_to_convergence',  conv_eco, ...
    'passive_safety_fraction', ps_frac_eco, ...
    'Q_stage',                 'zeros(6)', ...
    'R_empc',                  'eye(3)');

data = struct( ...
    'metadata', metadata, ...
    'std_mpc',  std_mpc, ...
    'empc',     empc, ...
    't_std',    t_std, ...
    't_eco',    t_eco, ...
    'X_std',    X_std', ...
    'X_eco',    X_eco', ...
    'U_std',    U_std', ...
    'U_eco',    U_eco', ...
    'ef_std',   ef_std, ...
    'ef_eco',   ef_eco, ...
    'is_safe_std', is_safe_std, ...
    'is_safe_eco', is_safe_eco);

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

%% ── Local functions ───────────────────────────────────────────────────────

function is_safe = check_passive_safety(X_log, n_logged, Ad_pow_safe, r_KOS, N_safe)
    is_safe = false(1, n_logged);
    for k = 1:n_logged
        xk = X_log(:, k);
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
