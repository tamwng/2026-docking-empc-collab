% sim_scenario1_docking.m
% Scenario 1: Final Docking Approach (ADRIOS / ClearSpace-1 analogue)
%
% Run 1 - Standard MPC (Q ≠ 0):
%   Stage cost ℓ = e'Qe + u'Ru  (Q, R from LQR design weights).
%   Terminal cost P = DARE solution (same Q, R), infinite-horizon LQR cost-to-go.
%
% Run 2 - EMPC, CLF terminal cost (Q = 0):
%   Stage cost ℓ = u'R_eco u  (pure fuel, R_eco = I).
%
% Run 3 - EMPC + hard terminal equality  x(N) = x_dock:
%   Same economic stage cost as Run 2 (Q=0, R=I) but with an explicit
%   equality constraint pinning the predicted terminal state to the docking
%   port at every receding step.

clear; clc;

%% Path setup 
script_dir = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(script_dir, '..', '..', 'tools', 'matlab2tikz-master', 'src')));

%% Constants 
constants;  

%% Shared scenario parameters 
dt     = 10;                    % [s]    close-range sampling time
N      = 24;                    % prediction horizon (min. feasible for Run 3 term.eq @150 m)
x0 = [0; 150; 0; 0; 0; 0];      % [m, m/s] V-bar hold at 150 m (start of final approach)
x_dock = zeros(6, 1);           % [m, m/s] docking port at origin
u_max  = 1e-2;                  % [m/s²] symmetric per-axis thrust bound
r_KOS  = 5;                     % [m]   keep-out sphere radius (capture ball)
N_safe = 554;                   % post-hoc free-drift horizon ≈ 1 orbital period
t_final = 1.5 * T;              % [s]   max run time (all runs converge well within this)
los_half_angle = pi/6;          % [rad] LoS approach-cone half-angle (30 deg)

%% CWH matrices (for DARE / CLF verification)
Ac_cwh = [0      0     0    1     0    0  ;
          0      0     0    0     1    0  ;
          0      0     0    0     0    1  ;
          3*n^2  0     0    0     2*n  0  ;
          0      0     0   -2*n   0    0  ;
          0      0    -n^2  0     0    0 ];
Bc_cwh = [zeros(3,3); eye(3)];
[Ad, Bd] = discretize(Ac_cwh, Bc_cwh, dt);

%% CLF terminal cost

Q_dare = diag([1e-2, 1e-2, 1e-2, 1e0, 1e0, 1e0]);
R_dare = diag([1e4,  1e4,  1e4]);
R_eco  = eye(3);   

[K_clf, P_clf] = lqr_controller(Ad, Bd, Q_dare, R_dare);

Acl_clf = Ad - Bd * K_clf;
M1      = Q_dare + K_clf' * (R_dare - R_eco) * K_clf;
M1      = (M1 + M1') / 2;
lam_M1  = min(eig(M1));
rho_Acl = max(abs(eig(Acl_clf)));

fprintf('CLF terminal cost verification (dt = %d s):\n', dt);
if lam_M1 > 1e-10 && rho_Acl < 1
    fprintf('  => P_clf satisfies Assumption 6. Valid CLF terminal cost.\n\n');
else
    warning('CLF condition not satisfied: check DARE weights.');
end

%% Cost structs 
cost_std  = struct('Q', Q_dare,   'R', R_dare, 'P', P_clf);  
cost_empc = struct('Q', zeros(6), 'R', R_eco,  'P', P_clf);    
cost_dock = struct('Q', zeros(6), 'R', R_eco,  'P', zeros(6)); 

%% Shared params base 
params_base.x0      = x0;
params_base.x_dock  = x_dock;
params_base.dt      = dt;
params_base.N       = N;
params_base.u_max   = u_max;
params_base.t_final = t_final;
params_base.r_KOS    = r_KOS;
params_base.N_safe   = N_safe;
params_base.conv_tol = 1.0;

%% Run 1 - Standard MPC

con1.u_max               = u_max;
con1.y_min_active        = true;
con1.los_cone.active     = true;
con1.los_cone.half_angle = los_half_angle;
con1.use_passive_safety  = false;

fprintf(' Run 1 - Standard MPC');
res1 = run_docking(cost_std, con1, params_base);

%% Run 2 - EMPC, CLF terminal cost  (Q = 0, R = I, P = P_clf)

con2.u_max               = u_max;
con2.y_min_active        = true;
con2.los_cone.active     = true;
con2.los_cone.half_angle = los_half_angle;
con2.use_passive_safety  = false;

fprintf(' Run 2 - EMPC  (Q = 0, R = I, P = P_clf)\n');
res2 = run_docking(cost_empc, con2, params_base);

%% Run 3 - EMPC + hard terminal equality  (Q = 0, R = I, P = 0, x(N)=x_dock)

con3.u_max               = u_max;
con3.y_min_active        = true;
con3.los_cone.active     = true;
con3.los_cone.half_angle = los_half_angle;
con3.terminal_eq         = true;
con3.use_passive_safety  = false;

fprintf(' Run 3 - EMPC + hard terminal equality\n');
res3 = run_docking(cost_dock, con3, params_base);

%% Postprocessing

% Δv and convergence summary
dv_saving_2vs1 = (res1.dv - res2.dv) / res1.dv * 100;
dv_saving_3vs1 = (res1.dv - res3.dv) / res1.dv * 100;

fprintf('\nDeltav summary (150 m -> 0 m docking, dt = %d s, N = %d):\n', dt, N);
fprintf('  %-36s  %8s  %8s  %6s\n', 'Run', 'dv [m/s]', 'steps', 'ps [%]');
fprintf('  %-36s  %8.4f  %8s  %6.1f\n', 'Run 1 - Standard MPC', ...
    res1.dv, num2str(res1.conv_step), 100*res1.ps_frac);
fprintf('  %-36s  %8.4f  %8s  %6.1f\n', 'Run 2 - EMPC, CLF terminal cost', ...
    res2.dv, num2str(res2.conv_step), 100*res2.ps_frac);
fprintf('  %-36s  %8.4f  %8s  %6.1f\n', 'Run 3 - EMPC + hard term. eq.', ...
    res3.dv, num2str(res3.conv_step), 100*res3.ps_frac);
fprintf('  Run 2 vs Run 1 dv saving: %.1f%%\n', dv_saving_2vs1);
fprintf('  Run 3 vs Run 1 dv saving: %.1f%%\n', dv_saving_3vs1);
fprintf('%s\n\n', repmat('-', 1, 67));

% Comparison table: Δv, convergence (steps / time), terminal distance
run_labels = {'Run 1 - Std MPC'; 'Run 2 - EMPC CLF'; 'Run 3 - EMPC+term.eq'};
dv_col     = [res1.dv;        res2.dv;        res3.dv];
step_col   = [res1.conv_step; res2.conv_step; res3.conv_step];
tmin_col   = step_col * dt / 60;                                  % convergence time [min]
tdist_col  = [norm(res1.X(1:3,end)); norm(res2.X(1:3,end)); ...   % terminal position dist [m]
              norm(res3.X(1:3,end))];

T_cmp = table(dv_col, step_col, tmin_col, tdist_col, ...
    'VariableNames', {'dv_ms', 'conv_steps', 'conv_time_min', 'term_dist_m'}, ...
    'RowNames', run_labels);
fprintf('%s\n', repmat('-', 1, 68));
fprintf('Comparison table\n');
fprintf('%s\n', repmat('-', 1, 68));
disp(T_cmp);
fprintf('%s\n\n', repmat('-', 1, 68));

% Colours 
c1 = [0.12 0.47 0.71];  
c2 = [0.00 0.68 0.68];   
c3 = [0.85 0.33 0.10];   

n_ctrl1 = res1.n_ctrl;
n_ctrl2 = res2.n_ctrl;
n_ctrl3 = res3.n_ctrl;

t1_min = res1.t_vec / 60;
t2_min = res2.t_vec / 60;
t3_min = res3.t_vec / 60;

%  Figure 1: Hill-frame trajectory (y-x plane), full approach + terminal zoom
% LoS cone cross-section (z = 0): docking axis along +y, edges x = ±y·tan(alpha).
alpha_deg = round(rad2deg(los_half_angle));
y_cone    = linspace(0, x0(2) * 1.05, 60);
edge_cone = y_cone * tan(los_half_angle);
cone      = struct('y', y_cone, 'edge', edge_cone, 'alpha_deg', alpha_deg);
cols      = {c1, c2, c3};
R         = {res1, res2, res3};

fig1 = figure('Name', 'S1 Docking - Hill-Frame Trajectory');
fig1.Position(3:4) = [1050, 460];

% Panel (a): full approach
ax1a = subplot(1, 2, 1);
draw_docking_scene(ax1a, R, x0, r_KOS, cols, cone, true);

% Panel (b): terminal zoom (~3× keep-out sphere around the port)
ax1b   = subplot(1, 2, 2);
draw_docking_scene(ax1b, R, x0, r_KOS, cols, cone, false);
zoom_r = 2 * r_KOS;
xlim(ax1b, [-2, zoom_r]);
ylim(ax1b, [-zoom_r, zoom_r]);
% No in-plot titles/sgtitle -- panel descriptions belong in the LaTeX
% figure caption instead.

% Figure 2: Cumulative Δv comparison
fig2 = figure('Name', 'S1 Docking - Cumulative Deltav');
ax2  = axes(fig2);
hold(ax2, 'on'); grid(ax2, 'on');

plot(ax2, t1_min(1:n_ctrl1), cumsum(vecnorm(res1.U(:,1:n_ctrl1),2,1))*dt, ...
    'Color', c1, 'LineWidth', 1.4, ...
    'DisplayName', sprintf('Run 1 - Std MPC  (%.4f m/s)', res1.dv));
plot(ax2, t2_min(1:n_ctrl2), cumsum(vecnorm(res2.U(:,1:n_ctrl2),2,1))*dt, ...
    'Color', c2, 'LineWidth', 1.4, ...
    'DisplayName', sprintf('Run 2 - EMPC, CLF $V_f$  (%.4f m/s)', res2.dv));
plot(ax2, t3_min(1:n_ctrl3), cumsum(vecnorm(res3.U(:,1:n_ctrl3),2,1))*dt, ...
    '--', 'Color', c3, 'LineWidth', 1.4, ...
    'DisplayName', sprintf('Run 3 - EMPC $+$ term.eq  (%.4f m/s)', res3.dv));

xlabel(ax2, 'Time [min]',                            'Interpreter', 'latex');
ylabel(ax2, 'Cumulative $\Delta v$ [m/s]',           'Interpreter', 'latex');
title(ax2,  'Cumulative $\Delta v$ - Docking Approach', 'Interpreter', 'latex');
legend(ax2, 'Location', 'northwest', 'Interpreter', 'latex');

% Figure 3: Thrust norm time history (log scale, 3 panels)
fig3 = figure('Name', 'S1 Docking - Thrust Profile');
fig3.Position(3:4) = [900, 350];

ax3a = subplot(1, 3, 1);
hold(ax3a, 'on'); grid(ax3a, 'on'); set(ax3a, 'YScale', 'log');
stem(ax3a, 1:n_ctrl1, max(res1.u_norm_seq, 1e-12), ...
    'Color', c1, 'LineWidth', 0.9, 'MarkerSize', 3, ...
    'DisplayName', sprintf('$\\Delta v=%.4f$ m/s', res1.dv));
xlabel(ax3a, 'Step $k$',             'Interpreter', 'latex');
ylabel(ax3a, '$\|u_k\|$ [m/s$^2$]', 'Interpreter', 'latex');
title(ax3a,  'Run 1 - Std MPC',      'Interpreter', 'latex');
legend(ax3a, 'Interpreter', 'latex', 'Location', 'northeast');

ax3b = subplot(1, 3, 2);
hold(ax3b, 'on'); grid(ax3b, 'on'); set(ax3b, 'YScale', 'log');
stem(ax3b, 1:n_ctrl2, max(res2.u_norm_seq, 1e-12), ...
    'Color', c2, 'LineWidth', 0.9, 'MarkerSize', 3, ...
    'DisplayName', sprintf('$\\Delta v=%.4f$ m/s', res2.dv));
xlabel(ax3b, 'Step $k$',                 'Interpreter', 'latex');
title(ax3b,  'Run 2 - EMPC, CLF $V_f$', 'Interpreter', 'latex');
legend(ax3b, 'Interpreter', 'latex', 'Location', 'northeast');

ax3c = subplot(1, 3, 3);
hold(ax3c, 'on'); grid(ax3c, 'on'); set(ax3c, 'YScale', 'log');
stem(ax3c, 1:n_ctrl3, max(res3.u_norm_seq, 1e-12), ...
    'Color', c3, 'LineWidth', 0.9, 'MarkerSize', 3, ...
    'DisplayName', sprintf('$\\Delta v=%.4f$ m/s', res3.dv));
xlabel(ax3c, 'Step $k$',                      'Interpreter', 'latex');
title(ax3c,  'Run 3 - EMPC $+$ term.~eq',    'Interpreter', 'latex');
legend(ax3c, 'Interpreter', 'latex', 'Location', 'northeast');

sgtitle(fig3, ...
    'Thrust profiles - Docking approach 150 m $\to$ 0 m  ($dt=10\,\mathrm{s}$, $N=24$)', ...
    'Interpreter', 'latex');

% Figure 4: Error norm decay (log scale)
fig4 = figure('Name', 'S1 Docking - Error decay');
ax4  = axes(fig4);
hold(ax4, 'on'); grid(ax4, 'on'); set(ax4, 'YScale', 'log');

plot(ax4, t1_min, res1.ef, '-',  'Color', c1, 'LineWidth', 1.4, ...
    'DisplayName', 'Run 1 - Std MPC');
plot(ax4, t2_min, res2.ef, '-',  'Color', c2, 'LineWidth', 1.4, ...
    'DisplayName', 'Run 2 - EMPC, CLF $V_f$');
plot(ax4, t3_min, res3.ef, '--', 'Color', c3, 'LineWidth', 1.4, ...
    'DisplayName', 'Run 3 - EMPC $+$ term.eq');
yline(ax4, params_base.conv_tol, '--k', 'LineWidth', 0.9, 'HandleVisibility', 'off');
text(ax4, 0, params_base.conv_tol*1.3, sprintf('conv. tol %.3g m', params_base.conv_tol), ...
    'FontSize', 8, 'Color', [0.3 0.3 0.3]);

xlabel(ax4, 'Time [min]',                        'Interpreter', 'latex');
ylabel(ax4, '$\|e_k\|$ [m]',                     'Interpreter', 'latex');
title(ax4,  'Error norm decay - Docking approach', 'Interpreter', 'latex');
legend(ax4, 'Location', 'northeast', 'Interpreter', 'latex');

% Figure 5: Δv bar chart
fig5 = figure('Name', 'S1 Docking - Deltav bar chart');
ax5  = axes(fig5);
dv_vals = [res1.dv, res2.dv, res3.dv];
b5 = bar(ax5, dv_vals);
b5.FaceColor = 'flat';
b5.CData     = [c1; c2; c3];
set(ax5, 'XTickLabel', {'Run 1  Std MPC', 'Run 2  EMPC CLF', 'Run 3  EMPC+term.eq'});
grid(ax5, 'on');
ylabel(ax5, 'Total $\Delta v$ [m/s]', 'Interpreter', 'latex');
title(ax5, sprintf('$\\Delta v$ - Docking 150 m $\\to$ 0 m  (Run 2: %.0f%% saving vs Run 1)', ...
    dv_saving_2vs1), 'Interpreter', 'latex');
for bi = 1:3
    text(ax5, bi, dv_vals(bi)*1.02, sprintf('%.4f', dv_vals(bi)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 9);
end

% Figure 6: 3-D Hill-frame trajectory with LoS approach cone
% Trajectories are planar (z = 0); the 3-D view exists to show the LoS
% constraint as an actual cone (radius = y·tan(alpha) about the +y axis).
fig6 = figure('Name', 'S1 Docking - 3D Hill Trajectory');
ax6  = axes(fig6);
hold(ax6, 'on'); grid(ax6, 'on');

% LoS cone surface, axis along +y (along-track), plotted (y, x, z)
ny = 32; nth = 30;
[YG, TH] = meshgrid(linspace(0, x0(2)*1.02, ny), linspace(0, 2*pi, nth));
RG = YG * tan(los_half_angle);
surf(ax6, YG, RG.*cos(TH), RG.*sin(TH), ...
    'FaceColor', [0.95 0.90 0.55], 'FaceAlpha', 0.15, ...
    'EdgeColor', [0.80 0.75 0.45], 'EdgeAlpha', 0.30, 'LineWidth', 0.3, ...
    'DisplayName', sprintf('LoS cone ($\\pm%d^\\circ$)', alpha_deg));

styles6 = {'-', '-', '--'};
names6  = {'Run 1 - Std MPC', 'Run 2 - EMPC, CLF $V_f$', 'Run 3 - EMPC $+$ term.eq'};
for i = 1:3
    plot3(ax6, R{i}.X(2,:), R{i}.X(1,:), R{i}.X(3,:), styles6{i}, ...
        'Color', cols{i}, 'LineWidth', 1.6, 'DisplayName', names6{i});
end

plot3(ax6, x0(2), x0(1), x0(3), 'o', 'Color', [0.4 0.4 0.4], ...
    'MarkerFaceColor', [0.4 0.4 0.4], 'MarkerSize', 7, ...
    'DisplayName', 'IC - V-bar hold (150 m)');
plot3(ax6, 0, 0, 0, 'pk', 'MarkerSize', 11, 'MarkerFaceColor', 'k', ...
    'DisplayName', 'Docking port (origin)');

xlabel(ax6, 'Along-track $y$ [m]', 'Interpreter', 'latex');
ylabel(ax6, 'Radial $x$ [m]',      'Interpreter', 'latex');
zlabel(ax6, 'Cross-track $z$ [m]', 'Interpreter', 'latex');
title(ax6,  'Hill-Frame Trajectory with LoS Cone (3-D)', 'Interpreter', 'latex');
legend(ax6, 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 8);
daspect(ax6, [1 1 1]);
view(ax6, -35, 22);

% PNG export
out_dir = 'results/figures/s1_docking';
[~, ~]  = mkdir(out_dir);
print(fig1, fullfile(out_dir, 'hill_trajectory'), '-dpng', '-r150');
print(fig6, fullfile(out_dir, 'hill_trajectory_3d'), '-dpng', '-r150');
print(fig2, fullfile(out_dir, 'cumulative_dv'),   '-dpng', '-r150');
print(fig3, fullfile(out_dir, 'thrust_profiles'), '-dpng', '-r150');
print(fig4, fullfile(out_dir, 'error_decay'),     '-dpng', '-r150');
print(fig5, fullfile(out_dir, 'dv_bar_chart'),    '-dpng', '-r150');
fprintf('PNG figures saved to %s/\n\n', out_dir);

% LaTeX comparison table (booktabs)
tex_labels = {'Standard MPC', 'EMPC (CLF $V_f$)', 'EMPC $+$ term.\ eq.'};
tex_path   = fullfile(out_dir, 'comparison_table.tex');
ftex = fopen(tex_path, 'w');
fprintf(ftex, '%% Auto-generated by sim_scenario1_docking.m, do not edit by hand.\n');
fprintf(ftex, '\\begin{tabular}{lrrrr}\n\\toprule\n');
fprintf(ftex, 'Run & $\\Delta v$ [m/s] & Steps & Time [min] & Term.\\ dist.\\ [m] \\\\\n\\midrule\n');
for i = 1:3
    fprintf(ftex, '%s & %.4f & %d & %.2f & %.3e \\\\\n', ...
        tex_labels{i}, dv_col(i), step_col(i), tmin_col(i), tdist_col(i));
end
fprintf(ftex, '\\bottomrule\n\\end{tabular}\n');
fclose(ftex);
fprintf('LaTeX comparison table saved to %s\n\n', tex_path);

% TikZ export
if exist('matlab2tikz', 'file')
    matlab2tikz(fullfile(out_dir, 'hill_trajectory.tikz'), ...
        'figurehandle', fig1, 'showInfo', false);
    matlab2tikz(fullfile(out_dir, 'hill_trajectory_3d.tikz'), ...
        'figurehandle', fig6, 'showInfo', false);
    matlab2tikz(fullfile(out_dir, 'cumulative_dv.tikz'), ...
        'figurehandle', fig2, 'showInfo', false);
    matlab2tikz(fullfile(out_dir, 'thrust_profiles.tikz'), ...
        'figurehandle', fig3, 'showInfo', false);
    matlab2tikz(fullfile(out_dir, 'error_decay.tikz'), ...
        'figurehandle', fig4, 'showInfo', false);
    fprintf('TikZ figures exported to %s/\n\n', out_dir);
else
    fprintf('TikZ export skipped (matlab2tikz not on path, run setup.m first).\n\n');
end

% JSON export
json_data.scenario = struct( ...
    'name',        'S1_docking', ...
    'description', '150m V-bar hold to docking port, ADRIOS/ClearSpace-1 analogue', ...
    'date',        char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')), ...
    'dt_s',        dt, ...
    'N',           N, ...
    'x0',          x0', ...
    'x_dock',      x_dock', ...
    'u_max',       u_max, ...
    'r_KOS_m',     r_KOS, ...
    'conv_tol_m',  params_base.conv_tol, ...
    'N_safe',      N_safe, ...
    'Q_dare',      diag(Q_dare)', ...
    'R_dare',      diag(R_dare)', ...
    'clf_lam_M1',  lam_M1, ...
    'clf_rho_Acl', rho_Acl);

json_data.run1 = make_run_entry(res1, 'Run1_StdMPC', ...
    'Q=Q_dare, R=R_dare, P=P_clf', false, N, dt, x0, x_dock, u_max);
json_data.run2 = make_run_entry(res2, 'Run2_EMPC_CLF', ...
    'Q=0, R=I, P=P_clf', false, N, dt, x0, x_dock, u_max);
json_data.run3 = make_run_entry(res3, 'Run3_EMPC_TermEq', ...
    'Q=0, R=I, P=0, terminal_eq=true', true, N, dt, x0, x_dock, u_max);

json_data.dv_comparison = struct( ...
    'run1_dv_ms',         res1.dv, ...
    'run2_dv_ms',         res2.dv, ...
    'run3_dv_ms',         res3.dv, ...
    'saving_run2vs1_pct', dv_saving_2vs1, ...
    'saving_run3vs1_pct', dv_saving_3vs1);

json_str = jsonencode(json_data, 'PrettyPrint', true);
fid = fopen('exports/scenarios/sim_scenario1_docking.json', 'w');
fprintf(fid, '%s', json_str);
fclose(fid);
fprintf('JSON exported to exports/scenarios/sim_scenario1_docking.json\n');

% Local functions

function draw_docking_scene(ax, R, x0, r_KOS, cols, cone, show_legend)
% Draws the Hill-frame (y-x) docking scene into axes AX: shaded LoS cone,
% the three trajectories, initial condition and docking port.
    hold(ax, 'on'); grid(ax, 'on'); axis(ax, 'equal');

    % LoS approach cone (z = 0 cross-section) as a shaded wedge
    fill(ax, [cone.y, fliplr(cone.y)], [cone.edge, fliplr(-cone.edge)], ...
        [0.98 0.96 0.75], 'FaceAlpha', 0.20, 'EdgeColor', [0.87 0.83 0.58], ...
        'LineWidth', 0.5, ...
        'DisplayName', sprintf('LoS-cone ($\\pm%d^\\circ$)', cone.alpha_deg));

    % Presentation-facing run labels (thesis/slide terminology, not the
    % code's internal run names), matching the Operation 1/2 convention.
    styles = {'-', '-', '--'};
    names  = {'Std. MPC', 'EMPC + term. cost', 'EMPC + term. eq.'};
    for i = 1:3
        plot(ax, R{i}.X(2,:), R{i}.X(1,:), styles{i}, 'Color', cols{i}, ...
            'LineWidth', 1.4, 'DisplayName', names{i});
    end

    plot(ax, x0(2), x0(1), '.k', 'MarkerSize', 18, 'DisplayName', 'IC');
    plot(ax, 0, 0, '+', 'Color', [0.2 0.6 0.2], 'MarkerSize', 10, ...
        'LineWidth', 1.5, 'DisplayName', 'target');

    xlabel(ax, 'Along-track $y$ [m]', 'Interpreter', 'latex');
    ylabel(ax, 'Radial $x$ [m]',      'Interpreter', 'latex');
    if show_legend
        legend(ax, 'Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 8);
    end
end

function s = make_run_entry(res, run_label, cost_label, has_terminal_eq, ...
                             N_used, dt_val, x0_val, x_dock_val, u_max_val)
    s.metadata = struct( ...
        'run',         run_label, ...
        'cost',        cost_label, ...
        'terminal_eq', has_terminal_eq, ...
        'N',           N_used, ...
        'dt_s',        dt_val, ...
        'x0',          x0_val', ...
        'x_dock',      x_dock_val', ...
        'u_max',       u_max_val);
    s.results = struct( ...
        'total_dv_ms',            res.dv, ...
        'n_steps_to_convergence', res.conv_step, ...
        'terminal_distance_m',    res.term_dist, ...
        'passive_safe_fraction',  res.ps_frac, ...
        'n_ps_rows_avg',          mean(res.n_ps_h), ...
        'n_ps_rows_max',          max([res.n_ps_h, 0]));
    s.turnpike = struct( ...
        'turnpike_steps',    res.turnpike_steps, ...
        'turnpike_fraction', res.turnpike_fraction, ...
        'tp_label',          res.tp_label, ...
        'u_norm_sequence',   res.u_norm_seq);
    s.t   = res.t_vec;
    s.X   = res.X';
    s.U   = res.U';
    s.ef  = res.ef;
end
