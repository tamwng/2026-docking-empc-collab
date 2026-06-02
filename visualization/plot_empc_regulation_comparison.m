% plot_empc_regulation_comparison.m
% Presentation figure: Hill frame trajectory comparison, Std MPC vs EMPC.
%
% Reads from exports/scenarios/sim_empc_regulation.json.
% Requires that file to be produced with Run 2 parameters:
%   dt=60 s, N=20, IC=[500;10000;200;0;0;0], no safety constraints.
%   Run sim_empc_regulation.m with those settings first.

clear; clc;
addpath('tools/matlab2tikz-master/src');

%% Load
raw    = jsondecode(fileread('exports/scenarios/sim_empc_regulation.json'));
X_std  = raw.X_std;              % [n_steps x 6]: col1=x(radial), col2=y(along-track)
X_eco  = raw.X_eco;
U_std  = raw.U_std;              % [n_steps x 3]
U_eco  = raw.U_eco;
t_min  = raw.t / 60;
dv_std = raw.metadata.dv_std_ms;
dv_eco = raw.metadata.dv_eco_ms;
u_max  = raw.metadata.u_max;

% Maximum achievable norm when all three axes saturate simultaneously
u_max_norm = sqrt(3) * u_max;   % [m/s^2]

n = size(X_std, 1);

%% Colors (colorblind-safe)
c_std = [0.00 0.45 0.70];   % blue   — Std MPC
c_eco = [0.85 0.33 0.10];   % orange — EMPC

%% Figure 1 — Hill frame trajectory
fig = figure('Name', 'Regulation — Hill Frame', ...
             'Units', 'centimeters', 'Position', [2 2 14 9]);
ax = axes(fig);
hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
axis(ax, 'equal');

plot(ax, X_std(:,2), X_std(:,1), 'Color', c_std, 'LineWidth', 1.8, ...
     'DisplayName', sprintf('Std MPC  \\Deltav = %.1f m/s', dv_std));
plot(ax, X_eco(:,2), X_eco(:,1), 'Color', c_eco, 'LineWidth', 1.8, ...
     'DisplayName', sprintf('EMPC  \\Deltav = %.1f m/s', dv_eco));

% IC and target markers
plot(ax, X_std(1,2), X_std(1,1), 'o', ...
     'Color', [0.4 0.4 0.4], 'MarkerFaceColor', [0.4 0.4 0.4], 'MarkerSize', 5, ...
     'HandleVisibility', 'off');
plot(ax, 0, 0, '+k', 'MarkerSize', 9, 'LineWidth', 1.8, 'HandleVisibility', 'off');

legend(ax, 'Location', 'northeast', 'FontSize', 9);
xlabel(ax, 'Along-track $y$ [m]', 'Interpreter', 'latex', 'FontSize', 10);
ylabel(ax, 'Radial $x$ [m]',      'Interpreter', 'latex', 'FontSize', 10);

%% TikZ export — trajectory
[~,~] = mkdir('results/figures/empc');
fig_out = figure('Visible', 'off');
ax_out  = copyobj(ax, fig_out);
ax_out.Position = [0.15 0.15 0.75 0.75];
matlab2tikz('results/figures/empc/regulation_trajectory_comparison.tikz', ...
    'figurehandle',    fig_out, ...
    'width',           '\figurewidth', ...
    'height',          '\figureheight', ...
    'showInfo',        false, ...
    'checkForUpdates', false);
close(fig_out);
fprintf('TikZ: results/figures/empc/regulation_trajectory_comparison.tikz\n');
fprintf('Std MPC dv = %.2f m/s,  EMPC dv = %.2f m/s\n', dv_std, dv_eco);

%% Figure 2 — Thrust magnitude over time
fig2 = figure('Name', 'Regulation — Thrust Profile', ...
              'Units', 'centimeters', 'Position', [18 2 14 9]);
ax2 = axes(fig2);
hold(ax2, 'on'); grid(ax2, 'on'); box(ax2, 'on');

plot(ax2, t_min, vecnorm(U_std, 2, 2) * 1000, 'Color', c_std, 'LineWidth', 1.6, ...
     'DisplayName', sprintf('Std MPC  \\Deltav = %.1f m/s', dv_std));
plot(ax2, t_min, vecnorm(U_eco, 2, 2) * 1000, 'Color', c_eco, 'LineWidth', 1.6, ...
     'DisplayName', sprintf('EMPC  \\Deltav = %.1f m/s', dv_eco));

% Dashed line: max achievable norm (all 3 axes saturated simultaneously)
yline(ax2, u_max_norm * 1000, '--k', 'LineWidth', 0.8, 'HandleVisibility', 'off');

legend(ax2, 'Location', 'northeast', 'FontSize', 9);
xlabel(ax2, 'Time [min]',                          'FontSize', 10);
ylabel(ax2, '$\|u_k\|$ [mm/s$^2$]', 'Interpreter', 'latex', 'FontSize', 10);

%% TikZ export — thrust profile
fig_out2 = figure('Visible', 'off');
ax_out2  = copyobj(ax2, fig_out2);
ax_out2.Position = [0.15 0.15 0.75 0.75];
matlab2tikz('results/figures/empc/regulation_thrust_profile.tikz', ...
    'figurehandle',    fig_out2, ...
    'width',           '\figurewidth', ...
    'height',          '\figureheight', ...
    'showInfo',        false, ...
    'checkForUpdates', false);
close(fig_out2);
fprintf('TikZ: results/figures/empc/regulation_thrust_profile.tikz\n');
