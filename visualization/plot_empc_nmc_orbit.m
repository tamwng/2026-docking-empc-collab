% plot_empc_nmc_orbit.m
% Presentation figure: NMC orbit acquisition and maintenance in the Hill frame.
%
% Reads from exports/scenarios/sim_empc_nmc.json (confirmed Run 1 data).
% Shows: hold point → injection arc → NMC ellipse for 3 orbits.

clear; clc;
addpath('tools/matlab2tikz-master/src');

%% Load
raw   = jsondecode(fileread('exports/scenarios/sim_empc_nmc.json'));
X     = raw.X;      % [n_steps x 6]: col1=x(radial), col2=y(along-track)
U     = raw.U;      % [n_steps x 3]
meta  = raw.metadata;

dt    = meta.dt;
T_orb = meta.period_min * 60;   % orbital period [s]
rho   = meta.rho_m;
phi_0 = meta.phi_0_rad;
y_c   = meta.y_c_m;
u_max = meta.u_max;

n_steps = size(X, 1);
t_min   = (0:n_steps-1) * dt / 60;

%% Injection / maintenance Δv
n_inject  = round(0.25 * T_orb / dt);
dv_cumsum = cumsum(vecnorm(U, 2, 2) * dt);
dv_inject = dv_cumsum(n_inject);
dv_maint  = dv_cumsum(end) - dv_cumsum(n_inject);

fprintf('Injection  dv: %.4f m/s  (first %d steps, %.0f min)\n', ...
        dv_inject, n_inject, n_inject * dt / 60);
fprintf('Maintenance dv: %.4f m/s\n', dv_maint);

%% Ideal NMC ellipse (one clean orbit for reference)
phi_ell = phi_0 + linspace(0, 2*pi, 500);
x_ell   =  rho * cos(phi_ell);
y_ell   =  y_c - 2*rho * sin(phi_ell);

%% Colors
c_arc  = [0.00 0.45 0.70];   % blue   — injection arc
c_nmc  = [0.85 0.33 0.10];   % orange — maintained NMC orbit
c_hold = [0.40 0.40 0.40];   % grey   — hold point

%% Figure 1 — Hill frame
fig = figure('Name', 'EMPC NMC — Orbit Acquisition', ...
             'Units', 'centimeters', 'Position', [2 2 13 11]);
ax = axes(fig);
hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
axis(ax, 'equal');

% Ideal ellipse (dashed background reference)
plot(ax, y_ell, x_ell, '--', 'Color', [0.75 0.75 0.75], 'LineWidth', 1.0, ...
     'HandleVisibility', 'off');

% Maintenance: actual trajectory after injection
plot(ax, X(n_inject:end, 2), X(n_inject:end, 1), ...
     'Color', c_nmc, 'LineWidth', 1.6, ...
     'DisplayName', sprintf('Maintenance  \\Deltav = %.3f m/s / 3 orbits', dv_maint));

% Injection arc
plot(ax, X(1:n_inject, 2), X(1:n_inject, 1), ...
     'Color', c_arc, 'LineWidth', 2.2, ...
     'DisplayName', sprintf('Injection  \\Deltav = %.3f m/s', dv_inject));

% Hold point marker
plot(ax, X(1,2), X(1,1), 'o', ...
     'Color', c_hold, 'MarkerFaceColor', c_hold, 'MarkerSize', 7, ...
     'HandleVisibility', 'off');

% Target (origin)
plot(ax, 0, 0, '+k', 'MarkerSize', 9, 'LineWidth', 1.8, ...
     'HandleVisibility', 'off');

% Hold point label inside the ellipse
text(ax, X(1,2) - 8, X(1,1) - 8, 'hold point', ...
     'Color', c_hold, 'FontSize', 8, ...
     'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');

legend(ax, 'Location', 'northwest', 'FontSize', 8);
xlabel(ax, 'Along-track $y$ [m]', 'Interpreter', 'latex', 'FontSize', 10);
ylabel(ax, 'Radial $x$ [m]',      'Interpreter', 'latex', 'FontSize', 10);

%% TikZ export — orbit
[~,~] = mkdir('results/figures/empc');
fig_out = figure('Visible', 'off');
ax_out  = copyobj(ax, fig_out);
ax_out.Position = [0.15 0.15 0.75 0.75];
matlab2tikz('results/figures/empc/nmc_orbit.tikz', ...
    'figurehandle',    fig_out, ...
    'width',           '\figurewidth', ...
    'height',          '\figureheight', ...
    'showInfo',        false, ...
    'checkForUpdates', false);
close(fig_out);
fprintf('TikZ: results/figures/empc/nmc_orbit.tikz\n');

%% Figure 2 — Control inputs over time
T_marks = (1:3) * T_orb / 60;   % orbit boundary times [min]

c_x = [0.00 0.60 0.90];
c_y = [0.90 0.40 0.00];
c_z = [0.20 0.80 0.20];

fig2 = figure('Name', 'EMPC NMC — Control Inputs', ...
              'Units', 'centimeters', 'Position', [17 2 14 9]);
ax2 = axes(fig2);
hold(ax2, 'on'); grid(ax2, 'on'); box(ax2, 'on');

plot(ax2, t_min, U(:,1)*1000, 'Color', c_x, 'LineWidth', 1.4, 'DisplayName', '$u_x$');
plot(ax2, t_min, U(:,2)*1000, 'Color', c_y, 'LineWidth', 1.4, 'DisplayName', '$u_y$');
plot(ax2, t_min, U(:,3)*1000, 'Color', c_z, 'LineWidth', 1.4, 'DisplayName', '$u_z$');
yline(ax2,  u_max*1000, '--k', 'LineWidth', 0.8, 'HandleVisibility', 'off');
yline(ax2, -u_max*1000, '--k', 'LineWidth', 0.8, 'HandleVisibility', 'off');
for tm = T_marks
    xline(ax2, tm, ':', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8, 'HandleVisibility', 'off');
end
% Zoom y-axis to where the action is — injection peak visible, maintenance detail too
u_peak = max(abs(U(:)), [], 'all') * 1000;
ylim(ax2, [-u_peak * 1.15,  u_peak * 1.15]);

legend(ax2, 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 9);
xlabel(ax2, 'Time [min]',                           'FontSize', 10);
ylabel(ax2, 'Acceleration [mm/s$^2$]', 'Interpreter', 'latex', 'FontSize', 10);

%% TikZ export — inputs
fig_out2 = figure('Visible', 'off');
ax_out2  = copyobj(ax2, fig_out2);
ax_out2.Position = [0.15 0.15 0.75 0.75];
matlab2tikz('results/figures/empc/nmc_inputs.tikz', ...
    'figurehandle',    fig_out2, ...
    'width',           '\figurewidth', ...
    'height',          '\figureheight', ...
    'showInfo',        false, ...
    'checkForUpdates', false);
close(fig_out2);
fprintf('TikZ: results/figures/empc/nmc_inputs.tikz\n');
