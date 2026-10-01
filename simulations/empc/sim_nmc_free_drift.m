% sim_nmc_free_drift.m
% Single Hill-frame plot: NMC 2:1 ellipse (rho = 75 m) vs. secular drift
% from the same initial position, with and without the drift-free condition
%   ydot_0 = -2 n x_0.
%
% Starting point: x0 = -rho, y0 = 0  (bottom of the NMC ellipse, phi = pi)
%   => with drift-free condition  ydot_0 = 2 n rho  -> closed orbit
%   => with ydot_0 = 0 (condition violated)          -> secular along-track drift

clear; clc;
script_dir = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(script_dir, '..', '..', 'tools', 'matlab2tikz-master', 'src')));
constants;   % loads n, T

rho = 75;          % [m]  radial semi-axis of NMC orbit
x0  = -rho;        % [m]  starting radial position (bottom of ellipse)
y0  =  0;          % [m]  starting along-track position

% Time grids
dt     = T / 500;
n_nmc  = round(T / dt);         % one full orbit (to show closure)
n_free = round(2 * T / dt);      % 2 full orbits of free drift
t_nmc  = (0 : n_nmc-1)  * dt;
t_free = (0 : n_free-1) * dt;

%% Analytical CWH free-drift solution
%   x(t) = (4-3c) x0 + (s/n) xd0 + 2(1-c)/n yd0
%   y(t) = 6(s-nt) x0 + y0 - 2(1-c)/n xd0 + (4s/n - 3t) yd0
%   c = cos(nt),  s = sin(nt)

cwh = @(t, xd0, yd0) deal( ...
    (4 - 3*cos(n*t))*x0 + sin(n*t)/n * xd0 + 2*(1-cos(n*t))/n * yd0, ...
    6*(sin(n*t) - n*t)*x0 + y0 - 2*(1-cos(n*t))/n * xd0 + (4*sin(n*t)/n - 3*t)*yd0 );

% Drift-free: ydot_0 = -2 n x0 = 2 n rho
[x_nmc,  y_nmc]  = cwh(t_nmc,  0, -2*n*x0);

% Free drift: ydot_0 = 0
[x_free, y_free] = cwh(t_free, 0, 0);

% Drift-free case over the same 2-orbit window, for the time-history plot
[x_nmc_th, y_nmc_th] = cwh(t_free, 0, -2*n*x0);

%% Plot
c_nmc  = [0.00 0.45 0.70];   % blue: bounded orbit
c_free = [0.85 0.33 0.10];   % red: secular drift
c_ic   = [0.18 0.49 0.20];   % green: start marker

figure('Name', 'NMC vs free drift', ...
       'Units', 'centimeters', 'Position', [2 4 26 11]);

%% Left: NMC orbit (axis equal, tight scale)
ax1 = subplot(1, 2, 1);
hold(ax1, 'on'); grid(ax1, 'on'); axis(ax1, 'equal');

plot(ax1, y_nmc, x_nmc, 'Color', c_nmc, 'LineWidth', 1.8, ...
     'DisplayName', ['NMC orbit, $\dot{y}_0 = -2nx_0$']);
plot(ax1, y0, x0, 'o', 'Color', c_ic, 'MarkerFaceColor', c_ic, ...
     'MarkerSize', 7, 'DisplayName', 'Initial state');
plot(ax1, 0, 0, '+k', 'MarkerSize', 9, 'LineWidth', 1.8, ...
     'DisplayName', 'Target');

xlabel(ax1, 'Along-track $y$ [m]', 'Interpreter', 'latex', 'FontSize', 11);
ylabel(ax1, 'Radial $x$ [m]',      'Interpreter', 'latex', 'FontSize', 11);
legend(ax1, 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 9);

%% Right: free drift (auto scale)
ax2 = subplot(1, 2, 2);
hold(ax2, 'on'); grid(ax2, 'on');

plot(ax2, y_free, x_free, 'Color', c_free, 'LineWidth', 1.8, ...
     'DisplayName', 'Free drift, $\dot{y}_0 = 0$');
plot(ax2, y0, x0, 'o', 'Color', c_ic, 'MarkerFaceColor', c_ic, ...
     'MarkerSize', 7, 'DisplayName', 'Initial state');
plot(ax2, 0, 0, '+k', 'MarkerSize', 9, 'LineWidth', 1.8, ...
     'DisplayName', 'Target');

xlabel(ax2, 'Along-track $y$ [m]', 'Interpreter', 'latex', 'FontSize', 11);
ylabel(ax2, 'Radial $x$ [m]',      'Interpreter', 'latex', 'FontSize', 11);
legend(ax2, 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 9);

%% TikZ export
fig = gcf;
out_dir = fullfile(script_dir, '..', '..', 'results', 'figures', 'nmc_free_drift');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

if exist('matlab2tikz', 'file')
    matlab2tikz(fullfile(out_dir, 'nmc_free_drift.tikz'), ...
        'figurehandle', fig, 'showInfo', false, ...
        'width',  '\linewidth', ...
        'height', '0.45\linewidth');
    fprintf('TikZ exported to %s\n', out_dir);
else
    fprintf('TikZ export skipped, run setup.m first.\n');
end
exportgraphics(fig, fullfile(out_dir, 'nmc_free_drift.png'), 'Resolution', 300);
fprintf('PNG exported to %s\n', out_dir);

%% Time-history figure (separate from the Hill-frame plots)
t_min = t_free / 60;   % [min]

fig2 = figure('Name', 'NMC vs free drift - time histories', ...
              'Units', 'centimeters', 'Position', [2 4 20 14]);

ax3 = subplot(2, 1, 1);
hold(ax3, 'on'); grid(ax3, 'on');
plot(ax3, t_min, x_nmc_th, 'Color', c_nmc,  'LineWidth', 1.8, ...
     'DisplayName', ['NMC orbit, $\dot{y}_0 = -2nx_0$']);
plot(ax3, t_min, x_free,   'Color', c_free, 'LineWidth', 1.8, ...
     'DisplayName', 'Free drift, $\dot{y}_0 = 0$');
ylabel(ax3, 'Radial $x$ [m]', 'Interpreter', 'latex', 'FontSize', 11);
legend(ax3, 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 9);

ax4 = subplot(2, 1, 2);
hold(ax4, 'on'); grid(ax4, 'on');
plot(ax4, t_min, y_nmc_th, 'Color', c_nmc,  'LineWidth', 1.8, ...
     'DisplayName', ['NMC orbit, $\dot{y}_0 = -2nx_0$']);
plot(ax4, t_min, y_free,   'Color', c_free, 'LineWidth', 1.8, ...
     'DisplayName', 'Free drift, $\dot{y}_0 = 0$');
xlabel(ax4, 'Time [min]', 'Interpreter', 'latex', 'FontSize', 11);
ylabel(ax4, 'Along-track $y$ [m]', 'Interpreter', 'latex', 'FontSize', 11);
legend(ax4, 'Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 9);

if exist('matlab2tikz', 'file')
    matlab2tikz(fullfile(out_dir, 'nmc_free_drift_timehist.tikz'), ...
        'figurehandle', fig2, 'showInfo', false, ...
        'width',  '\linewidth', ...
        'height', '0.7\linewidth');
    fprintf('TikZ exported to %s\n', out_dir);
else
    fprintf('TikZ export skipped, run setup.m first.\n');
end
exportgraphics(fig2, fullfile(out_dir, 'nmc_free_drift_timehist.png'), 'Resolution', 300);
fprintf('PNG exported to %s\n', out_dir);