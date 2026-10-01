% plot_rpo.m
% Dedicated figure export script for sim_rpo results.
% Run AFTER sim_rpo.m has been executed.
% Exports TikZ .tex files to results/figures/rpo/ (thesis figure convention).

clear; clc;
%%Paths 
% Get the directory where plot_rpo.m lives, then navigate from there
base_dir = fileparts(mfilename('fullpath'));   % .../2026-docking-empc-collab/visualization
repo_dir = fullfile(base_dir, '..');           % .../2026-docking-empc-collab

addpath(fullfile(repo_dir, 'tools', 'matlab2tikz-master', 'src'));

output_dir = fullfile(repo_dir, 'results', 'figures', 'rpo');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

%% Load simulation output 
load(fullfile(repo_dir, 'exports', 'scenarios', 'sim_rpo.mat'));
%% Global LaTeX-consistent style
set(groot, 'defaultTextInterpreter',          'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter',        'latex');
set(groot, 'defaultAxesFontSize',             10);
set(groot, 'defaultLineLineWidth',            1.5);

figW = '\figwidth';    % resolved in LaTeX via \newlength
figH = '\figheight';

%% ── Figure 1: 2D Phase Plot (bounded orbit) ─────────────────────────────
fig1 = figure('Units','centimeters','Position',[0 0 12 8]);
plot(Rho(2,:)/1e3, Rho(1,:)/1e3, 'b');
hold on;
plot(Rho(2,1)/1e3, Rho(1,1)/1e3, 'go', 'MarkerSize', 6, ...
     'DisplayName', 'Start');
xlabel('Along-track $y$ [km]');
ylabel('Radial $x$ [km]');
title('Bounded Relative Orbit -- CWH Model');
legend('Trajectory', 'Start', 'Location', 'best');
grid on; axis equal;

matlab2tikz(fullfile(output_dir, 'rpo_bounded_2d.tex'), ...
    'width', figW, 'height', figH, 'showInfo', false, ...
    'extraAxisOptions', {'ylabel style={font=\small}', ...
                         'xlabel style={font=\small}'});
close(fig1);

%% ── Figure 2: Secular drift (drift-free condition OFF) ───────────────────
% Recompute quickly with y_dot0 = 0
x0 = -1000e3; y0 = 0; z0 = 500e3;
Rho0_drift = [x0; y0; z0; 0; 0; 0];   % y_dot0 = 0, no drift-free condition

% Integrate inline
n_steps_drift = length(t);
Rho_drift = zeros(6, n_steps_drift);
Rho_drift(:,1) = Rho0_drift;
dt = t(2) - t(1);

for k = 1:n_steps_drift-1
    f  = @(rho) clohessy_wiltshire(rho, n);
    xk = Rho_drift(:,k);
    k1 = f(xk);
    k2 = f(xk + 0.5*dt*k1);
    k3 = f(xk + 0.5*dt*k2);
    k4 = f(xk + dt*k3);
    Rho_drift(:,k+1) = xk + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
end

fig2 = figure('Units','centimeters','Position',[0 0 12 8]);
plot(Rho_drift(2,:)/1e3, Rho_drift(1,:)/1e3, 'r');
hold on;
plot(Rho_drift(2,1)/1e3, Rho_drift(1,1)/1e3, 'go', 'MarkerSize', 6);
xlabel('Along-track $y$ [km]');
ylabel('Radial $x$ [km]');
title('Secular Along-Track Drift -- CWH Model ($\dot{y}_0 = 0$)');
grid on;

matlab2tikz(fullfile(output_dir, 'rpo_drift_2d.tex'), ...
    'width', figW, 'height', figH, 'showInfo', false);
close(fig2);

%% ── Figure 3: Time histories ─────────────────────────────────────────────
fig3 = figure('Units','centimeters','Position',[0 0 14 10]);

subplot(3,1,1);
plot(t/60, Rho(1,:)/1e3, 'b');
ylabel('$x$ [km]'); grid on;

subplot(3,1,2);
plot(t/60, Rho(2,:)/1e3, 'b');
ylabel('$y$ [km]'); grid on;

subplot(3,1,3);
plot(t/60, Rho(3,:)/1e3, 'b');
ylabel('$z$ [km]'); grid on;
xlabel('Time [min]');

matlab2tikz(fullfile(output_dir, 'rpo_time_history.tex'), ...
    'width', figW, 'height', figH, 'showInfo', false);
close(fig3);

disp('All figures exported to results/figures/rpo/');