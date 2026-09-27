% sim_hillframe_docking.m
%
% Operation 3 (docking): Hill-frame trajectory comparison at three
% representative eccentricities (e = 0, 0.1, 0.3), same controller/cost/
% constraints/scenario as sim_closed_loop_robustness.m -- just three fixed
% points instead of the full 16-point sweep, run for trajectory shape
% rather than summary statistics.
%
% Output: JSON (raw x/y trajectories) to exports/scenarios/validity_sweep/,
% figure to results/figures/validity_sweep/.

clear; clc; close all;

addpath('src/dynamics'); addpath('src/utils'); addpath('src/control');
addpath('tools/matlab2tikz-master/src');
constants;

out_json = 'exports/scenarios/validity_sweep';
out_fig  = 'results/figures/validity_sweep';
if ~exist(out_json, 'dir'), mkdir(out_json); end
if ~exist(out_fig,  'dir'), mkdir(out_fig);  end

%% Docking scenario (identical setup to sim_closed_loop_robustness.m)
dt     = 10;   N = 24;
x0     = [0; 150; 0; 0; 0; 0];
x_dock = zeros(6, 1);
u_max  = 1e-2;

Ac = [0 0 0 1 0 0;
      0 0 0 0 1 0;
      0 0 0 0 0 1;
      3*n^2 0 0 0 2*n 0;
      0 0 0 -2*n 0 0;
      0 0 -n^2 0 0 0];
Bc = [zeros(3,3); eye(3)];
[Ad, Bd] = discretize(Ac, Bc, dt);
Q = diag([1e-2 1e-2 1e-2 1 1 1]);  R = diag([1e4 1e4 1e4]);
[~, P] = lqr_controller(Ad, Bd, Q, R);

cost = struct('Q', Q, 'R', R, 'P', P);
con  = struct('u_max', u_max, 'y_min_active', true);
con.los_cone.active = true;  con.los_cone.half_angle = pi/6;

params = struct('x0', x0, 'x_dock', x_dock, 'dt', dt, 'N', N, 'u_max', u_max, ...
                't_final', 1.5*T, 'conv_tol', 1.0, 'r_KOS', 5, 'verbose', false);

p_off = struct('mu', mu, 'use_zonal', false, 'use_drag', false);
plant = struct('a', a, 'i', deg2rad(51.6), 'p_t', p_off, 'p_c', p_off);

%% Run at the three representative eccentricities
e_vals  = [0, 0.1, 0.3];
labels  = arrayfun(@(e) sprintf('$e=%.2g$', e), e_vals, 'UniformOutput', false);
colors  = [0.0000 0.4470 0.7410;   % e=0    MATLAB default blue
           0.8500 0.3250 0.0980;   % e=0.1  MATLAB default orange
           0.6350 0.0780 0.1840];  % e=0.3  MATLAB default dark red
traj = cell(1, numel(e_vals));

for ie = 1:numel(e_vals)
    plant.e = e_vals(ie);
    res = run_docking_truth(cost, con, params, plant);
    traj{ie} = res.X_hill;   % [6 x n], rows: x,y,z,vx,vy,vz (Hill frame)
    fprintf('e=%.2g: term_dist=%.3f m, dv=%.4f m/s, n_steps=%d\n', ...
        e_vals(ie), res.term_dist, res.dv, size(res.X_hill, 2));
end

%% Export JSON (raw trajectory data, reproducible without resimulating)
data = struct('metadata', struct( ...
    'study',      'hillframe_docking', ...
    'controller', 'CW-based standard MPC (Q~=0, DARE terminal)', ...
    'plant',      'nonlinear Keplerian truth (perturbations off)', ...
    'scenario',   '150 m V-bar docking, N=24, dt=10 s', ...
    'e_vals',     e_vals, ...
    'a_km',       a/1000, ...
    'date',       datestr(now, 'yyyy-mm-dd HH:MM:SS')));
for ie = 1:numel(e_vals)
    fname = matlab.lang.makeValidName(sprintf('e_%g', e_vals(ie)));
    data.(fname) = struct('x_radial_m', traj{ie}(1,:), 'y_alongtrack_m', traj{ie}(2,:));
end
fid = fopen(fullfile(out_json, 'hillframe_docking.json'), 'w');
fprintf(fid, '%s', jsonencode(data, 'PrettyPrint', true));
fclose(fid);

%% Figure -- Hill-frame trajectory overlay (along-track y horizontal, radial
% x vertical, matching sim_scenario1_docking.m / sim_scenario4_closing.m),
% full approach + terminal zoom, same two-panel layout as S1's own
% hill_trajectory.tikz.
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

fig = figure('Name', 'Hill-frame trajectory vs eccentricity - Docking');
fig.Position(3:4) = [1050, 460];

panel_axes = gobjects(1,2);
for p = 1:2
    ax = subplot(1, 2, p);
    panel_axes(p) = ax;
    hold(ax, 'on'); grid(ax, 'on'); axis(ax, 'equal');

    for ie = 1:numel(e_vals)
        X = traj{ie};
        plot(ax, X(2,:), X(1,:), '-', 'Color', colors(ie,:), 'LineWidth', 1.4, ...
            'DisplayName', labels{ie});
        plot(ax, X(2,1), X(1,1), 'o', 'Color', colors(ie,:), ...
            'MarkerFaceColor', colors(ie,:), 'MarkerSize', 6, 'HandleVisibility', 'off');
        plot(ax, X(2,end), X(1,end), '^', 'Color', colors(ie,:), ...
            'MarkerFaceColor', colors(ie,:), 'MarkerSize', 6, 'HandleVisibility', 'off');
    end
    plot(ax, 0, 0, 'pk', 'MarkerSize', 10, 'MarkerFaceColor', 'k', ...
        'DisplayName', 'Docking port (origin)');

    xlabel(ax, 'Along-track $y$ [m]');
    ylabel(ax, 'Radial $x$ [m]');
    if p == 1
        legend(ax, 'Location', 'northwest', 'FontSize', 9);
    else
        % Along-track (horizontal) window slid +1 m so it spans -0.5 to
        % +2.5 m; radial (vertical) stays centered at +-1.5 m.
        xlim(ax, [-0.5, 2.5]); ylim(ax, [-1.5, 1.5]);
    end
end

tikz_path = fullfile(out_fig, 'hillframe_docking.tikz');
matlab2tikz(tikz_path, 'figurehandle', fig, 'showInfo', false);

% legend(...,'FontSize',...) is silently dropped by this matlab2tikz
% version (same issue as in sim_closed_loop_robustness.m) -- patch directly.
txt = fileread(tikz_path);
txt = regexprep(txt, 'legend style=\{', 'legend style={font=\\large, ', 'once');
fid = fopen(tikz_path, 'w'); fprintf(fid, '%s', txt); fclose(fid);

fprintf('Exported JSON and figure.\n');
