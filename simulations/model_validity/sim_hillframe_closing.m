% sim_hillframe_closing.m
%
% Operation 1 (far-range closing): Hill-frame trajectory comparison at three
% representative eccentricities (e = 0, 0.1, 0.3), same controller/cost/
% constraints/scenario as sim_closing_robustness.m -- just three fixed
% points instead of the full 16-point sweep.
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

%% Closing scenario (identical setup to sim_closing_robustness.m)
dt     = 120;  N = 20;
x0     = [0; 10000; 0; 0; 0; 0];
x_sw   = [0;   300; 0; 0; 0; 0];
u_max  = 1e-2;

Ac = [0 0 0 1 0 0;
      0 0 0 0 1 0;
      0 0 0 0 0 1;
      3*n^2 0 0 0 2*n 0;
      0 0 0 -2*n 0 0;
      0 0 -n^2 0 0 0];
Bc = [zeros(3,3); eye(3)];
[Ad, Bd] = discretize(Ac, Bc, dt);
Q_dare = diag([1e-2 1e-2 1e-2 1 1 1]);  R_dare = diag([1e4 1e4 1e4]);
[~, P_clf] = lqr_controller(Ad, Bd, Q_dare, R_dare);

cost = struct('Q', zeros(6), 'R', eye(3), 'P', P_clf);
con  = struct('u_max', u_max, 'y_min_active', false);
con.los_cone.active = false;

params = struct('x0', x0, 'x_switch', x_sw, 'dt', dt, 'N', N, 'u_max', u_max, ...
                't_final', 5*T, 'conv_tol', 50, 'verbose', false);

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
    res = run_closing_truth(cost, con, params, plant);
    traj{ie} = res.X_hill;   % [6 x n], rows: x,y,z,vx,vy,vz (Hill frame)
    fprintf('e=%.2g: reached=%d, term_dist=%.1f m, dv=%.4f m/s, n_steps=%d\n', ...
        e_vals(ie), res.reached_ball, res.term_dist, res.dv, size(res.X_hill, 2));
end

%% Export JSON (raw trajectory data, reproducible without resimulating)
data = struct('metadata', struct( ...
    'study',      'hillframe_closing', ...
    'controller', 'CW-EMPC pure-fuel + CLF terminal cost', ...
    'plant',      'nonlinear Keplerian truth (perturbations off)', ...
    'scenario',   '10 km -> 300 m closing, N=20, dt=120 s', ...
    'e_vals',     e_vals, ...
    'a_km',       a/1000, ...
    'date',       datestr(now, 'yyyy-mm-dd HH:MM:SS')));
for ie = 1:numel(e_vals)
    fname = matlab.lang.makeValidName(sprintf('e_%g', e_vals(ie)));
    data.(fname) = struct('x_radial_m', traj{ie}(1,:), 'y_alongtrack_m', traj{ie}(2,:));
end
fid = fopen(fullfile(out_json, 'hillframe_closing.json'), 'w');
fprintf(fid, '%s', jsonencode(data, 'PrettyPrint', true));
fclose(fid);

%% Figure -- Hill-frame trajectory overlay (along-track y horizontal, radial
% x vertical, matching sim_scenario1_docking.m / sim_scenario4_closing.m).
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

fig = figure('Name', 'Hill-frame trajectory vs eccentricity - Closing');
ax  = axes(fig);
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
plot(ax, x_sw(2), x_sw(1), 's', 'Color', [0.2 0.6 0.2], ...
    'MarkerFaceColor', [0.2 0.6 0.2], 'MarkerSize', 7, 'DisplayName', 'Handover point');
plot(ax, 0, 0, '+k', 'MarkerSize', 9, 'LineWidth', 1.5, 'DisplayName', 'Target');

xlabel(ax, 'Along-track $y$ [m]');
ylabel(ax, 'Radial $x$ [m]');
% Explicit, sparse ticks -- pgfplots' auto ticking on the 0-10000 m
% along-track range was falling back to a x10^4 scale factor with
% fractional labels (0.05, 0.1, ...), which is what needed cleaning up.
xticks(ax, 0:2000:10000);
yticks(ax, -1000:1000:3000);
legend(ax, 'Location', 'northeast');

tikz_path = fullfile(out_fig, 'hillframe_closing.tikz');
matlab2tikz(tikz_path, 'figurehandle', fig, 'showInfo', false);

% legend(...,'FontSize',...) is silently dropped by this matlab2tikz
% version (same issue as in the other hillframe/robustness scripts).
txt = fileread(tikz_path);
txt = regexprep(txt, 'legend style=\{', 'legend style={font=\\large, ', 'once');
fid = fopen(tikz_path, 'w'); fprintf(fid, '%s', txt); fclose(fid);

fprintf('Exported JSON and figure.\n');
