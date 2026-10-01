% sim_hillframe_nmc.m
%
% Operation 2 (NMC fly-around): Hill-frame trajectory comparison at three
% representative eccentricities (e = 0, 0.1, 0.3), same controller/terminal
% ingredient/scenario as sim_nmc_robustness_empc.m -- just three fixed
% points instead of the full 16-point sweep.
%
% At e=0.3 the controller hits the QP-infeasibility cliff (Section
% "Operation 2: NMC fly-around"): quadprog falls back to zero thrust on
% infeasible steps rather than stopping, so the raw trajectory keeps
% propagating under natural (uncontrolled) drift out to ~hundreds of km,
% which would swamp the shared axis scale. The plotted line is therefore
% truncated at the first infeasible step -- the point where the controller
% actually loses the orbit -- while the full untruncated trajectory is
% still saved to JSON.
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

%% S2 discretisation + CW-NMC orbit Pi* (identical setup to
%% sim_nmc_robustness_empc.m)
P  = 92;
dt = T / P;

Ac = [0 0 0 1 0 0;
      0 0 0 0 1 0;
      0 0 0 0 0 1;
      3*n^2 0 0 0 2*n 0;
      0 0 0 -2*n 0 0;
      0 0 -n^2 0 0 0];
Bc = [zeros(3,3); eye(3)];
[A_d, B_d] = discretize(Ac, Bc, dt);

b_nmc = 75;
x_star = zeros(6, P);
x_star(:,1) = [b_nmc; 0; 0; 0; -2*n*b_nmc; 0];
for k = 1:P-1, x_star(:,k+1) = A_d * x_star(:,k); end
assert(norm(A_d*x_star(:,end) - x_star(:,1)) < 1e-6, 'Pi* does not close under A_d');

u_max = 1e-2;  N = P;  n_sim = 4*P;  n_x = 6;  n_u = 3;
Sx = zeros(n_x*N, n_x);  Su = zeros(n_x*N, n_u*N);
Ad_pow = eye(n_x);
for i = 1:N
    Ad_pow = Ad_pow * A_d;
    Sx((i-1)*n_x+1:i*n_x, :) = Ad_pow;
    for j = 1:i
        Su((i-1)*n_x+1:i*n_x, (j-1)*n_u+1:j*n_u) = A_d^(i-j) * B_d;
    end
end
pos_idx = zeros(1, 3*N);
for ki = 1:N, pos_idx((ki-1)*3+1:ki*3) = (ki-1)*6+1:(ki-1)*6+3; end

matrices = struct('A_d', A_d, 'B_d', B_d, 'Sx', Sx, 'Su', Su, ...
    'Sx_pos_all', Sx(pos_idx,:), 'Su_pos_all', Su(pos_idx,:), ...
    'H', (2*eye(3*N)+2*eye(3*N))/2, 'A_u', [eye(3*N); -eye(3*N)], ...
    'b_u', u_max*ones(6*N,1));

x0 = [0; 300; 0; 0; 0; 0];
params = struct('x0', x0, 'rho_min', 50, 'rho_max', 200, 'u_max', u_max, ...
                'N', N, 'P', P, 'n', n, 'dt', dt, 'n_sim', n_sim, 'verbose', false);

run_cfg = struct('use_nmc_manifold', true, 'use_periodic_terminal', true, ...
                 'use_band', false, 'x_star', x_star);

p_off = struct('mu', mu, 'use_zonal', false, 'use_drag', false);
plant = struct('a', a, 'i', deg2rad(51.6), 'p_t', p_off, 'p_c', p_off);

%% Run at the three representative eccentricities
e_vals  = [0, 0.1, 0.3];
labels  = arrayfun(@(e) sprintf('$e=%.2g$', e), e_vals, 'UniformOutput', false);
colors  = [0.0000 0.4470 0.7410;   % e=0    MATLAB default blue
           0.8500 0.3250 0.0980;   % e=0.1  MATLAB default orange
           0.6350 0.0780 0.1840];  % e=0.3  MATLAB default dark red
traj_full  = cell(1, numel(e_vals));   % untruncated, for JSON
traj_plot  = cell(1, numel(e_vals));   % truncated for the figure
cut_idx    = zeros(1, numel(e_vals));

% Range-based cutoff, not literal QP infeasibility: by the time e=0.3 first
% goes infeasible (empirically step ~194 of 368) the trajectory has already
% ballooned to >10 km, which swamps the shared axis and hides the Pi*
% ellipse and the e=0/e=0.1 curves entirely (verified by rendering it).
% range_cut sits just above e=0.1's own max excursion (587 m, from the
% robustness sweep), so it only ever triggers on the e=0.3 divergence.
range_cut = 800;

for ie = 1:numel(e_vals)
    plant.e = e_vals(ie);
    res = run_flyaround_truth(run_cfg, matrices, params, plant);
    traj_full{ie} = res.x_log;   % [6 x n_sim+1]

    first_infeas = find(res.exitflag_log <= 0, 1);
    first_far    = find(res.range_log > range_cut, 1);
    cut_idx(ie)  = min([size(res.x_log, 2), first_infeas, first_far], [], 'omitnan');
    traj_plot{ie} = res.x_log(:, 1:cut_idx(ie));

    fprintf('e=%.2g: n_infeasible=%d, first_infeasible_step=%s, first_range>%dm=%s, cut at %d, max_range=%.0f m\n', ...
        e_vals(ie), res.n_infeasible, mat2str(first_infeas), range_cut, mat2str(first_far), ...
        cut_idx(ie), max(res.range_log));
end

%% Export JSON (full untruncated trajectory, reproducible without resimulating)
data = struct('metadata', struct( ...
    'study',      'hillframe_nmc', ...
    'controller', 'economic fly-around EMPC (S2 Run 2: pure fuel + periodic terminal eq x(N)=Pi*)', ...
    'plant',      'nonlinear Keplerian truth (perturbations off)', ...
    'scenario',   sprintf('NMC b=%.0f m, V-bar 300 m IC, dt=T/%d, N=P=%d, %d orbits', b_nmc, P, P, n_sim/P), ...
    'e_vals',     e_vals, ...
    'plot_truncated_at_step', cut_idx, ...
    'plot_truncation_rule', sprintf('first of: QP infeasible, or range > %d m', range_cut), ...
    'a_km',       a/1000, ...
    'date',       datestr(now, 'yyyy-mm-dd HH:MM:SS')));
for ie = 1:numel(e_vals)
    fname = matlab.lang.makeValidName(sprintf('e_%g', e_vals(ie)));
    data.(fname) = struct('x_radial_m', traj_full{ie}(1,:), 'y_alongtrack_m', traj_full{ie}(2,:));
end
fid = fopen(fullfile(out_json, 'hillframe_nmc.json'), 'w');
fprintf(fid, '%s', jsonencode(data, 'PrettyPrint', true));
fclose(fid);

%% Figure -- Hill-frame trajectory overlay (along-track y horizontal, radial
% x vertical, matching sim_scenario1_docking.m / sim_scenario4_closing.m).
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

fig = figure('Name', 'Hill-frame trajectory vs eccentricity - NMC fly-around');
ax  = axes(fig);
hold(ax, 'on'); grid(ax, 'on'); axis(ax, 'equal');

x_star_cl = [x_star, x_star(:,1)];   % close the curve for plotting
plot(ax, x_star_cl(2,:), x_star_cl(1,:), '--', 'Color', [0.4 0.4 0.4], ...
    'LineWidth', 1.2, 'DisplayName', sprintf('$\\Pi^*$ ($b=%d$ m)', b_nmc));

for ie = 1:numel(e_vals)
    X = traj_plot{ie};
    plot(ax, X(2,:), X(1,:), '-', 'Color', colors(ie,:), 'LineWidth', 1.4, ...
        'DisplayName', labels{ie});
    plot(ax, X(2,1), X(1,1), 'o', 'Color', colors(ie,:), ...
        'MarkerFaceColor', colors(ie,:), 'MarkerSize', 6, 'HandleVisibility', 'off');
    plot(ax, X(2,end), X(1,end), '^', 'Color', colors(ie,:), ...
        'MarkerFaceColor', colors(ie,:), 'MarkerSize', 6, 'HandleVisibility', 'off');
end
plot(ax, 0, 0, '.k', 'MarkerSize', 14, 'DisplayName', 'Target');

xlabel(ax, 'Along-track $y$ [m]');
ylabel(ax, 'Radial $x$ [m]');
legend(ax, 'Location', 'northeast');

tikz_path = fullfile(out_fig, 'hillframe_nmc.tikz');
matlab2tikz(tikz_path, 'figurehandle', fig, 'showInfo', false);

% legend(...,'FontSize',...) is silently dropped by this matlab2tikz
% version (same issue as in the other hillframe/robustness scripts).
txt = fileread(tikz_path);
txt = regexprep(txt, 'legend style=\{', 'legend style={font=\\large, ', 'once');
fid = fopen(tikz_path, 'w'); fprintf(fid, '%s', txt); fclose(fid);

fprintf('Exported JSON and figure.\n');
