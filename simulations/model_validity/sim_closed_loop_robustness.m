% sim_closed_loop_robustness.m
%
% Output: JSON to exports/scenarios/validity_sweep/, figure to
% results/figures/validity_sweep/.

clear; clc; close all;

addpath('src/dynamics'); addpath('src/utils'); addpath('src/control');
constants;

out_json = 'exports/scenarios/validity_sweep';
out_fig  = 'results/figures/validity_sweep';
if ~exist(out_json, 'dir'), mkdir(out_json); end
if ~exist(out_fig,  'dir'), mkdir(out_fig);  end

%% Docking scenario
dt     = 10;   N = 24;
x0     = [0; 150; 0; 0; 0; 0];
x_dock = zeros(6, 1);
u_max  = 1e-2;

% Controller CW model + standard-MPC cost (Q~=0 tracking, DARE terminal)
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

%% Sweep
e_grid = 0 : 0.02 : 0.30;
nE = numel(e_grid);
term_dist = zeros(1, nE);
dv        = zeros(1, nE);
conv_step = nan(1, nE);
los_frac  = zeros(1, nE);
ol_horizon_err = zeros(1, nE);   % open-loop CW error over N*dt, same IC

opts_ode = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);
n_mm = sqrt(mu / a^3);

fprintf('Closed-loop robustness sweep over %d eccentricities...\n', nE);
t0 = tic;
for ie = 1:nE
    plant.e = e_grid(ie);

    % closed loop (truth plant, CW-MPC)
    res = run_docking_truth(cost, con, params, plant);
    term_dist(ie) = res.term_dist;
    dv(ie)        = res.dv;
    if ~isnan(res.conv_step), conv_step(ie) = res.conv_step; end
    los_frac(ie)  = res.los_viol_frac;

    % open-loop CW prediction error over one horizon (no control)
    X_t0 = oe2eci(a, plant.e, plant.i, 0, 0, 0, mu);
    combo = @(t, y) [ truth_propagator_at(t, y(1:6), p_off); ...
                      clohessy_wiltshire(y(7:12), n_mm) ];
    [~, Yc] = ode113(combo, [0, N*dt], [X_t0; x0], opts_ode);
    X_c0 = hill2eci(X_t0, x0);
    combo_c = @(t, y) [ truth_propagator_at(t, y(1:6), p_off); ...
                        truth_propagator_at(t, y(7:12), p_off) ];
    [tt, Yt] = ode113(combo_c, [0, N*dt], [X_t0; X_c0], opts_ode);
    % truth relative at horizon end vs CW propagation of the same IC
    rho_truth_end = eci2hill(Yt(end,1:6).', Yt(end,7:12).');
    rho_cw_end    = Yc(end, 7:12).';
    ol_horizon_err(ie) = norm(rho_truth_end(1:3) - rho_cw_end(1:3));
end
fprintf('Sweep complete in %.1fs.\n', toc(t0));

dv_infl = (dv - dv(1)) / dv(1) * 100;   % Delta-v inflation vs circular [%]

%% Console table
fprintf('\n  %5s  %10s  %8s  %8s  %10s\n', 'e', 'miss[m]', 'dv[m/s]', 'infl[%]', 'OLerr[m]');
for ie = 1:nE
    fprintf('  %5.2f  %10.3f  %8.4f  %8.1f  %10.2f\n', ...
        e_grid(ie), term_dist(ie), dv(ie), dv_infl(ie), ol_horizon_err(ie));
end

%% Export JSON
metadata = struct('study','layer2_closed_loop_robustness', ...
    'controller','CW-based standard MPC (Q~=0, DARE terminal)', ...
    'plant','nonlinear Keplerian truth (perturbations off)', ...
    'scenario','150 m V-bar docking, N=24, dt=10 s', ...
    'a_km', a/1000, 'date', datestr(now,'yyyy-mm-dd HH:MM:SS'));
data = struct('metadata', metadata, 'e_grid', e_grid, ...
    'term_dist_m', term_dist, 'dv_ms', dv, 'dv_inflation_pct', dv_infl, ...
    'conv_step', conv_step, 'los_viol_frac', los_frac, ...
    'open_loop_horizon_err_m', ol_horizon_err);
fid = fopen(fullfile(out_json, 'closed_loop_robustness.json'), 'w');
fprintf(fid, '%s', jsonencode(data, 'PrettyPrint', true));
fclose(fid);

%% Figure
figure('Position', [100 100 1000 420]);

% Panel 1: closed-loop terminal miss + Delta-v vs eccentricity
subplot(1,2,1);
yyaxis left;
plot(e_grid, term_dist, '-o', 'LineWidth', 1.5); hold on;
yline(params.conv_tol, '--', 'conv tol', 'HandleVisibility','off');
ylabel('terminal miss [m]'); ylim([0, max(params.conv_tol*1.5, max(term_dist)*1.2)]);
yyaxis right;
plot(e_grid, dv, '-s', 'LineWidth', 1.5);
ylabel('total \Deltav [m/s]');
xlabel('eccentricity [-]'); grid on;
title('Closed-loop docking (truth plant)');

% Panel 2: the contrast -- open-loop horizon error vs closed-loop miss
subplot(1,2,2);
semilogy(e_grid, ol_horizon_err, '-^', 'LineWidth', 1.5); hold on;
semilogy(e_grid, max(term_dist, 1e-3), '-o', 'LineWidth', 1.5);
grid on; xlabel('eccentricity [-]'); ylabel('position error [m]');
legend('open-loop CW error (1 horizon, no control)', ...
       'closed-loop terminal miss (with MPC)', 'Location','northwest');
title('What feedback buys');

sgtitle('Layer 2 - CW-MPC robustness to eccentricity');
saveas(gcf, fullfile(out_fig, 'closed_loop_robustness.png'));

fprintf('Exported JSON and figure.\n');