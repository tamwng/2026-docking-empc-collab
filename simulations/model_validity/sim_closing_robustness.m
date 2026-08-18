% sim_closing_robustness.m
%
% Controller: CW-based EMPC, pure fuel + CLF terminal cost (Run C style).
% Plant: nonlinear Keplerian truth (perturbations off) at eccentricity e.
%
% Output: JSON + figure in validity_sweep/.

clear; clc; close all;

addpath('src/dynamics'); addpath('src/utils'); addpath('src/control');
constants;

out_json = 'exports/scenarios/validity_sweep';
out_fig  = 'results/figures/validity_sweep';
if ~exist(out_json, 'dir'), mkdir(out_json); end
if ~exist(out_fig,  'dir'), mkdir(out_fig);  end

%% Closing scenario 
dt     = 120;  N = 20;
x0     = [0; 10000; 0; 0; 0; 0];
x_sw   = [0;   300; 0; 0; 0; 0];
u_max  = 1e-2;

% Controller CW model + CLF terminal cost (pure-fuel EMPC)
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

%% Sweep
e_grid = 0 : 0.02 : 0.30;
nE = numel(e_grid);
term_dist = zeros(1, nE);
dv        = zeros(1, nE);
t_hand_h  = nan(1, nE);
reached   = false(1, nE);

fprintf('Closing robustness sweep over %d eccentricities...\n', nE);
t0 = tic;
for ie = 1:nE
    plant.e = e_grid(ie);
    res = run_closing_truth(cost, con, params, plant);
    term_dist(ie) = res.term_dist;
    dv(ie)        = res.dv;
    reached(ie)   = res.reached_ball;
    if res.reached_ball, t_hand_h(ie) = res.t_vec(res.conv_step)/3600; end
    fprintf('  e=%.2f: reached=%d, term=%.1f m, dv=%.4f, t=%.2f h\n', ...
        e_grid(ie), reached(ie), term_dist(ie), dv(ie), ...
        (res.conv_step * dt)/3600 * (res.reached_ball));
end
fprintf('Sweep complete in %.1fs.\n', toc(t0));

dv_infl = (dv - dv(1)) / dv(1) * 100;

%% Export JSON
metadata = struct('study','layer2_closing_robustness', ...
    'controller','CW-EMPC pure-fuel + CLF terminal cost', ...
    'plant','nonlinear Keplerian truth (perturbations off)', ...
    'scenario','10 km -> 300 m closing, N=20, dt=120 s, handover ball 50 m', ...
    'a_km', a/1000, 'date', datestr(now,'yyyy-mm-dd HH:MM:SS'));
data = struct('metadata', metadata, 'e_grid', e_grid, ...
    'term_dist_m', term_dist, 'dv_ms', dv, 'dv_inflation_pct', dv_infl, ...
    'time_to_handover_h', t_hand_h, 'reached_ball', double(reached));
fid = fopen(fullfile(out_json, 'closing_robustness.json'), 'w');
fprintf(fid, '%s', jsonencode(data, 'PrettyPrint', true));
fclose(fid);

%% Figure
figure('Position', [100 100 1000 420]);

subplot(1,2,1);
yyaxis left;
plot(e_grid, dv, '-o', 'LineWidth', 1.5); ylabel('total \Deltav [m/s]');
yyaxis right;
plot(e_grid, t_hand_h, '-s', 'LineWidth', 1.5); ylabel('time to handover [h]');
xlabel('eccentricity [-]'); grid on;
title('Closing: \Deltav & time to handover');

subplot(1,2,2);
plot(e_grid, dv_infl, '-o', 'LineWidth', 1.5); hold on;
bad = ~reached;
if any(bad)
    plot(e_grid(bad), dv_infl(bad), 'rx', 'MarkerSize', 10, 'LineWidth', 2, ...
        'DisplayName', 'handover NOT reached');
    legend('Location','northwest');
end
grid on; xlabel('eccentricity [-]'); ylabel('\Deltav inflation vs e=0 [%]');
title('Closing: fuel penalty vs eccentricity');

sgtitle('Layer 2 (closing) — CW-EMPC robustness, 10 km \rightarrow 300 m');
saveas(gcf, fullfile(out_fig, 'closing_robustness.png'));

fprintf('Exported JSON and figure.\n');