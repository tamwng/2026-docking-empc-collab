% sim_closing_robustness_drag.m
%
% Controller: CW-EMPC pure-fuel + CLF terminal (unchanged).
% Plant: truth with zonal J2-J4 + drag, NRLMSISE-00 DEFAULT space weather.
%
% Output: JSON + figure in validity_sweep/.

clear; clc; close all;

addpath('src/dynamics'); addpath('src/utils'); addpath('src/control');
constants;
warning('off', 'aero:atmosnrlmsise00:setf107af107aph');

out_json = 'exports/scenarios/validity_sweep';
out_fig  = 'results/figures/validity_sweep';
if ~exist(out_json, 'dir'), mkdir(out_json); end
if ~exist(out_fig,  'dir'), mkdir(out_fig);  end

%% Closing scenario (as sim_closing_robustness.m)
dt = 120;  N = 20;
x0   = [0; 10000; 0; 0; 0; 0];
x_sw = [0;   300; 0; 0; 0; 0];
u_max = 1e-2;

Ac = [0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1;
      3*n^2 0 0 0 2*n 0; 0 0 0 -2*n 0 0; 0 0 -n^2 0 0 0];
Bc = [zeros(3,3); eye(3)];
[Ad, Bd] = discretize(Ac, Bc, dt);
Q_dare = diag([1e-2 1e-2 1e-2 1 1 1]);  R_dare = diag([1e4 1e4 1e4]);
[~, P_clf] = lqr_controller(Ad, Bd, Q_dare, R_dare);

cost = struct('Q', zeros(6), 'R', eye(3), 'P', P_clf);
con  = struct('u_max', u_max, 'y_min_active', false);
con.los_cone.active = false;

params = struct('x0', x0, 'x_switch', x_sw, 'dt', dt, 'N', N, 'u_max', u_max, ...
                't_final', 5*T, 'conv_tol', 50, 'verbose', false);

%% Truth plant: circular (e=0), J2-J4 + drag ON
epoch = [2026 8 4 12 0 0];
Cd_t = 2.2;  A_t = 1.0;  m_t = 100;              % target Cd*A/m = 0.022
p_t = struct('mu', mu, 'epoch', epoch, 'use_zonal', true, 'zonal_degree', 4, ...
             'use_drag', true, 'Cd', Cd_t, 'A', A_t, 'm', m_t, 'omega_earth', omega_earth);
plant = struct('a', a, 'e', 0, 'i', deg2rad(51.6), 'p_t', p_t);

%% Sweep chaser ballistic-coefficient ratio (differential drag)
bc_ratio = [1 1.5 2 3 5];
nB = numel(bc_ratio);
dv        = zeros(1, nB);
term_dist = zeros(1, nB);
t_hand_h  = nan(1, nB);
reached   = false(1, nB);

fprintf('Closing drag sweep over %d BC ratios (J2+drag, e=0)...\n', nB);
t0 = tic;
for ib = 1:nB
    p_c = p_t;  p_c.Cd = Cd_t * bc_ratio(ib);    % scale chaser Cd*A/m
    plant.p_c = p_c;
    res = run_closing_truth(cost, con, params, plant);
    dv(ib)        = res.dv;
    term_dist(ib) = res.term_dist;
    reached(ib)   = res.reached_ball;
    if res.reached_ball, t_hand_h(ib) = res.t_vec(res.conv_step)/3600; end
    fprintf('  ratio=%.1f: reached=%d, term=%.1f m, dv=%.4f\n', ...
        bc_ratio(ib), reached(ib), term_dist(ib), dv(ib));
end
fprintf('Sweep complete in %.1fs.\n', toc(t0));

dv_infl = (dv - dv(1)) / dv(1) * 100;

%% Export JSON
metadata = struct('study','layer2_closing_robustness_drag', ...
    'controller','CW-EMPC pure-fuel + CLF terminal', ...
    'plant','truth J2-J4 + drag, e=0, NRLMSISE-00 default SW', ...
    'scenario','10 km -> 300 m closing, N=20, dt=120 s', ...
    'a_km', a/1000, 'date', datestr(now,'yyyy-mm-dd HH:MM:SS'));
data = struct('metadata', metadata, 'bc_ratio', bc_ratio, ...
    'dv_ms', dv, 'dv_inflation_pct', dv_infl, 'term_dist_m', term_dist, ...
    'time_to_handover_h', t_hand_h, 'reached_ball', double(reached));
fid = fopen(fullfile(out_json, 'closing_robustness_drag.json'), 'w');
fprintf(fid, '%s', jsonencode(data, 'PrettyPrint', true));
fclose(fid);

%% Figure
figure('Position', [100 100 640 420]);
yyaxis left;
plot(bc_ratio, dv, '-o', 'LineWidth', 1.5); ylabel('total \Deltav [m/s]');
yyaxis right;
plot(bc_ratio, dv_infl, '-s', 'LineWidth', 1.5); ylabel('\Deltav inflation vs matched BC [%]');
xlabel('drag ratio (Cd A/m)_c / (Cd A/m)_t'); grid on;
title('Closing (J2+drag, e=0): differential-drag fuel penalty');
saveas(gcf, fullfile(out_fig, 'closing_robustness_drag.png'));

fprintf('Exported JSON and figure.\n');