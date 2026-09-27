% sim_closing_regularizer_comparison.m
%
% Compares two ways of stabilising the economic (pure-fuel) closing MPC,
% swept over eccentricity against the truth plant:
%
%   'terminal'  Q=0, P=P_clf        stability from CLF terminal cost only
%               (the current thesis closing controller, no stage distortion)
%   'reg-*'     Q=eps*I6, P=0        stability from a strictly-dissipative
%               rotated stage cost (Gruene regularised economic MPC), no
%               terminal ingredient
%   'both'      Q=eps*I6, P=P_clf    regulariser + terminal cost
% Output: JSON + figure in validity_sweep/.

clear; clc; close all;

addpath('src/dynamics'); addpath('src/utils'); addpath('src/control');
addpath('tools/matlab2tikz-master/src');
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

Ac = [0 0 0 1 0 0;
      0 0 0 0 1 0;
      0 0 0 0 0 1;
      3*n^2 0 0 0 2*n 0;
      0 0 0 -2*n 0 0;
      0 0 -n^2 0 0 0];
Bc = [zeros(3,3); eye(3)];
[Ad, Bd] = discretize(Ac, Bc, dt);

% CLF terminal cost from the same DARE weights as sim_closing_robustness.m
Q_dare = diag([1e-2 1e-2 1e-2 1 1 1]);  R_dare = diag([1e4 1e4 1e4]);
[~, P_clf] = lqr_controller(Ad, Bd, Q_dare, R_dare);

con  = struct('u_max', u_max, 'y_min_active', false);
con.los_cone.active = false;

params = struct('x0', x0, 'x_switch', x_sw, 'dt', dt, 'N', N, 'u_max', u_max, ...
                't_final', 5*T, 'conv_tol', 50, 'verbose', false);

p_off = struct('mu', mu, 'use_zonal', false, 'use_drag', false);
plant = struct('a', a, 'i', deg2rad(51.6), 'p_t', p_off, 'p_c', p_off);

%% Controller configurations
configs = { ...
    struct('key','terminal', 'label','terminal cost ($Q=0$, $P_{clf}$)', ...
           'Q', zeros(6),      'P', P_clf), ...
    struct('key','reg1e-10',  'label','regulariser $\epsilon=10^{-10}$ (no term.)', ...
           'Q', 1e-10*eye(6),  'P', zeros(6)), ...
    struct('key','reg1e-8',   'label','regulariser $\epsilon=10^{-8}$ (no term.)', ...
           'Q', 1e-8*eye(6),   'P', zeros(6)), ...
    struct('key','both',      'label','$\epsilon=10^{-8}$ + terminal cost', ...
           'Q', 1e-8*eye(6),   'P', P_clf) };
nC = numel(configs);

e_grid = 0 : 0.02 : 0.30;
nE = numel(e_grid);

dv        = zeros(nC, nE);
term_dist = zeros(nC, nE);
reached   = false(nC, nE);

for ic = 1:nC
    cost = struct('Q', configs{ic}.Q, 'R', eye(3), 'P', configs{ic}.P);
    fprintf('=== config: %s ===\n', configs{ic}.key);
    for ie = 1:nE
        plant.e = e_grid(ie);
        res = run_closing_truth(cost, con, params, plant);
        dv(ic,ie)        = res.dv;
        term_dist(ic,ie) = res.term_dist;
        reached(ic,ie)   = res.reached_ball;
        fprintf('  e=%.2f: reached=%d, term=%.1f m, dv=%.4f\n', ...
            e_grid(ie), reached(ic,ie), term_dist(ic,ie), dv(ic,ie));
    end
end

% Fuel inflation vs each config's own e=0 baseline
dv_infl = (dv - dv(:,1)) ./ dv(:,1) * 100;

%% Nominal (e=0) readout: the regulariser's economic distortion
fprintf('\nNominal (e=0) fuel comparison:\n');
dv0_ref = dv(1,1);   % terminal-cost config is the economic reference
for ic = 1:nC
    fprintf('  %-10s : dv=%.4f  (%+.2f%% vs terminal-cost baseline)\n', ...
        configs{ic}.key, dv(ic,1), (dv(ic,1)-dv0_ref)/dv0_ref*100);
end

%% Export JSON
metadata = struct('study','closing_regularizer_comparison', ...
    'controller','CW economic closing MPC, terminal-cost vs regularised', ...
    'configs','terminal | reg1e-10 | reg1e-8 | both', ...
    'plant','nonlinear Keplerian truth (perturbations off)', ...
    'scenario','10 km -> 300 m closing, N=20, dt=120 s, handover ball 50 m', ...
    'a_km', a/1000, 'date', datestr(now,'yyyy-mm-dd HH:MM:SS'));
data = struct('metadata', metadata, 'e_grid', e_grid, ...
    'config_keys', {{configs{1}.key, configs{2}.key, configs{3}.key, configs{4}.key}}, ...
    'dv_ms', dv, 'dv_inflation_pct', dv_infl, ...
    'term_dist_m', term_dist, 'reached_ball', double(reached));
fid = fopen(fullfile(out_json, 'closing_regularizer_comparison.json'), 'w');
fprintf(fid, '%s', jsonencode(data, 'PrettyPrint', true));
fclose(fid);

%% Figure
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

figure('Position', [100 100 1000 420]);
mk = {'-o','-s','-^','-d'};

subplot(1,2,1);
for ic = 1:nC
    plot(e_grid, dv(ic,:), mk{ic}, 'LineWidth', 1.5, ...
        'DisplayName', configs{ic}.label); hold on;
end
bad = ~reached;
if any(bad(:))
    [ir, je] = find(bad);
    plot(e_grid(je), arrayfun(@(k) dv(ir(k),je(k)), 1:numel(je)), ...
        'rx', 'MarkerSize', 10, 'LineWidth', 2, 'HandleVisibility', 'off');
end
grid on; xlabel('eccentricity $e$ [-]'); ylabel('total $\Delta v$ [m/s]');
legend('Location', 'northwest');

subplot(1,2,2);
for ic = 1:nC
    plot(e_grid, dv_infl(ic,:), mk{ic}, 'LineWidth', 1.5, ...
        'DisplayName', configs{ic}.label); hold on;
end
grid on; xlabel('eccentricity $e$ [-]'); ylabel('$\Delta v$ inflation vs own $e{=}0$ [\%]');
legend('Location', 'northwest');

saveas(gcf, fullfile(out_fig, 'closing_regularizer_comparison.png'));
matlab2tikz(fullfile(out_fig, 'closing_regularizer_comparison.tikz'), ...
    'width', '\figurewidth', 'height', '\figureheight', 'showInfo', false);

fprintf('\nExported JSON and figure.\n');
