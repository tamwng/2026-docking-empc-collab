% sim_nmc_robustness_empc_drag.m
%
% Controller: CW condensed QP, pure fuel, periodic terminal equality (Run 2).
% Plant: truth J2-J4 + drag, e=0, NRLMSISE-00 DEFAULT space weather.
%
% Output: JSON + figure in validity_sweep/.

clear; clc; close all;

addpath('src/dynamics'); addpath('src/utils'); addpath('src/control');
addpath('tools/matlab2tikz-master/src');
constants;
warning('off', 'aero:atmosnrlmsise00:setf107af107aph');

out_json = 'exports/scenarios/validity_sweep';
out_fig  = 'results/figures/validity_sweep';
if ~exist(out_json, 'dir'), mkdir(out_json); end
if ~exist(out_fig,  'dir'), mkdir(out_fig);  end

%% S2 discretisation + CW-NMC orbit Pi* 
P  = 92;
dt = T / P;
Ac = [0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1;
      3*n^2 0 0 0 2*n 0; 0 0 0 -2*n 0 0; 0 0 -n^2 0 0 0];
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
    'H', 2*eye(3*N), 'A_u', [eye(3*N); -eye(3*N)], 'b_u', u_max*ones(6*N,1));

x0 = [0; 300; 0; 0; 0; 0];
params = struct('x0', x0, 'rho_min', 50, 'rho_max', 200, 'u_max', u_max, ...
                'N', N, 'P', P, 'n', n, 'dt', dt, 'n_sim', n_sim, 'verbose', false);
run_cfg = struct('use_nmc_manifold', true, 'use_periodic_terminal', true, ...
                 'use_band', false, 'x_star', x_star);

%% Truth plant: circular (e=0), J2-J4 + drag ON
epoch = [2026 8 4 12 0 0];
Cd_t = 2.2;  A_t = 1.0;  m_t = 100;
p_t = struct('mu', mu, 'epoch', epoch, 'use_zonal', true, 'zonal_degree', 4, ...
             'use_drag', true, 'Cd', Cd_t, 'A', A_t, 'm', m_t, 'omega_earth', omega_earth);
plant = struct('a', a, 'e', 0, 'i', deg2rad(51.6), 'p_t', p_t);

%% Sweep chaser ballistic-coefficient ratio
bc_ratio = [1 1.5 2 3 5];
nB = numel(bc_ratio);
dv_last_orbit = zeros(1, nB);
dv_total      = zeros(1, nB);
range_last    = zeros(1, nB);
n_infeas      = zeros(1, nB);

fprintf('EMPC NMC drag sweep over %d BC ratios (J2+drag, e=0)...\n', nB);
t0 = tic;
for ib = 1:nB
    p_c = p_t;  p_c.Cd = Cd_t * bc_ratio(ib);
    plant.p_c = p_c;
    res = run_flyaround_truth(run_cfg, matrices, params, plant);
    dv_last_orbit(ib) = sum(res.dv_log(max(1, n_sim-P+1):end));
    dv_total(ib)      = res.dv_total;
    range_last(ib)    = max(res.range_log(max(1, n_sim-P+1):end));
    n_infeas(ib)      = res.n_infeasible;
    fprintf('  ratio=%.1f: dv/orbit=%.4f, dv_tot=%.4f, range_last=%.0f m, infeas=%d\n', ...
        bc_ratio(ib), dv_last_orbit(ib), dv_total(ib), range_last(ib), n_infeas(ib));
end
fprintf('Sweep complete in %.1fs.\n', toc(t0));

%% Export JSON
metadata = struct('study','layer2_nmc_robustness_empc_drag', ...
    'controller','economic fly-around EMPC (S2 Run 2: pure fuel + periodic terminal eq x(N)=Pi*)', ...
    'plant','truth J2-J4 + drag, e=0, NRLMSISE-00 default SW', ...
    'scenario', sprintf('NMC b=%.0f m, V-bar 300 m IC, dt=T/%d, N=P=%d, %d orbits', b_nmc, P, P, n_sim/P), ...
    'a_km', a/1000, 'date', datestr(now,'yyyy-mm-dd HH:MM:SS'));
data = struct('metadata', metadata, 'bc_ratio', bc_ratio, ...
    'dv_last_orbit_ms', dv_last_orbit, 'dv_total_ms', dv_total, ...
    'range_last_orbit_m', range_last, 'n_infeasible', n_infeas);
fid = fopen(fullfile(out_json, 'nmc_robustness_empc_drag.json'), 'w');
fprintf(fid, '%s', jsonencode(data, 'PrettyPrint', true));
fclose(fid);

%% Figure
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

figure('Position', [100 100 1000 420]);
bc_label = '$\frac{\beta_c}{\beta_t}$, $\beta = C_dA/m$ [-]';

subplot(1,2,1);
plot(bc_ratio, dv_last_orbit, '-o', 'LineWidth', 1.5); hold on;
plot(bc_ratio, dv_total, '--s', 'LineWidth', 1.5);
xlabel(bc_label); ylabel('$\Delta v$ [m/s]'); grid on;
legend('maintenance $\Delta v$ (last orbit)', 'total $\Delta v$', 'Location', 'northwest');

subplot(1,2,2);
plot(bc_ratio, range_last, '-^', 'LineWidth', 1.5);
xlabel(bc_label); ylabel('max range, final orbit [m]'); grid on;

saveas(gcf, fullfile(out_fig, 'nmc_robustness_empc_drag.png'));
matlab2tikz(fullfile(out_fig, 'nmc_robustness_empc_drag.tikz'), ...
    'width', '\figurewidth', 'height', '\figureheight', 'showInfo', false);

fprintf('Exported JSON and figure.\n');