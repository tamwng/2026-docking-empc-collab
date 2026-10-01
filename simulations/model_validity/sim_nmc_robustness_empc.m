% sim_nmc_robustness_empc.m
%
% Controller: CW condensed QP, pure fuel, periodic terminal equality (Run 2).
% Plant: nonlinear Keplerian truth (perturbations off) at eccentricity e.
%
% Output: JSON + figure in validity_sweep/.

clear; clc; close all;

addpath('src/dynamics'); addpath('src/utils'); addpath('src/control');
addpath('tools/matlab2tikz-master/src');
constants;

out_json = 'exports/scenarios/validity_sweep';
out_fig  = 'results/figures/validity_sweep';
if ~exist(out_json, 'dir'), mkdir(out_json); end
if ~exist(out_fig,  'dir'), mkdir(out_fig);  end

%% S2 discretisation: dt = T/P so the CW-NMC closes exactly under A_d
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

%% Prescribed CW-NMC orbit Pi* (2:1 ellipse, b_nmc radial semi-axis)
b_nmc = 75;
x_star = zeros(6, P);
x_star(:,1) = [b_nmc; 0; 0; 0; -2*n*b_nmc; 0];
for k = 1:P-1, x_star(:,k+1) = A_d * x_star(:,k); end
assert(norm(A_d*x_star(:,end) - x_star(:,1)) < 1e-6, 'Pi* does not close under A_d');

%% Condensed QP matrices (as sim_scenario2_flyaround.m)
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

x0 = [0; 300; 0; 0; 0; 0];   % V-bar hold IC (as S2)
params = struct('x0', x0, 'rho_min', 50, 'rho_max', 200, 'u_max', u_max, ...
                'N', N, 'P', P, 'n', n, 'dt', dt, 'n_sim', n_sim, 'verbose', false);

% Run 2 config: periodic terminal equality, no band (Pi* pins the orbit)
run_cfg = struct('use_nmc_manifold', true, 'use_periodic_terminal', true, ...
                 'use_band', false, 'x_star', x_star);

p_off = struct('mu', mu, 'use_zonal', false, 'use_drag', false);
plant = struct('a', a, 'i', deg2rad(51.6), 'p_t', p_off, 'p_c', p_off);

%% Sweep eccentricity
e_grid = 0 : 0.02 : 0.30;
nE = numel(e_grid);
dv_total  = zeros(1, nE);
n_infeas  = zeros(1, nE);
dv_last_orbit = zeros(1, nE);   % maintenance rate = dv over final orbit
max_range     = zeros(1, nE);   % orbit integrity: bounded (~150 m) vs lost (km)
range_last    = zeros(1, nE);   % max range in the final orbit

fprintf('EMPC NMC robustness sweep (%d eccentricities, %d steps each)...\n', nE, n_sim);
t0 = tic;
for ie = 1:nE
    plant.e = e_grid(ie);
    res = run_flyaround_truth(run_cfg, matrices, params, plant);
    dv_total(ie)  = res.dv_total;
    n_infeas(ie)  = res.n_infeasible;
    dv_last_orbit(ie) = sum(res.dv_log(max(1, n_sim-P+1):end));
    max_range(ie)     = max(res.range_log);
    range_last(ie)    = max(res.range_log(max(1, n_sim-P+1):end));
    fprintf('  e=%.2f: dv/orbit=%.4f, dv_tot=%.4f, range_last=%.0f m, maxrange=%.0f m, infeas=%d\n', ...
        e_grid(ie), dv_last_orbit(ie), dv_total(ie), ...
        range_last(ie), max_range(ie), n_infeas(ie));
end
fprintf('Sweep complete in %.1fs.\n', toc(t0));

% Operational limit: last eccentricity before the controller breaks
% (QP stays feasible AND the orbit stays bounded, say < 5x the 150 m envelope)
ok = (n_infeas == 0) & (range_last < 750);
e_limit = e_grid(find(~ok, 1) - 1);
if isempty(e_limit), e_limit = e_grid(end); end
fprintf('Operational eccentricity limit (feasible + bounded orbit): e <= %.2f\n', e_limit);

%% Export JSON
metadata = struct('study','layer2_nmc_robustness_empc', ...
    'controller','economic fly-around EMPC (S2 Run 2: pure fuel + periodic terminal eq x(N)=Pi*)', ...
    'plant','nonlinear Keplerian truth (perturbations off)', ...
    'scenario', sprintf('NMC b=%.0f m, V-bar 300 m IC, dt=T/%d, N=P=%d, %d orbits', b_nmc, P, P, n_sim/P), ...
    'a_km', a/1000, 'date', datestr(now,'yyyy-mm-dd HH:MM:SS'));
data = struct('metadata', metadata, 'e_grid', e_grid, ...
    'dv_total_ms', dv_total, 'dv_last_orbit_ms', dv_last_orbit, ...
    'max_range_m', max_range, 'range_last_orbit_m', range_last, ...
    'n_infeasible', n_infeas, 'e_operational_limit', e_limit);
fid = fopen(fullfile(out_json, 'nmc_robustness_empc.json'), 'w');
fprintf(fid, '%s', jsonencode(data, 'PrettyPrint', true));
fclose(fid);

%% Figure
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

figure('Position', [100 100 1000 380]);

% Panel (a): maintenance Delta-v per orbit (log). No e_limit threshold line
% -- the feasibility cliff is shown directly in panel (b) instead.
ax1 = subplot(1,2,1);
semilogy(e_grid, max(dv_last_orbit, 1e-4), '-o', 'LineWidth', 1.5);
xlabel('eccentricity $e$ [-]'); ylabel('maint.\ $\Delta v$ / orbit [m/s]'); grid on;
set(gca, 'YMinorGrid', 'off');   % dense minor gridlines over 4+ decades clutter a short panel

% Panel (b): QP infeasibility onset
ax2 = subplot(1,2,2);
plot(e_grid, n_infeas, '-s', 'LineWidth', 1.5);
xlabel('eccentricity $e$ [-]'); ylabel('\# infeasible QP steps'); grid on;

% Equal-width panels, explicit generous gap -- same numbers as the other
% two-panel figures in this section.
pos1 = get(ax1, 'Position'); pos2 = get(ax2, 'Position');
left0 = 0.02; pw = 0.346; gap = 0.253;
pos1(1) = left0;            pos1(3) = pw;
pos2(1) = left0 + pw + gap; pos2(3) = pw; pos2(2) = pos1(2); pos2(4) = pos1(4);
set(ax1, 'Position', pos1); set(ax2, 'Position', pos2);

saveas(gcf, fullfile(out_fig, 'nmc_robustness_empc.png'));
matlab2tikz(fullfile(out_fig, 'nmc_robustness_empc.tikz'), ...
    'width', '\figurewidth', 'height', '\figureheight', 'showInfo', false);

fprintf('Exported JSON and figure.\n');