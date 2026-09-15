% sim_nmc_terminal_comparison.m
%
% Compares two NMC terminal INGREDIENTS for the economic (pure-fuel)
% fly-around controller, swept over eccentricity against the truth plant:
%
%   'phase'    periodic terminal equality x(N)=Pi*(k+N mod P)  (S2 Run 2)
%   'fixed'    single fixed point on Pi*, x(N)=x_star(:,1)
% Output: JSON + figure in validity_sweep/.

clear; clc; close all;

addpath('src/dynamics'); addpath('src/utils'); addpath('src/control');
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
    'H', 2*eye(3*N), 'A_u', [eye(3*N); -eye(3*N)], ...
    'b_u', u_max*ones(6*N,1));

x0 = [0; 300; 0; 0; 0; 0];   % V-bar hold IC (as S2)
params = struct('x0', x0, 'rho_min', 50, 'rho_max', 200, 'u_max', u_max, ...
                'N', N, 'P', P, 'n', n, 'dt', dt, 'n_sim', n_sim, 'verbose', false);

p_off = struct('mu', mu, 'use_zonal', false, 'use_drag', false);
plant = struct('a', a, 'i', deg2rad(51.6), 'p_t', p_off, 'p_c', p_off);

%% Two terminal ingredients (both pure-fuel, no band)
modes = { ...
    struct('key','phase',    'label','phase-synchronised x(N)=\Pi^*(k)'), ...
    struct('key','fixed',    'label','fixed point x(N)=\Pi^*(1)') };
nM = numel(modes);

e_grid = 0 : 0.02 : 0.30;
nE = numel(e_grid);

dv_last_orbit = zeros(nM, nE);   % maintenance rate over final orbit
range_last    = zeros(nM, nE);   % orbit integrity, final orbit
n_infeas      = zeros(nM, nE);   % QP infeasible steps

for im = 1:nM
    run_cfg = struct('use_nmc_manifold', true, 'use_periodic_terminal', true, ...
                     'use_band', false, 'x_star', x_star, ...
                     'terminal_mode', modes{im}.key, 'fixed_idx', 1);
    fprintf('=== terminal mode: %s ===\n', modes{im}.key);
    for ie = 1:nE
        plant.e = e_grid(ie);
        res = run_flyaround_truth(run_cfg, matrices, params, plant);
        dv_last_orbit(im,ie) = sum(res.dv_log(max(1, n_sim-P+1):end));
        range_last(im,ie)    = max(res.range_log(max(1, n_sim-P+1):end));
        n_infeas(im,ie)      = res.n_infeasible;
        fprintf('  e=%.2f: dv/orbit=%.4f, range_last=%.0f m, infeas=%d\n', ...
            e_grid(ie), dv_last_orbit(im,ie), range_last(im,ie), n_infeas(im,ie));
    end
end

%% Nominal (e=0) readout: why fixed is wrong even without eccentricity
fprintf('\nNominal (e=0) comparison:\n');
for im = 1:nM
    fprintf('  %-10s : dv/orbit=%.4f, range_last=%.0f m, infeas=%d\n', ...
        modes{im}.key, dv_last_orbit(im,1), range_last(im,1), n_infeas(im,1));
end

%% Export JSON
metadata = struct('study','nmc_terminal_comparison', ...
    'controller','economic fly-around EMPC (pure fuel), 2 terminal ingredients', ...
    'modes','phase | fixed', ...
    'plant','nonlinear Keplerian truth (perturbations off)', ...
    'scenario', sprintf('NMC b=%.0f m, V-bar 300 m IC, dt=T/%d, N=P=%d, %d orbits', ...
                        b_nmc, P, P, n_sim/P), ...
    'a_km', a/1000, 'date', datestr(now,'yyyy-mm-dd HH:MM:SS'));
data = struct('metadata', metadata, 'e_grid', e_grid, ...
    'mode_keys', {cellfun(@(m) m.key, modes, 'UniformOutput', false)}, ...
    'dv_last_orbit_ms', dv_last_orbit, 'range_last_orbit_m', range_last, ...
    'n_infeasible', n_infeas);
fid = fopen(fullfile(out_json, 'nmc_terminal_comparison.json'), 'w');
fprintf(fid, '%s', jsonencode(data, 'PrettyPrint', true));
fclose(fid);

%% Figure
figure('Position', [100 100 1000 420]);
mk = {'-o','-s','-^'};

subplot(1,2,1);
for im = 1:nM
    semilogy(e_grid, max(dv_last_orbit(im,:), 1e-4), mk{im}, 'LineWidth', 1.5, ...
        'DisplayName', modes{im}.label); hold on;
end
grid on; xlabel('eccentricity [-]'); ylabel('maintenance \Deltav / orbit [m/s]');
legend('Location','northwest'); title('Maintenance cost by terminal ingredient');

subplot(1,2,2);
for im = 1:nM
    semilogy(e_grid, max(range_last(im,:), 1), mk{im}, 'LineWidth', 1.5, ...
        'DisplayName', modes{im}.label); hold on;
end
yline(150, '--', '2:1 envelope (\approx150 m)', 'HandleVisibility','off');
grid on; xlabel('eccentricity [-]'); ylabel('max range, final orbit [m]');
legend('Location','northwest'); title(sprintf('Orbit integrity (b = %.0f m)', b_nmc));

sgtitle('NMC terminal ingredients - phase-sync vs fixed point');
saveas(gcf, fullfile(out_fig, 'nmc_terminal_comparison.png'));

fprintf('\nExported JSON and figure.\n');
