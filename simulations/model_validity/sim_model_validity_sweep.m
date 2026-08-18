% sim_model_validity_sweep.m
% Question: over the (eccentricity, separation) plane, how far do the
% linear relative-motion models (CW, general-elliptic) drift from the exact
% nonlinear Keplerian relative motion within a fixed prediction window?
%
% Perturbations are OFF, so relative_motion_nonlinear is EXACT (validated by
% tests/validate_truth_vs_nonlinear.m) and serves as the reference. This
% isolates the two Keplerian modelling assumptions:
%   - CW      : circular target orbit  -> breaks with eccentricity
%   - general : |rho| << r_t linearization -> breaks with separation
%
% For each grid point a single ode113 call propagates a combined 24-state
%   [ X_t(6) ; Rho_nl(6) ; Rho_gen(6) ; Rho_cw(6) ]
% so every model sees the identical target trajectory and time grid.
%
% Output: JSON grid to exports/scenarios/validity_sweep/, contour figures to
% results/figures/validity_sweep/.

clear; clc; close all;

addpath('src/dynamics'); addpath('src/utils');
constants;

out_json = 'exports/scenarios/validity_sweep';
out_fig  = 'results/figures/validity_sweep';
if ~exist(out_json, 'dir'), mkdir(out_json); end
if ~exist(out_fig,  'dir'), mkdir(out_fig);  end

%% Sweep grid
e_grid   = 0 : 0.02 : 0.30;                       % eccentricity
sep_grid = logspace(1, 6, 16);                    % separation 10 m -> 1000 km
n_orbits = 3;                                     % prediction window [orbits]

% Fixed target geometry (only e is swept here; altitude fixed at baseline a)
i_t    = deg2rad(51.6);
RAAN_t = 0;  argp_t = 0;  nu_t = 0;               % start at perigee
n_mm   = sqrt(mu / a^3);                          % mean motion for CW

opts   = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);

nE = numel(e_grid);  nS = numel(sep_grid);

% Result matrices: max normalized position error over the window [-]
err_cw   = zeros(nE, nS);
err_gen  = zeros(nE, nS);
% Secular signature: max normalized ALONG-TRACK (y) error [-]
drift_cw  = zeros(nE, nS);
drift_gen = zeros(nE, nS);

fprintf('Sweeping %d x %d = %d grid points...\n', nE, nS, nE*nS);
t0 = tic;

for ie = 1:nE
    e_t = e_grid(ie);

    % Target ECI IC and its period (depends on a only -> constant here)
    X_t0 = oe2eci(a, e_t, i_t, RAAN_t, argp_t, nu_t, mu);
    T_orb = 2*pi*sqrt(a^3/mu);
    t_span = [0, n_orbits*T_orb];

    for is = 1:nS
        sep = sep_grid(is);

        % CW-bounded relative IC scaled to this separation (no secular drift
        % under CW): x0 = z0 = sep/sqrt(2), vy0 = -2 n x0, rest zero.
        s   = sep / sqrt(2);
        Rho0 = [s; 0; s; 0; -2*n_mm*s; 0];

        p_off = struct('mu', mu, 'use_zonal', false, 'use_drag', false);

        odefun = @(t, y) [ ...
            truth_propagator_at(t, y(1:6), p_off); ...              % target
            relative_motion_nonlinear(y(7:12),  y(1:6), mu); ...    % ref
            relative_motion_general(  y(13:18), y(1:6), mu); ...    % general
            clohessy_wiltshire(       y(19:24), n_mm) ];            % CW

        y0 = [X_t0; Rho0; Rho0; Rho0];
        [~, Y] = ode113(odefun, t_span, y0, opts);

        ref = Y(:, 7:12);
        gen = Y(:, 13:18);
        cw  = Y(:, 19:24);

        pe_cw  = vecnorm(cw(:,1:3)  - ref(:,1:3), 2, 2);
        pe_gen = vecnorm(gen(:,1:3) - ref(:,1:3), 2, 2);

        err_cw(ie, is)  = max(pe_cw)  / sep;
        err_gen(ie, is) = max(pe_gen) / sep;
        drift_cw(ie, is)  = max(abs(cw(:,2)  - ref(:,2)))  / sep;
        drift_gen(ie, is) = max(abs(gen(:,2) - ref(:,2))) / sep;
    end
    fprintf('  e = %.2f done (%.1fs)\n', e_t, toc(t0));
end
fprintf('Sweep complete in %.1fs.\n', toc(t0));

%% Export JSON 
metadata = struct( ...
    'study',        'layer1_open_loop_model_error', ...
    'reference',    'relative_motion_nonlinear (perturbations off)', ...
    'a_km',         a/1000, ...
    'i_deg',        rad2deg(i_t), ...
    'n_orbits',     n_orbits, ...
    'date',         datestr(now, 'yyyy-mm-dd HH:MM:SS'));

data = struct( ...
    'metadata',   metadata, ...
    'e_grid',     e_grid, ...
    'sep_grid_m', sep_grid, ...
    'err_cw',     err_cw, ...        % max ||pos err|| / sep
    'err_gen',    err_gen, ...
    'drift_cw',   drift_cw, ...      % max |along-track err| / sep
    'drift_gen',  drift_gen);

fid = fopen(fullfile(out_json, 'model_error_e_sep.json'), 'w');
fprintf(fid, '%s', jsonencode(data, 'PrettyPrint', true));
fclose(fid);

%% Contour figures 
[SEP, E] = meshgrid(sep_grid, e_grid);

figure('Position', [100 100 1000 420]);

subplot(1,2,1);
contourf(SEP, E, log10(err_cw), 20, 'LineColor', 'none'); hold on;
[C,h] = contour(SEP, E, err_cw, [0.01 0.01], 'w-', 'LineWidth', 2);
clabel(C, h, 'Color', 'w');
set(gca, 'XScale', 'log'); colorbar;
xlabel('separation [m]'); ylabel('eccentricity [-]');
title('CW model: log_{10}(max rel. pos. error)');

subplot(1,2,2);
contourf(SEP, E, log10(err_gen), 20, 'LineColor', 'none'); hold on;
[C,h] = contour(SEP, E, err_gen, [0.01 0.01], 'w-', 'LineWidth', 2);
clabel(C, h, 'Color', 'w');
set(gca, 'XScale', 'log'); colorbar;
xlabel('separation [m]'); ylabel('eccentricity [-]');
title('General model: log_{10}(max rel. pos. error)');

sgtitle(sprintf('Open-loop model error vs nonlinear (%.0f orbits)', n_orbits));
saveas(gcf, fullfile(out_fig, 'model_error_e_sep.png'));

fprintf('Exported JSON and figure.\n');
