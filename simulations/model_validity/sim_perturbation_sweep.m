% sim_perturbation_sweep.m
% Question: with J2-J4 and drag ON, how far do the Keplerian relative models
% (CW, general, nonlinear) drift from the TRUTH within a fixed window, over
% the (altitude, differential-ballistic-coefficient) plane?
%
% Setup (chosen to isolate perturbations):
%   e = 0            -> CW not already killed by eccentricity
%   separation 100 m -> linearization assumption comfortably valid
%   models are fed the TRUE (perturbed) target trajectory, so the only error
%   source is the unmodeled DIFFERENTIAL perturbation on the chaser.
%
% Axes:
%   altitude  -> sets drag magnitude (strong low, negligible high); J2 ~ flat
%   drag ratio-> (Cd*A/m)_chaser / (Cd*A/m)_target : differential drag driver
%
% NRLMSISE-00 space weather: Aerospace Toolbox DEFAULT (F10.7/Ap unset).
%
% Combined 30-state ode113 per grid point:
%   [ X_t(6) ; X_c(6) ; Rho_nl(6) ; Rho_gen(6) ; Rho_cw(6) ]
% truth for X_t/X_c (own ballistic coeff each); Keplerian models on true X_t.

clear; clc; close all;

addpath('src/dynamics'); addpath('src/utils');
addpath('tools/matlab2tikz-master/src');
constants;

warning('off', 'aero:atmosnrlmsise00:setf107af107aph');

out_json = 'exports/scenarios/validity_sweep';
out_fig  = 'results/figures/validity_sweep';
if ~exist(out_json, 'dir'), mkdir(out_json); end
if ~exist(out_fig,  'dir'), mkdir(out_fig);  end

%% Sweep grid
alt_grid   = linspace(300, 800, 21) * 1e3;     % altitude [m]
bc_ratio   = linspace(1, 5, 21);               % (Cd*A/m)_c / (Cd*A/m)_t
n_orbits   = 5;                                % prediction window [orbits]

sep    = 100;                                  % fixed separation [m]
e_t    = 0;                                    % circular target
i_t    = deg2rad(51.6);
epoch  = [2026 8 4 12 0 0];                    % UTC epoch

% Target ballistic properties (chaser scales Cd*A/m by bc_ratio)
Cd_t = 2.2;  A_t = 1.0;  m_t = 100;            % -> Cd*A/m = 0.022 m^2/kg

opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);

nA = numel(alt_grid);  nB = numel(bc_ratio);
err_nl  = zeros(nA, nB);
err_gen = zeros(nA, nB);
err_cw  = zeros(nA, nB);

% Capture time series at the worst corner (lowest alt, highest ratio)
corner_series = struct();

fprintf('Perturbation sweep %d x %d = %d points...\n', nA, nB, nA*nB);
t0 = tic;

for ia = 1:nA
    a_t   = R_earth + alt_grid(ia);
    n_mm  = sqrt(mu / a_t^3);
    T_orb = 2*pi*sqrt(a_t^3/mu);
    t_span = linspace(0, n_orbits*T_orb, 400);

    X_t0 = oe2eci(a_t, e_t, i_t, 0, 0, 0, mu);

    % CW-bounded relative IC scaled to 100 m (no CW secular drift)
    s    = sep / sqrt(2);
    Rho0 = [s; 0; s; 0; -2*n_mm*s; 0];
    X_c0 = hill2eci(X_t0, Rho0);

    % Target perturbation params (fixed BC)
    p_t = struct('mu', mu, 'epoch', epoch, 'use_zonal', true, ...
                 'zonal_degree', 4, 'use_drag', true, ...
                 'Cd', Cd_t, 'A', A_t, 'm', m_t, 'omega_earth', omega_earth);

    for ib = 1:nB
        % Chaser BC: scale Cd*A/m by bc_ratio (higher = more drag)
        p_c = p_t;
        p_c.Cd = Cd_t * bc_ratio(ib);          % simplest way to scale Cd*A/m

        odefun = @(t, y) [ ...
            truth_propagator_at(t, y(1:6),   p_t); ...            % target truth
            truth_propagator_at(t, y(7:12),  p_c); ...            % chaser truth
            relative_motion_nonlinear(y(13:18), y(1:6), mu); ...  % nonlinear
            relative_motion_general(  y(19:24), y(1:6), mu); ...  % general
            clohessy_wiltshire(       y(25:30), n_mm) ];          % CW

        y0 = [X_t0; X_c0; Rho0; Rho0; Rho0];
        [tt, Y] = ode113(odefun, t_span, y0, opts);

        % Reference: differenced truth
        Nt = numel(tt);
        ref = zeros(Nt, 6);
        for k = 1:Nt
            ref(k, :) = eci2hill(Y(k, 1:6).', Y(k, 7:12).').';
        end

        pe_nl  = vecnorm(Y(:, 13:15) - ref(:, 1:3), 2, 2);
        pe_gen = vecnorm(Y(:, 19:21) - ref(:, 1:3), 2, 2);
        pe_cw  = vecnorm(Y(:, 25:27) - ref(:, 1:3), 2, 2);

        err_nl(ia, ib)  = max(pe_nl)  / sep;
        err_gen(ia, ib) = max(pe_gen) / sep;
        err_cw(ia, ib)  = max(pe_cw)  / sep;

        if ia == 1 && ib == nB
            corner_series.t_orbits = tt / T_orb;
            corner_series.nl  = pe_nl;
            corner_series.gen = pe_gen;
            corner_series.cw  = pe_cw;
        end
    end
    fprintf('  alt = %3.0f km done (%.1fs)\n', alt_grid(ia)/1e3, toc(t0));
end
fprintf('Sweep complete in %.1fs.\n', toc(t0));

%% Export JSON 
metadata = struct( ...
    'study',        'layer1_perturbation_model_error', ...
    'reference',    'differenced truth (J2-J4 + drag, NRLMSISE-00 default SW)', ...
    'models_fed',   'true perturbed target trajectory', ...
    'e',            e_t, ...
    'separation_m', sep, ...
    'i_deg',        rad2deg(i_t), ...
    'n_orbits',     n_orbits, ...
    'date',         datestr(now, 'yyyy-mm-dd HH:MM:SS'));

data = struct( ...
    'metadata',      metadata, ...
    'alt_grid_km',   alt_grid/1e3, ...
    'bc_ratio',      bc_ratio, ...
    'err_nonlinear', err_nl, ...
    'err_general',   err_gen, ...
    'err_cw',        err_cw, ...
    'corner_series', corner_series);   % worst-corner time history, for the right-panel figure

fid = fopen(fullfile(out_json, 'model_error_alt_bc.json'), 'w');
fprintf(fid, '%s', jsonencode(data, 'PrettyPrint', true));
fclose(fid);

%% Figures
% Heatmap via imagesc, exported by matlab2tikz as vector \fill rectangles
% ('imagesAsPng' false -- no PNG) instead of contourf -- see
% sim_model_validity_sweep.m for why.
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

bc_label = '$\frac{\beta_c}{\beta_t}$, $\beta = C_dA/m$ [-]';

% Two separate figures/tikz exports, not subplots of one figure: a
% colorbar-bearing image axis sharing a tikzpicture with an unrelated plain
% axis is a confirmed matlab2tikz bug in this vendored version -- the plain
% axis's (log, 1e-10 to 1e5) y-range leaks into the colorbar's tick labels,
% verified by test-compiling the combined export through pdflatex.
figure('Position', [100 100 520 420]);
imagesc(bc_ratio, alt_grid/1e3, log10(err_nl)); set(gca, 'YDir', 'normal');
cb = colorbar; cb.Label.Interpreter = 'latex'; cb.Label.String = '$\log_{10}$(max rel.\ pos.\ error) [-]';
xlabel(bc_label); ylabel('altitude [km]');
saveas(gcf, fullfile(out_fig, 'model_error_alt_bc_heatmap.png'));
matlab2tikz(fullfile(out_fig, 'model_error_alt_bc_heatmap.tikz'), ...
    'width', '\figurewidth', 'height', '\figureheight', 'showInfo', false, ...
    'imagesAsPng', false);

figure('Position', [650 100 520 420]);
semilogy(corner_series.t_orbits, corner_series.nl,  'LineWidth', 1.5); hold on;
semilogy(corner_series.t_orbits, corner_series.gen, '--', 'LineWidth', 1.5);
semilogy(corner_series.t_orbits, corner_series.cw,  ':', 'LineWidth', 1.5);
grid on; xlabel('orbits'); ylabel('position error [m] ($\log_{10}$ scale)');
% The error starts at ~1e-10 (models coincide at t=0) then climbs past 1
% within the first ~5% of orbit 1 -- letting the log axis autoscale down
% to 1e-10 wastes ~9 empty decades. Floor it just below where the climb
% actually starts so the full 5-orbit growth trend fills the plot.
% ylim must be set BEFORE legend(): setting it after silently breaks
% matlab2tikz's exported legend position (verified -- it fell back to a
% southwest placement regardless of the requested Location).
ylim([1e-1, 1e5]);
legend('nonlinear', 'general', 'CW', 'Location', 'southeast');
saveas(gcf, fullfile(out_fig, 'model_error_alt_bc_corner.png'));
matlab2tikz(fullfile(out_fig, 'model_error_alt_bc_corner.tikz'), ...
    'width', '\figurewidth', 'height', '\figureheight', 'showInfo', false);

fprintf('Exported JSON and figure.\n');
