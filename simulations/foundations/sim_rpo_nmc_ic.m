% sim_rpo_nmc_ic.m
% Relative motion (CWH) demo of the NMC drift-free condition, for
% orbital_viz.m. NO CONTROLLER: pure open-loop orbital dynamics, exactly
% like sim_rpo.m, but exported twice from the same chaser starting point:
%
%   IC MET      ydot0 = -2*n*x0  -> bounded 2:1 NMC ellipse (Pi*, b=75 m)
%   IC VIOLATED ydot0 = 0        -> secular along-track drift on top of it
%
% Both cases share the same target orbit and chaser starting POSITION
% (x0, y0, z0); only the initial along-track velocity differs, isolating
% the drift-free condition (Schaub & Junkins, Ch. 14) as the sole cause of
% the different long-term behaviour. Reference case for
% simulations/empc/sim_nmc_free_drift.m (same condition, analytical/TikZ
% only); this script instead exports the orbital_viz.m JSON schema
% { metadata, t, X_t, Rho } over several orbits so the divergence is
% visible in playback.
%
% Output: exports/scenarios/rpo_nmc_ic_met.json
%         exports/scenarios/rpo_nmc_ic_violated.json

clear; clc;
script_dir = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(script_dir, '..', '..', 'src')));
constants;   % loads a, e, i, mu, n, T, R_earth, etc.

%% Shared setup
dynamics_model = 'CWH';
n_orbits       = 3;
dt             = T / 300;                 % 300 steps/orbit -> smooth playback
n_steps        = n_orbits * round(T / dt);
t              = (0:n_steps-1) * dt;

% Target: circular orbit, equatorial (same convention as sim_rpo.m)
rx0_t = a;  ry0_t = 0;  rz0_t = 0;
vc    = sqrt(mu / a);
X_t0  = [rx0_t; ry0_t; rz0_t; 0; vc; 0];

f_target = @(X) two_body(X, mu);
X_t = rk4_integrator(f_target, X_t0, dt, n_steps);

% Chaser: bottom of a 2:1 NMC ellipse, radial semi-axis b_nmc (S2 scale)
% Cross-track (z) motion is a pure SHM at the orbital rate n, fully
% decoupled from the in-plane radial/along-track dynamics (Schaub 14.13);
% adding it tilts the chaser's relative orbit out of the target's plane
% for visual differentiation in the ECI panel, without altering the
% drift-free condition or the met/violated comparison at all.
b_nmc  = 75;             % [m]
x0     = -b_nmc;
y0     = 0;
z0     = 0.5 * b_nmc;    % [m]  out-of-plane amplitude
x_dot0 = 0;
z_dot0 = n * z0;         % quarter-period phase shift -> circular-looking tilt

cases = { ...
    struct('key', 'ic_met',      'y_dot0', -2*n*x0, ...
           'label', 'IC met: ydot0 = -2 n x0 (bounded NMC ellipse)'), ...
    struct('key', 'ic_violated', 'y_dot0', 0, ...
           'label', 'IC violated: ydot0 = 0 (secular along-track drift)') };

f_rel = @(Rho, ~) clohessy_wiltshire(Rho, n);

out_dir = fullfile(script_dir, '..', '..', 'exports', 'scenarios');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

for ci = 1:numel(cases)
    cs = cases{ci};
    Rho0 = [x0; y0; z0; x_dot0; cs.y_dot0; z_dot0];

    Rho = zeros(6, n_steps);
    Rho(:, 1) = Rho0;
    for k = 1:n_steps-1
        xk = Rho(:, k);
        k1 = f_rel(xk);
        k2 = f_rel(xk + 0.5*dt*k1);
        k3 = f_rel(xk + 0.5*dt*k2);
        k4 = f_rel(xk + dt*k3);
        Rho(:, k+1) = xk + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
    end

    fprintf('%s\n', cs.label);
    fprintf('  final range: %.1f m   (bounded case stays ~%.0f m, drift case grows)\n', ...
        norm(Rho(1:3, end)), b_nmc*2);

    metadata = struct( ...
        'a_km',           a / 1000, ...
        'e',              e, ...
        'i_deg',          rad2deg(i), ...
        'period_min',     T / 60, ...
        'altitude_km',    (a - R_earth) / 1000, ...
        'dt',             dt, ...
        'n_steps',        n_steps, ...
        'n_orbits',       n_orbits, ...
        'dynamics_model', dynamics_model, ...
        'scenario',       cs.key, ...
        'b_nmc_m',        b_nmc, ...
        'x0_m',           x0, 'y0_m', y0, 'z0_m', z0, ...
        'x_dot0_mps',     x_dot0, 'y_dot0_mps', cs.y_dot0, 'z_dot0_mps', z_dot0, ...
        'drift_free_ydot0_mps', -2*n*x0, ...
        'note',           cs.label);

    data = struct( ...
        'metadata', metadata, ...
        't',        t, ...
        'X_t',      X_t', ...     % target ECI state, N x 6
        'Rho',      Rho');        % relative state in Hill frame, N x 6

    json_str = jsonencode(data, 'PrettyPrint', true);
    fname    = fullfile(out_dir, sprintf('rpo_nmc_%s.json', cs.key));
    fid      = fopen(fname, 'w');
    fprintf(fid, '%s', json_str);
    fclose(fid);
    fprintf('  exported -> %s\n\n', fname);
end

disp('Both scenarios exported. Load either into orbital_viz.m to compare.');
