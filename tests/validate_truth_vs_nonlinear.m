% validate_truth_vs_nonlinear.m
%
% REGRESSION ANCHOR for the Aerospace-Toolbox truth model.
%
% With all perturbations OFF, the absolute truth propagator is exact two-body
% motion. Differencing two such trajectories into the Hill frame (eci2hill.m)
% must therefore reproduce the EXACT nonlinear Keplerian relative dynamics
% (relative_motion_nonlinear.m) -- same physics, two independent code paths.
%
% This test drives both paths from a single ode113 call over a combined
% 18-state so they share the identical target trajectory and time grid:
%
%   state = [ X_t (6) ; X_c (6) ; Rho (6) ]
%     X_t, X_c : absolute ECI, propagated by truth_propagator_at (pert OFF)
%     Rho      : relative Hill, propagated by relative_motion_nonlinear
%
% PASS: max ||eci2hill(X_t,X_c) - Rho|| stays at integration-tolerance level.
% If it does not, the bug is in eci2hill/hill2eci -- caught BEFORE any sweep
% result is trusted.
%
% This must pass before truth_propagator_at is used with perturbations ON.

clear; clc; close all;

addpath('src/dynamics'); addpath('src/utils');
constants;

fprintf('====================================================\n');
fprintf(' Truth-model validation: eci2hill vs nonlinear EoM\n');
fprintf('====================================================\n\n');

%% Target: elliptic orbit (exercises radial rate + varying omega)
e_t = 0.2;
a_t = a;
r_p = a_t*(1 - e_t);                          % perigee radius
v_p = sqrt(mu*(1 + e_t)/(a_t*(1 - e_t)));     % perigee speed
X_t0 = [r_p; 0; 0; 0; v_p; 0];

%% Chaser: seeded from a Hill-frame relative IC (100s of metres + cross-track)
Rho0 = [200; -350; 120; 0.10; -0.05; 0.08];   % [m, m/s]
X_c0 = hill2eci(X_t0, Rho0);

% Consistency: eci2hill of that chaser must return Rho0 exactly
Rho0_check = eci2hill(X_t0, X_c0);
fprintf('Seed consistency (hill2eci->eci2hill) : %.3e m\n', ...
        max(abs(Rho0_check - Rho0)));

%% Combined 18-state propagation
p = struct('mu', mu);          % perturbations default OFF -> pure two-body
p.use_zonal = false;
p.use_drag  = false;

odefun = @(t, s) [ truth_propagator_at(t, s(1:6),  p);      % target absolute
                   truth_propagator_at(t, s(7:12), p);      % chaser absolute
                   relative_motion_nonlinear(s(13:18), s(1:6), mu) ]; % relative

s0     = [X_t0; X_c0; Rho0];
t_span = [0, 3*T];             % three target orbital periods
opts   = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

[t, S] = ode113(odefun, t_span, s0, opts);

%% Compare differenced truth against the directly-integrated relative model
N = numel(t);
err = zeros(N, 1);
Rho_diff = zeros(N, 6);
for k = 1:N
    Rho_diff(k, :) = eci2hill(S(k, 1:6).', S(k, 7:12).').';
    err(k) = norm(Rho_diff(k, :).' - S(k, 13:18).');
end

pos_err = vecnorm(Rho_diff(:, 1:3) - S(:, 13:15), 2, 2);
vel_err = vecnorm(Rho_diff(:, 4:6) - S(:, 16:18), 2, 2);

fprintf('\nOver %.1f orbits (%d output steps):\n', 3, N);
fprintf('  max position error : %.3e m\n',   max(pos_err));
fprintf('  max velocity error : %.3e m/s\n', max(vel_err));

%% Verdict
tol_pos = 1e-3;    % 1 mm over 3 orbits -- generous vs 1e-12 integrator tol
if max(pos_err) < tol_pos
    fprintf('\n  PASS: differencing reproduces nonlinear EoM.\n');
else
    fprintf(2, '\n  FAIL: eci2hill/hill2eci inconsistent with EoM.\n');
end

%% Error growth (informative)
figure;
subplot(2,1,1);
semilogy(t/T, pos_err, 'LineWidth', 1.2); grid on;
xlabel('orbits'); ylabel('position error [m]');
title('Truth differencing vs nonlinear relative EoM');
subplot(2,1,2);
semilogy(t/T, vel_err, 'LineWidth', 1.2); grid on;
xlabel('orbits'); ylabel('velocity error [m/s]');
