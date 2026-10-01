% sim_adrios_full_mission.m
% Operator-in-the-loop sequential simulation of the full ADRIOS mission.
%
% Sequence: Closing (S4) -> NMC Flyaround (S2) -> Docking from FKP (S1)
%
% Before each phase the operator selects a controller configuration:
%   S4 - [A1] energy+ball  [B] dissipative  [C] energy+CLF (default)
%   S2 - [2]  periodic terminal equality (default)  [5] periodic stage cost
%   S1 - [1]  standard MPC  [2] EMPC+CLF (default)
%
% Each phase gate then pauses for a GO/HOLD decision before proceeding.
% Phase configs mirror the standalone sim_scenario{4,2,1}_*.m scripts.
% Docking enforces the LoS approach cone (30°, +V-bar) and y ≥ 0.

clear; clc;

%% Path setup
script_dir = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(script_dir, '..', '..', 'tools', 'matlab2tikz-master', 'src')));

%% Constants
constants;   % loads n, T, mu, etc.

%% Mission parameters
dt_S4    = 120;        % [s] closing sample time
dt_S2    = T / 92;    % [s] flyaround sample time
dt_S1    = 10;         % [s] docking sample time
N_S4     = 20;         % default MPC horizon, closing (B-noball uses 40)
P_S2     = 92;         % orbital period in steps
N_S2     = 92;         % MPC horizon, flyaround (= P)
N_S1     = 24;         % MPC horizon, docking (min. feasible for Run 3 term.eq @~150 m FKP)
u_max    = 1e-2;       % [m/s²] thrust bound, all phases
r_KOS    = 50;         % [m] keep-out sphere, closing phase
b_target = 75;         % [m] NMC orbit radial semi-axis target
los_deg  = 30;         % [deg] docking LoS approach-cone half-angle
rho_c    = 1e-11;      % state regularisation weight (B-noball option)
eps_r5   = 1e-4;       % augmented stage cost weight (Run 5 option)
P_avail  = 1e-5;       % [m²/s⁴] avg-power budget (tracking only)

%% Shared CWH continuous matrices
Ac_cwh = [0      0     0    1     0    0  ;
          0      0     0    0     1    0  ;
          0      0     0    0     0    1  ;
          3*n^2  0     0    0     2*n  0  ;
          0      0     0   -2*n   0    0  ;
          0      0    -n^2  0     0    0 ];
Bc_cwh = [zeros(3,3); eye(3)];

%% Precompute CLF terminal costs for S4 and S1
% Both are computed upfront regardless of which option is chosen later.

Q_dare = diag([1e-2, 1e-2, 1e-2, 1e0, 1e0, 1e0]);
R_dare = diag([1e4,  1e4,  1e4]);

[Ad_S4, Bd_S4] = discretize(Ac_cwh, Bc_cwh, dt_S4);
[~, P_clf_S4]  = lqr_controller(Ad_S4, Bd_S4, Q_dare, R_dare);

[Ad_S1, Bd_S1] = discretize(Ac_cwh, Bc_cwh, dt_S1);
[~, P_clf_S1]  = lqr_controller(Ad_S1, Bd_S1, Q_dare, R_dare);

%% =========================================================================
%%  PHASE 1 - Closing  (10 km -> 300 m)
%% =========================================================================

fprintf('\n%s\n', repmat('=', 1, 61));
fprintf('  PHASE 1 - Closing  (10 km -> 300 m)\n');
fprintf('%s\n', repmat('=', 1, 61));
fprintf('\n  Select controller:\n');
fprintf('    [1] Run A1  - pure energy  l=||u||2,  terminal ball ON  (r=200 m)\n');
fprintf('    [2] Run B   - dissipative  l=||u||2+rho_c||e||2,  no ball  (N=40)\n');
fprintf('    [3] Run C   - pure energy  l=||u||2,  CLF terminal cost  [default]\n\n');
s4_choice = select_option(3, 3);

% Base constraint fields shared by all S4 options
con_S4.u_max              = u_max;
con_S4.y_min_active       = false;
con_S4.los_cone.active    = false;
con_S4.use_avg_power      = false;
con_S4.avg_power.P_avail  = P_avail;
con_S4.use_passive_safety = false;

params_S4.x0       = [0; 10000; 0; 0; 0; 0];
params_S4.x_switch = [0;   300; 0; 0; 0; 0];
params_S4.dt       = dt_S4;
params_S4.u_max    = u_max;
params_S4.t_final  = 5 * T;
params_S4.r_KOS    = r_KOS;

switch s4_choice
    case 1   % Run A1: pure energy, terminal ball ON
        label_S4                          = 'Run A1 (energy + ball)';
        cost_S4                           = struct('Q', zeros(6), 'R', eye(3), 'P', zeros(6));
        con_S4.use_terminal_ball          = true;
        con_S4.terminal_ball.r_switch     = 200;   % [m]
        params_S4.N                       = N_S4;

    case 2   % Run B: strict dissipativity, no terminal ball, N=40
        label_S4                          = 'Run B (dissipative, N=40)';
        cost_S4                           = struct('Q', rho_c*eye(6), 'R', eye(3), 'P', zeros(6));
        con_S4.use_terminal_ball          = false;
        params_S4.N                       = 40;

    case 3   % Run C: pure energy + CLF terminal cost (default)
        label_S4                          = 'Run C (energy + CLF Vf)';
        cost_S4                           = struct('Q', zeros(6), 'R', eye(3), 'P', P_clf_S4);
        con_S4.use_terminal_ball          = false;
        params_S4.N                       = N_S4;
end

fprintf('\n  Running S4 [%s]...\n\n', label_S4);
res_S4 = run_closing(cost_S4, con_S4, params_S4);

% Extract terminal state and timing
if ~isnan(res_S4.conv_step)
    x_end_S4 = res_S4.X(:, res_S4.conv_step);
    steps_S4 = res_S4.conv_step;
    time_S4  = res_S4.t_vec(res_S4.conv_step) / 60;
else
    x_end_S4 = res_S4.X(:, end);
    steps_S4 = res_S4.n_ctrl;
    time_S4  = res_S4.n_ctrl * dt_S4 / 60;
    fprintf('\n*** WARNING: Phase 1 DID NOT CONVERGE in %d steps ***\n', res_S4.n_ctrl);
end

dv_S4      = res_S4.dv;
dv_mission = dv_S4;

% Gate 1
fprintf('\n');
fprintf('%s\n', repmat('=', 1, 60));
fprintf('PHASE 1 COMPLETE - CLOSING (10 km -> 300 m)\n');
fprintf('%s\n', repmat('=', 1, 60));
fprintf('Controller:    %s\n', label_S4);
if isnan(res_S4.conv_step)
fprintf('*** DID NOT CONVERGE ***\n');
end
fprintf('Deltav consumed:       %7.4f m/s\n', dv_S4);
fprintf('Steps / time:      %4d steps  (%5.1f min)\n', steps_S4, time_S4);
fprintf('Passive safety:    %5.1f%%\n', 100*res_S4.ps_frac);
fprintf('Terminal state:   [x,y,z] = [%6.2f,%7.2f,%5.2f] m\n', ...
    x_end_S4(1), x_end_S4(2), x_end_S4(3));
fprintf('Mission Deltav so far:  %7.4f m/s\n', dv_mission);
fprintf('%s\n\n', repmat('=', 1, 60));

while true
    choice = input('  -> Type GO to proceed to Flyaround, or HOLD to pause: ', 's');
    if strcmpi(strtrim(choice), 'GO'), break; end
    fprintf('  Holding at 300 m. Type GO to continue.\n');
end

%% =========================================================================
%%  PHASE 2 - NMC Flyaround acquisition
%% =========================================================================

fprintf('\n%s\n', repmat('=', 1, 61));
fprintf('  PHASE 2 - NMC Flyaround Acquisition\n');
fprintf('%s\n', repmat('=', 1, 61));
fprintf('\n  Select controller:\n');
fprintf('    [1] Run 2   - periodic terminal equality  x(N)=Pi*(k+N-1 mod P)  [default]\n');
fprintf('    [2] Run 5   - periodic augmented stage cost  l_aug=||u||2+eps||x-x*_k||2\n\n');
s2_choice = select_option(2, 1);

% Discretize at dt_S2 and verify NMC closure
[A_d_S2, B_d_S2] = discretize(Ac_cwh, Bc_cwh, dt_S2);
x_test_S2 = [100; 0; 0; 0; -2*n*100; 0];
closure_S2 = norm(A_d_S2^P_S2 * x_test_S2 - x_test_S2);
fprintf('\nNMC closure residual: %.2e  (expect < 1e-4)\n', closure_S2);
assert(closure_S2 < 1e-4, 'NMC orbit does not close: check dt_S2 and P_S2');

% Prescribed NMC orbit Pi* at b_target
n_x_S2 = 6;
n_u_S2 = 3;
x_star_S2 = zeros(n_x_S2, P_S2);
x_star_S2(:,1) = [b_target; 0; 0; 0; -2*n*b_target; 0];
for k_s = 1:P_S2-1
    x_star_S2(:, k_s+1) = A_d_S2 * x_star_S2(:, k_s);
end

% Condensed prediction matrices Sx [6N×6], Su [6N×3N]
Sx_S2  = zeros(n_x_S2*N_S2, n_x_S2);
Su_S2  = zeros(n_x_S2*N_S2, n_u_S2*N_S2);
Ad_pow = eye(n_x_S2);
for i = 1:N_S2
    Ad_pow = Ad_pow * A_d_S2;
    Sx_S2((i-1)*n_x_S2+1 : i*n_x_S2, :) = Ad_pow;
    for j = 1:i
        Su_S2((i-1)*n_x_S2+1 : i*n_x_S2, (j-1)*n_u_S2+1 : j*n_u_S2) = ...
            A_d_S2^(i-j) * B_d_S2;
    end
end

% Position-row extracts
pos_idx_S2 = zeros(1, 3*N_S2);
for ki = 1:N_S2
    pos_idx_S2((ki-1)*3+1 : ki*3) = (ki-1)*6+1 : (ki-1)*6+3;
end
Sx_pos_S2 = Sx_S2(pos_idx_S2, :);
Su_pos_S2 = Su_S2(pos_idx_S2, :);

% Base input-bound matrices (shared by both options)
A_u_S2 = [eye(3*N_S2); -eye(3*N_S2)];
b_u_S2 = u_max * ones(6*N_S2, 1);

% Base matrices struct (pure fuel H; Run 5 overwrites H, adds F_mat and f_ref_mat)
mats_S2.A_d        = A_d_S2;
mats_S2.B_d        = B_d_S2;
mats_S2.Sx         = Sx_S2;
mats_S2.Su         = Su_S2;
mats_S2.Sx_pos_all = Sx_pos_S2;
mats_S2.Su_pos_all = Su_pos_S2;
mats_S2.A_u        = A_u_S2;
mats_S2.b_u        = b_u_S2;

params_S2.x0      = x_end_S4;
params_S2.rho_min = 50;
params_S2.rho_max = 200;
params_S2.u_max   = u_max;
params_S2.N       = N_S2;
params_S2.P       = P_S2;
params_S2.n       = n;
params_S2.dt      = dt_S2;
params_S2.n_sim   = 6 * P_S2;

switch s2_choice
    case 1   % Run 2: periodic terminal equality (default)
        label_S2 = 'Run 2 (periodic terminal equality)';
        mats_S2.H                        = 2 * eye(3*N_S2);
        run_cfg_S2.use_nmc_manifold      = true;
        run_cfg_S2.use_periodic_terminal = true;
        run_cfg_S2.x_star                = x_star_S2;
        run_cfg_S2.n_sim_override        = [];
        run_cfg_S2.use_band              = false;

    case 2   % Run 5: periodic augmented stage cost (no terminal constraint)
        label_S2  = sprintf('Run 5 (periodic l_aug, eps=%.0e)', eps_r5);
        Qbar_r5   = kron(eye(N_S2), eps_r5 * eye(6));
        H_r5      = 2 * (Su_S2' * Qbar_r5 * Su_S2 + eye(3*N_S2));
        mats_S2.H          = (H_r5 + H_r5') / 2;
        mats_S2.F_mat      = 2 * Su_S2' * Qbar_r5 * Sx_S2;     % 3N × 6
        mats_S2.f_ref_mat  = -2 * Su_S2' * Qbar_r5;            % 3N × 6N
        run_cfg_S2.use_nmc_manifold      = false;
        run_cfg_S2.use_periodic_terminal = false;
        run_cfg_S2.x_star                = x_star_S2;   % needed for f_ref_mat evaluation
        run_cfg_S2.n_sim_override        = [];
        run_cfg_S2.use_band              = false;
end

fprintf('\n  Running S2 [%s]...\n\n', label_S2);
res_S2 = run_flyaround(run_cfg_S2, mats_S2, params_S2);

b_established = res_S2.b_est;
dv_S2         = res_S2.dv_total;

% Phase error over last orbit
n_sim_S2     = res_S2.n_sim;
phase_err_S2 = zeros(1, n_sim_S2 + 1);
for k_pe = 1:n_sim_S2 + 1
    col = mod(k_pe - 1, P_S2) + 1;
    phase_err_S2(k_pe) = norm(res_S2.x_log(:, k_pe) - x_star_S2(:, col));
end
last_orb  = max(1, n_sim_S2 - P_S2 + 1);
phase_err_max = max(phase_err_S2(last_orb:end));

% Find last V-bar crossing in x_log: search final half-orbit for min|x| with y>0.
% The NMC has vx=+nb at y=+2b; killing that velocity gives the FKP.
last_half_S2  = max(1, n_sim_S2 - round(P_S2/2));
x_srch_S2     = abs(res_S2.x_log(1, last_half_S2:end));
y_srch_S2     = res_S2.x_log(2, last_half_S2:end);
x_srch_S2(y_srch_S2 < b_established * 0.5) = inf;   % mask behind-target half
[~, k_rel_S2] = min(x_srch_S2);
k_vbar_S2     = last_half_S2 + k_rel_S2 - 1;         % column index in x_log
x_depart      = res_S2.x_log(:, k_vbar_S2);
fprintf('  V-bar crossing detected at step %d of %d (y=%.1f m, vx=%.4f m/s).\n', ...
    k_vbar_S2-1, n_sim_S2, x_depart(2), x_depart(4));

% Departure burn: cancel actual radial velocity at V-bar crossing
dv_departure = abs(x_depart(4));
x_FKP        = [0; x_depart(2); 0; 0; 0; 0];

dv_mission = dv_S4 + dv_S2 + dv_departure;

% Gate 2
fprintf('\n');
fprintf('%s\n', repmat('=', 1, 60));
fprintf('PHASE 2 COMPLETE - NMC FLYAROUND\n');
fprintf('%s\n', repmat('=', 1, 60));
fprintf('Controller:    %s\n', label_S2);
fprintf('Deltav consumed (EMPC):    %7.4f m/s\n', dv_S2);
fprintf('Orbit radius achieved:  %5.1f m  (target: 75.0 m)\n', b_established);
fprintf('Phase error (max, last orbit):  %6.2f m\n', phase_err_max);
fprintf('-- Departure burn --\n');
fprintf('DV departure:           %7.4f m/s  -> FKP\n', dv_departure);
fprintf('FKP state:              [0, %5.1f, 0] m (V-bar)\n', x_FKP(2));
fprintf('Mission Deltav so far:       %7.4f m/s\n', dv_mission);
fprintf('%s\n\n', repmat('=', 1, 60));

while true
    choice = input('  -> Type GO to proceed to Docking, or HOLD to pause: ', 's');
    if strcmpi(strtrim(choice), 'GO'), break; end
    fprintf('  Holding at FKP. Type GO to continue.\n');
end

%% =========================================================================
%%  PHASE 3 - Docking from FKP
%%  LoS approach cone (30°, +V-bar) + y >= 0 active (matches standalone S1).
%% =========================================================================

fprintf('\n%s\n', repmat('=', 1, 61));
fprintf('  PHASE 3 - Docking from FKP  (%.1f m -> 0 m)\n', x_FKP(2));
fprintf('%s\n', repmat('=', 1, 61));
fprintf('\n  Select controller:\n');
fprintf('    [1] Run 2   - EMPC, CLF terminal cost   Q=0, R=I, P=P_clf         [default]\n');
fprintf('    [2] Run 3   - EMPC, hard terminal eq.   Q=0, R=I, P=0, x(N)=dock\n');
fprintf('    (Both produce near-identical V-bar trajectories; mechanism differs.)\n\n');
s1_choice = select_option(2, 1);

% Base constraint struct, mirrors standalone sim_scenario1_docking.m:
% LoS approach cone (30°, +V-bar) + y ≥ 0 half-space; no approach corridor.
con_S1.u_max               = u_max;
con_S1.y_min_active        = true;
con_S1.los_cone.active     = true;
con_S1.los_cone.half_angle = deg2rad(los_deg);   % along-track approach cone
con_S1.use_passive_safety  = false;

params_S1.x0       = x_FKP;
params_S1.x_dock   = zeros(6, 1);
params_S1.dt       = dt_S1;
params_S1.N        = N_S1;
params_S1.u_max    = u_max;
params_S1.t_final  = 3 * T;
params_S1.r_KOS    = 5;
params_S1.conv_tol = 1.0;   % [m] position-error tol (matches standalone)

switch s1_choice
    case 1   % Run 2: EMPC, CLF terminal cost (default)
        label_S1 = 'Run 2 (EMPC, CLF Vf)';
        cost_S1  = struct('Q', zeros(6), 'R', eye(3), 'P', P_clf_S1);
        % terminal_eq not set, defaults to false in run_docking

    case 2   % Run 3: EMPC, hard terminal equality x(N) = x_dock
        label_S1          = 'Run 3 (EMPC, hard term. eq.)';
        cost_S1           = struct('Q', zeros(6), 'R', eye(3), 'P', zeros(6));
        con_S1.terminal_eq = true;
end

fprintf('\n  Running S1 [%s]...\n\n', label_S1);
res_S1 = run_docking(cost_S1, con_S1, params_S1);

if ~isnan(res_S1.conv_step)
    steps_S1 = res_S1.conv_step;
    time_S1  = res_S1.t_vec(res_S1.conv_step) / 60;
else
    steps_S1 = res_S1.n_ctrl;
    time_S1  = res_S1.n_ctrl * dt_S1 / 60;
    fprintf('\n*** WARNING: Phase 3 DID NOT CONVERGE in %d steps ***\n', res_S1.n_ctrl);
end

dv_S1            = res_S1.dv;
dv_total_mission = dv_S4 + dv_S2 + dv_departure + dv_S1;

% Gate 3 / Mission summary
fprintf('\n');
fprintf('%s\n', repmat('=', 1, 60));
fprintf('PHASE 3 COMPLETE - DOCKING APPROACH\n');
fprintf('%s\n', repmat('=', 1, 60));
fprintf('Controller:    %s\n', label_S1);
if isnan(res_S1.conv_step)
fprintf('*** DID NOT CONVERGE ***\n');
end
fprintf('Deltav consumed:       %7.4f m/s\n', dv_S1);
fprintf('Steps / time:      %4d steps  (%5.1f min)\n', steps_S1, time_S1);
fprintf('Terminal distance:  %6.2f m from docking port\n', res_S1.term_dist);
fprintf('%s\n', repmat('=', 1, 60));
fprintf('MISSION COMPLETE - ADRIOS FULL SEQUENCE\n');
fprintf('S4 [%-18s]:  %7.4f m/s\n', label_S4, dv_S4);
fprintf('S2 [%-18s]:  %7.4f m/s\n', label_S2, dv_S2);
fprintf('Departure burn:            %7.4f m/s\n', dv_departure);
fprintf('S1 [%-18s]:  %7.4f m/s\n', label_S1, dv_S1);
fprintf('%s\n', repmat('-', 1, 40));
fprintf('TOTAL Deltav:              %7.4f m/s\n', dv_total_mission);
fprintf('%s\n\n', repmat('=', 1, 60));

%% =========================================================================
%%  JSON EXPORT for orbital_viz.m
%% =========================================================================
%
% Format required by orbital_viz.m:
%   t     - time vector [N×1, s]
%   X_t   - target ECI state [N×6, m and m/s]
%   Rho   - chaser LVLH (Hill-frame) state [N×6, m and m/s]
%   metadata.period_min  (used for PERIOD telemetry in the viz)
%
% Phases are concatenated in sequence. The time vector is non-uniform
% (120 s per step in S4, ~60 s in S2, 10 s in S1).  The viz animation
% pace is set by t(2)-t(1) = dt_S4; use the speed slider at ~100-200x
% to play through the full mission in reasonable wall-clock time.

n_run_S4 = size(res_S4.X, 2);
n_run_S1 = size(res_S1.X, 2);

% -- Phase 1: full S4 trajectory --
Rho_S4 = res_S4.X;                          % 6 × n_run_S4
t_S4   = (0 : n_run_S4-1) * dt_S4;         % 1 × n_run_S4

% -- Phase 2: S2 up to V-bar departure crossing (k_vbar_S2) --
% Skip col 1 (= S4 terminal state, already the last frame of Rho_S4).
% Truncating at k_vbar_S2 removes the return leg to R-bar and eliminates
% the spatial jump in the visualization.
n_steps_S2 = k_vbar_S2 - 1;
Rho_S2 = res_S2.x_log(:, 2:k_vbar_S2);    % 6 × n_steps_S2
t_S2   = t_S4(end) + (1:n_steps_S2) * dt_S2;

% -- Departure burn: single frame at x_FKP --
t_burn  = t_S2(end) + dt_S2;               % one S2 step after flyaround ends
Rho_dep = x_FKP;                            % 6 × 1

% -- Phase 3: S1 (skip col 1 = x_FKP, already inserted as Rho_dep) --
if n_run_S1 > 1
    Rho_S1 = res_S1.X(:, 2:end);           % 6 × (n_run_S1-1)
    t_S1   = t_burn + (1 : n_run_S1-1) * dt_S1;
else
    Rho_S1 = zeros(6, 0);
    t_S1   = zeros(1, 0);
end

% -- Unified arrays --
Rho_all = [Rho_S4, Rho_S2, Rho_dep, Rho_S1];  % 6 × N_total
t_all   = [t_S4,   t_S2,   t_burn,  t_S1];     % 1 × N_total
N_all   = length(t_all);

% Phase index per frame (1=closing, 2=flyaround, 3=departure burn, 4=docking)
phase_idx = [ones(1, n_run_S4), ...
             2*ones(1, n_steps_S2), ...
             3, ...
             4*ones(1, max(0, n_run_S1-1))];

% -- Target ECI trajectory --
% Circular Keplerian orbit (consistent with CWH assumption):
%   RAAN=0, AoP=0, i=51.6 deg (ISS), starting at ascending node.
X_t_all = zeros(N_all, 6);
for k_viz = 1:N_all
    theta   = n * t_all(k_viz);   % true anomaly [rad]
    pos_eci = a * [cos(theta); sin(theta)*cos(i); sin(theta)*sin(i)];
    vel_eci = a*n * [-sin(theta); cos(theta)*cos(i); cos(theta)*sin(i)];
    X_t_all(k_viz, :) = [pos_eci; vel_eci]';
end

% -- Build and write JSON --
viz_data.metadata = struct( ...
    'mission',         'ADRIOS Full Sequence', ...
    'date',            char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')), ...
    'a_km',            a / 1000, ...
    'e',               0, ...
    'i_deg',           rad2deg(i), ...
    'period_min',      T / 60, ...
    'altitude_km',     alt / 1000, ...
    'controller_S4',   label_S4, ...
    'controller_S2',   label_S2, ...
    'controller_S1',   label_S1, ...
    'dv_S4_ms',        dv_S4, ...
    'dv_S2_ms',        dv_S2, ...
    'dv_departure_ms', dv_departure, ...
    'dv_S1_ms',        dv_S1, ...
    'dv_total_ms',     dv_total_mission, ...
    'b_nmc_m',         b_established, ...
    'los_cone_deg',    los_deg, ...
    'n_frames',        N_all, ...
    'dt_hint_s',       dt_S4);   % animation pacing hint (first step)

viz_data.t     = t_all';       % N × 1
viz_data.X_t   = X_t_all;     % N × 6  (target ECI)
viz_data.Rho   = Rho_all';    % N × 6  (chaser LVLH / Hill frame)
viz_data.phase = phase_idx';   % N × 1  (phase tag per frame)

json_str = jsonencode(viz_data, 'PrettyPrint', true);
out_path = fullfile(script_dir, '..', '..', 'exports', 'scenarios', 'sim_adrios_full_mission.json');
fid = fopen(out_path, 'w');
fprintf(fid, '%s', json_str);
fclose(fid);
fprintf('Viz JSON exported to exports/scenarios/sim_adrios_full_mission.json\n');
fprintf('Open orbital_viz.m, click LOAD JSON, and set speed ~100-200x.\n\n');

%% =========================================================================
%%  Local functions
%% =========================================================================

function choice = select_option(n_options, default_val)
% Reads a menu selection; returns default_val on empty (Enter) input.
while true
    raw = strtrim(input(sprintf('  Enter choice (1-%d, default=%d): ', ...
        n_options, default_val), 's'));
    if isempty(raw)
        choice = default_val;
        return;
    end
    val = str2double(raw);
    if ~isnan(val) && val >= 1 && val <= n_options && val == floor(val)
        choice = val;
        return;
    end
    fprintf('  Invalid: enter a number between 1 and %d.\n', n_options);
end
end