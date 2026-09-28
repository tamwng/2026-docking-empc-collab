% sim_scenario4_closing.m
% Scenario 4: Three EMPC stabilisation mechanisms on the closing phase.
%
% Run A1     - pure fuel (Q=0, P=0), terminal ball ON: converges.
% Run A2     - same, ball OFF: u≡0, chaser stuck
% Run B-noball - fuel + small state penalty ℓ=‖u‖²+ρc‖e‖², no ball: converges
%               via strict dissipativity
% Run B-ball - same cost, ball ON: confirms ball is inactive.
% Run C      - pure fuel, CLF terminal cost Vf=e'P_clf e, no ball: asymptotically
%               stable without strict dissipativity (Amrit et al. 2011).
% Run MPC    - standard regulation MPC (Q,R≠0, LQR terminal cost): fuel-hungry baseline.

clear; clc;

%% Path setup
script_dir = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(script_dir, '..', '..', 'tools', 'matlab2tikz-master', 'src')));

%% Constants
constants;

%% Shared scenario parameters
dt        = 120;
N         = 20;
x0        = [0; 10000; 0; 0; 0; 0];
x_switch  = [0;   300; 0; 0; 0; 0];  % [m] handover point (V-bar, 300 m)
u_max     = 1e-2;                     % [m/s²] thrust bound
r_KOS     = 50;                       % [m] keep-out sphere radius
P_avail   = 1e-5;                     % [m²/s⁴] avg-power budget
rho_c     = 1e-11;                     % state regularisation weight
t_final   = 5 * T;                    % [s] ~7.7 h (5 orbital periods)

%% Cost structs
cost_A = struct('Q', zeros(6),    'R', eye(3), 'P', zeros(6));
cost_B = struct('Q', rho_c*eye(6),'R', eye(3), 'P', zeros(6));

%% Precompute CLF terminal cost (Run C)
Ac_cwh = [0      0     0    1     0    0  ;
          0      0     0    0     1    0  ;
          0      0     0    0     0    1  ;
          3*n^2  0     0    0     2*n  0  ;
          0      0     0   -2*n   0    0  ;
          0      0    -n^2  0     0    0 ];
Bc_cwh = [zeros(3,3); eye(3)];
[Ad_clf, Bd_clf] = discretize(Ac_cwh, Bc_cwh, dt);

Q_dare = diag([1e-2, 1e-2, 1e-2, 1e0, 1e0, 1e0]);
R_dare = diag([1e4,  1e4,  1e4]);
[K_clf, P_clf] = lqr_controller(Ad_clf, Bd_clf, Q_dare, R_dare);
R_eco_clf      = eye(3);

% Verify Amrit et al. (2011) Assumption 6 at dt=120s.
% M1 = P - Acl'*P*Acl - K'*R_eco*K = Q_dare + K'*(R_dare - R_eco)*K >= 0.
Acl_clf = Ad_clf - Bd_clf * K_clf;
M1      = Q_dare + K_clf' * (R_dare - R_eco_clf) * K_clf;
M1      = (M1 + M1') / 2;
lam_M1  = min(eig(M1));
rho_Acl = max(abs(eig(Acl_clf)));
if lam_M1 > 1e-10 && rho_Acl < 1
    fprintf('  => P_clf satisfies Assumption 6. Terminal cost is a valid CLF.  ✓\n\n');
else
    warning('CLF condition not satisfied: check DARE weights.');
end

cost_C = struct('Q', zeros(6), 'R', eye(3), 'P', P_clf);

%% Shared params base
params_base.x0       = x0;
params_base.x_switch = x_switch;
params_base.dt       = dt;
params_base.N        = N;
params_base.u_max    = u_max;
params_base.t_final  = t_final;
params_base.r_KOS    = r_KOS;

%% Run A1: pure energy, terminal ball ON

con_A1.u_max                  = u_max;
con_A1.y_min_active           = false;
con_A1.los_cone.active        = false;
con_A1.use_avg_power          = false;
con_A1.avg_power.P_avail      = P_avail;
con_A1.use_terminal_ball      = true;
con_A1.terminal_ball.r_switch = 200;   % [m] L-inf position ball
con_A1.use_passive_safety     = false;

params_A1 = params_base;

fprintf(' Run A1 - pure energy ℓ=‖u‖², ball ON\n');
resA1 = run_closing(cost_A, con_A1, params_A1);

%% Run A2: pure energy, terminal ball OFF

con_A2.u_max              = u_max;
con_A2.y_min_active       = false;
con_A2.los_cone.active    = false;
con_A2.use_avg_power      = false;
con_A2.avg_power.P_avail  = P_avail;
con_A2.use_terminal_ball  = false;
con_A2.use_passive_safety = false;

params_A2                  = params_base;
params_A2.n_steps_override = 60;

fprintf(' Run A2 - pure energy ℓ=‖u‖², ball OFF\n');
resA2 = run_closing(cost_A, con_A2, params_A2);

max_disp_A2 = max(vecnorm(resA2.X(1:3,:) - x0(1:3), 2, 1));
max_u_A2    = max([resA2.u_norm_seq, 0]);
fprintf('A2 check: max position displacement: %.4f m   (expect ≈ 0)\n', max_disp_A2);
fprintf('          max ||u||:                 %.2e m/s² (expect ≈ 0)\n\n', max_u_A2);

%% Run B-noball: strict dissipativity, no terminal ball
% Adding ρc‖e‖² restores strict dissipativity (storage fn λ≡0, α_ℓ(r)=ρc·r²).

N_B_noball = N;

con_B_noball.u_max              = u_max;
con_B_noball.y_min_active       = false;
con_B_noball.los_cone.active    = false;
con_B_noball.use_avg_power      = false;
con_B_noball.avg_power.P_avail  = P_avail;
con_B_noball.use_terminal_ball  = false;
con_B_noball.use_passive_safety = false;

params_B_noball   = params_base;
params_B_noball.N = N_B_noball;

fprintf(' Run B-noball - ρc=%.0e, N=%d, ball OFF\n', rho_c, N_B_noball);
resB_noball = run_closing(cost_B, con_B_noball, params_B_noball);

if isnan(resB_noball.conv_step)
    N_B_noball        = 40;
    params_B_noball.N = N_B_noball;
    fprintf('\n[Retry] N=20 did not converge; retrying with N=%d...\n\n', N_B_noball);
    fprintf('%s\n', repmat('=', 1, 48));
    fprintf(' Run B-noball (retry) - N=%d, ball OFF\n', N_B_noball);
    fprintf('%s\n', repmat('=', 1, 48));
    resB_noball = run_closing(cost_B, con_B_noball, params_B_noball);
end
fprintf('B-noball: used N=%d, converged=%s\n', N_B_noball, ...
    mat2str(~isnan(resB_noball.conv_step)));

%% Run B-ball: strict dissipativity + terminal ball ON
% terminal constraint never binds.  Active-fraction ≈ 0 confirms this.

r_switch_ball_B = 200;

con_B_ball = con_B_noball;
con_B_ball.use_terminal_ball      = true;
con_B_ball.terminal_ball.r_switch = r_switch_ball_B;

params_B_ball = params_base;   % N=20

fprintf(' Run B-ball - ρc=%.0e, N=20, ball ON r=%d m\n', rho_c, r_switch_ball_B);
resB_ball = run_closing(cost_B, con_B_ball, params_B_ball);

%% Run C: pure energy + CLF terminal cost, no ball
% Terminal cost Vf = e'P_clf e (P_clf from DARE) is a valid local CLF for
% the error dynamics.

con_C.u_max              = u_max;
con_C.y_min_active       = false;
con_C.los_cone.active    = false;
con_C.use_avg_power      = true;
con_C.avg_power.P_avail  = P_avail;
con_C.use_terminal_ball  = false;
con_C.use_passive_safety = false;

params_C = params_base;


fprintf(' Run C - energy ℓ=‖u‖², CLF terminal cost V_f=e''P_clf e\n');
resC = run_closing(cost_C, con_C, params_C);

%% Run MPC: Standard regulation MPC (comparison baseline)

cost_std = struct('Q', Q_dare, 'R', R_dare, 'P', P_clf);

con_MPC.u_max              = u_max;
con_MPC.y_min_active       = false;
con_MPC.los_cone.active    = false;
con_MPC.use_avg_power      = true;
con_MPC.avg_power.P_avail  = P_avail;
con_MPC.use_terminal_ball  = false;
con_MPC.use_passive_safety = false;

params_MPC = params_base;

fprintf(' Run MPC - standard regulation ℓ=e''Qe+u''Ru  (Q,R = DARE weights, P = LQR cost-to-go)\n');
resMPC = run_closing(cost_std, con_MPC, params_MPC);


%% Postprocessing

% Ball-role evidence (trajectory comparison)

n_B_cmp  = min(resB_noball.n_ctrl, resB_ball.n_ctrl);
dX_B_max = max(abs(resB_noball.X(1:3, 1:n_B_cmp+1) - resB_ball.X(1:3, 1:n_B_cmp+1)), [], 'all');

fprintf('\nBall-role evidence:\n');
fprintf('  A2 total Δv (no ball, Q=P=0):     %.2e m/s  -> ball is load-bearing in A1\n', resA2.dv);
fprintf('  ||pos(B-noball) - pos(B-ball)||:  %.2e m   -> ball is slack in B-ball\n\n', dX_B_max);

% MPC comparison: reuse the first-class Run MPC (resMPC) computed above
X_mpc  = resMPC.X;
U_mpc  = resMPC.U;
t_mpc  = resMPC.t_vec;
dv_mpc = resMPC.dv;

fprintf('\n%s\n', repmat('-', 1, 63));
fprintf('Δv summary\n');
fprintf('%s\n', repmat('-', 1, 63));
fprintf('  %-32s  %10s\n', 'Run', 'Δv [m/s]');
fprintf('  %-32s  %10.4f\n', 'A1 - energy, ball ON',               resA1.dv);
fprintf('  %-32s  %10.4f\n', sprintf('B-noball (N=%d)',N_B_noball),    resB_noball.dv);
fprintf('  %-32s  %10.4f\n', 'B-ball  - strict dissip.',              resB_ball.dv);
fprintf('  %-32s  %10.4f\n', 'C  - energy + CLF terminal cost',       resC.dv);
fprintf('  %-32s  %10.4f\n', 'MPC - standard regulation',             dv_mpc);
fprintf('  A1 vs MPC savings: %.4f m/s  (%.1f%%)\n', ...
    dv_mpc - resA1.dv, 100*(dv_mpc - resA1.dv)/dv_mpc);
fprintf('  C  vs MPC savings: %.4f m/s  (%.1f%%)\n', ...
    dv_mpc - resC.dv,  100*(dv_mpc - resC.dv) /dv_mpc);
fprintf('%s\n\n', repmat('-', 1, 63));

%% Figure: Hill-frame - A1, B-noball, C, MPC (full view + switch-point zoom)
c_A1  = [0.85 0.33 0.10];
c_Bnb = [0.20 0.63 0.17];
c_C   = [0.00 0.68 0.68];
c_mpc = [0.12 0.47 0.71];

conv_tol_s4 = 50;         % [m] run_closing convergence threshold
zoom_hw     = 150;        % [m] half-width of the switch-point zoom window


lbl_A1  = 'Term.\ ball';
lbl_Bnb = '$\rho_c$ stage-cost reg.';
lbl_C   = 'Term.\ cost';
lbl_MPC = 'Std.\ MPC';

fig4 = figure('Name', 'S4 - Hill Frame: A1, B-noball, C, MPC (full + zoom)');
fig4.Position(3:4) = [980, 420];

% Panel A: full trajectory (IC -> switch), carries the one shared legend
ax4a = subplot(1, 2, 1);
hold(ax4a,'on'); grid(ax4a,'on'); axis(ax4a,'equal');
plot(ax4a, resA1.X(2,:),       resA1.X(1,:),       '-', 'Color', c_A1, 'LineWidth', 1.4, ...
    'DisplayName', lbl_A1);
plot(ax4a, resB_noball.X(2,:), resB_noball.X(1,:), '-', 'Color', c_Bnb, 'LineWidth', 1.4, ...
    'DisplayName', lbl_Bnb);
plot(ax4a, resC.X(2,:),        resC.X(1,:),        '-', 'Color', c_C,  'LineWidth', 1.4, ...
    'DisplayName', lbl_C);
plot(ax4a, resMPC.X(2,:),      resMPC.X(1,:),      '-', 'Color', c_mpc, 'LineWidth', 1.4, ...
    'DisplayName', lbl_MPC);
plot(ax4a, x0(2), x0(1), 'o', 'Color', [0.4 0.4 0.4], ...
    'MarkerFaceColor', [0.4 0.4 0.4], 'MarkerSize', 6, 'DisplayName', 'IC (10 km)');
plot(ax4a, x_switch(2), x_switch(1), 's', 'Color', [0.2 0.6 0.2], ...
    'MarkerFaceColor', [0.2 0.6 0.2], 'MarkerSize', 7, 'DisplayName', 'Switch point');
plot(ax4a, 0, 0, '+k', 'MarkerSize', 9, 'LineWidth', 1.5, 'DisplayName', 'Target');
% Dashed rectangle showing the zoom window
rectangle(ax4a, 'Position', [x_switch(2)-zoom_hw, x_switch(1)-zoom_hw, 2*zoom_hw, 2*zoom_hw], ...
    'EdgeColor', [0.4 0.4 0.4], 'LineStyle', '--', 'LineWidth', 0.8);
xlabel(ax4a, 'Along-track $y$ [m]', 'Interpreter', 'latex', 'FontSize', 11);
ylabel(ax4a, 'Radial $x$ [m]',      'Interpreter', 'latex', 'FontSize', 11);
legend(ax4a, 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 8);

% Panel B: zoom on the switch point
ax4b = subplot(1, 2, 2);
hold(ax4b,'on'); grid(ax4b,'on'); axis(ax4b,'equal');
plot(ax4b, resA1.X(2,:),       resA1.X(1,:),       '-', 'Color', c_A1, 'LineWidth', 1.4);
plot(ax4b, resB_noball.X(2,:), resB_noball.X(1,:), '-', 'Color', c_Bnb, 'LineWidth', 1.4);
plot(ax4b, resC.X(2,:),        resC.X(1,:),        '-', 'Color', c_C,  'LineWidth', 1.4);
plot(ax4b, resMPC.X(2,:),      resMPC.X(1,:),      '-', 'Color', c_mpc, 'LineWidth', 1.4);
plot(ax4b, x_switch(2), x_switch(1), 's', 'Color', [0.2 0.6 0.2], ...
    'MarkerFaceColor', [0.2 0.6 0.2], 'MarkerSize', 8);
xlim(ax4b, x_switch(2) + [-zoom_hw, zoom_hw]);
ylim(ax4b, x_switch(1) + [-zoom_hw, zoom_hw]);
xlabel(ax4b, 'Along-track $y$ [m]', 'Interpreter', 'latex', 'FontSize', 11);
ylabel(ax4b, 'Radial $x$ [m]',      'Interpreter', 'latex', 'FontSize', 11);

out_dir_hill = 'results/figures/scenario4';
[~, ~] = mkdir(out_dir_hill);
if exist('matlab2tikz', 'file')
    matlab2tikz(fullfile(out_dir_hill, 'B_hill_frame.tikz'), ...
        'figurehandle', fig4, 'showInfo', false);
end
exportgraphics(fig4, fullfile(out_dir_hill, 'B_hill_frame_slide.png'), 'Resolution', 300);

%% JSON export
json_data.runA1 = make_run_entry(resA1, 'A1', 'pure_energy', true, 0, N, ...
    dt, x0, x_switch, u_max, con_A1.terminal_ball.r_switch, false, resA2.dv, true);

json_data.runA2 = make_run_entry(resA2, 'A2', 'pure_energy', false, 0, N, ...
    dt, x0, x_switch, u_max, NaN, false, NaN, true);

json_data.runB_noball = make_run_entry(resB_noball, 'B_noball', ...
    sprintf('energy_rho_c_%.0e', rho_c), false, rho_c, N_B_noball, ...
    dt, x0, x_switch, u_max, NaN, false, NaN, true);

json_data.runB_ball = make_run_entry(resB_ball, 'B_ball', ...
    sprintf('energy_rho_c_%.0e', rho_c), true, rho_c, N, ...
    dt, x0, x_switch, u_max, r_switch_ball_B, false, dX_B_max, true);

json_data.runC = make_run_entry(resC, 'C', 'energy_clf_terminal', false, 0, N, ...
    dt, x0, x_switch, u_max, NaN, false, NaN, false);

json_data.runMPC = make_run_entry(resMPC, 'MPC', 'standard_regulation_Qdare_Rdare', ...
    false, 0, N, dt, x0, x_switch, u_max, NaN, false, NaN, false);

json_str = jsonencode(json_data, 'PrettyPrint', true);
fid = fopen('exports/scenarios/sim_scenario4_closing.json', 'w');
fprintf(fid, '%s', json_str);
fclose(fid);
fprintf('JSON exported to exports/scenarios/sim_scenario4_closing.json\n');

%%  Local functions

function s = make_run_entry(res, run_label, cost_label, ball_on, rho_c_val, N_used, ...
                             dt_val, x0_val, x_switch_val, u_max_val, r_ball_val, ...
                             use_ps, ball_evidence, v_f_zero)
    if nargin < 14, v_f_zero = true; end
    s.metadata = struct( ...
        'run',                run_label, ...
        'cost',               cost_label, ...
        'ball_on',            ball_on, ...
        'rho_c',              rho_c_val, ...
        'N',                  N_used, ...
        'dt',                 dt_val, ...
        'x0',                 x0_val', ...
        'x_switch',           x_switch_val', ...
        'u_max',              u_max_val, ...
        'r_ball',             r_ball_val, ...
        'use_passive_safety', use_ps, ...
        'V_f_zero',           v_f_zero);
    s.results = struct( ...
        'total_dv_ms',            res.dv, ...
        'n_steps_to_convergence', res.conv_step, ...
        'passive_safe_fraction',  res.ps_frac, ...
        'terminal_distance_m',    res.term_dist, ...
        'ball_evidence',          ball_evidence, ...
        'n_ps_rows_avg',          mean(res.n_ps_h), ...
        'n_ps_rows_max',          max([res.n_ps_h, 0]));
    s.turnpike = struct( ...
        'u_thresh',          1e-4 * u_max_val, ...
        'turnpike_steps',    res.turnpike_steps, ...
        'turnpike_fraction', res.turnpike_fraction, ...
        'tp_label',          res.tp_label, ...
        'u_norm_sequence',   res.u_norm_seq);
    s.t       = res.t_vec;
    s.X       = res.X';
    s.U       = res.U';
    s.ef      = res.ef;
    s.is_safe = res.is_safe;
end
