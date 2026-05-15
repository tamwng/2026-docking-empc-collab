% sim_empc_regulation.m
% Economic MPC regulation to origin — comparison against standard regulation MPC.
%
%   Structure:
%   min  Σᵢ₌₀ᴺ⁻¹ uᵢ'R_eco uᵢ  +  xₙ'P xₙ
%
% We still need a terminal cost P. Ideally P would be the infinite-horizon
% economic cost-to-go V∞(x) = min Σᵢ₌₀^∞ uᵢ'R_eco uᵢ, but computing
% this requires solving an LQR/DARE with Q = 0 in the stage cost.
%
% Fix: compute P from DARE with Q_dare > 0, R_dare > 0 — the same matrices
% used in the standard regulation MPC.
%
% See Amrit, Rawlings & Angeli
% (2011) for a terminal cost that guarantees practical stability of EMPC.

clear; clc;


%% Constants
constants;

%% Simulation parameters

dt      = 60;           % [s]   sample time
t_final = 2 * T;        % [s]  
n_steps = round(t_final / dt);
t       = (0:n_steps-1) * dt;
t_min   = t / 60;

%% MPC parameters
N = 20;   % [-]  prediction horizon

%%  CONSTRAINTS 
con.u_max        = 1e-2;   % [m/s^2]
con.y_min_active = false;
con.los_cone.active     = false;
con.los_cone.half_angle = deg2rad(20);
con.los_cone.n_faces    = 10;

%%  CWH STATE SPACE
Ac = [0      0     0    1     0    0  ;
      0      0     0    0     1    0  ;
      0      0     0    0     0    1  ;
      3*n^2  0     0    0     2*n  0  ;
      0      0     0   -2*n   0    0  ;
      0      0    -n^2  0     0    0 ];
Bc = [zeros(3,3); eye(3)];

[Ad, Bd] = discretize(Ac, Bc, dt);

%%  INITIAL CONDITIONS
X0 = [500; 10000; 200; 0; 0; 0];


%%  TERMINAL COST 
%
%  Q_dare and R_dare are used ONLY for the DARE to compute P.
%  They are NOT the stage cost of the economic MPC.

Q_dare = diag([1e-2, 1e-2, 1e-2, 1e0, 1e0, 1e0]);
R_dare = diag([1e4,  1e4,  1e4]);

[~, P] = lqr_controller(Ad, Bd, Q_dare, R_dare);


%%  STANDARD REGULATION MPC (as before)

fprintf('Running Standard Regulation MPC...\n');

X_std = zeros(6, n_steps);   U_std = zeros(3, n_steps);
X_std(:,1) = X0;

for k = 1:n_steps-1
    [u_opt, ~, ~]   = mpc_regulation(X_std(:,k), Ad, Bd, Q_dare, R_dare, P, N, con);
    U_std(:,k)      = u_opt;
    X_std(:,k+1)    = Ad * X_std(:,k) + Bd * u_opt;
end
U_std(:,end) = mpc_regulation(X_std(:,end), Ad, Bd, Q_dare, R_dare, P, N, con);


%%  ECONOMIC MPC 
%   Stage cost: ℓ(x,u) = u'R_eco u  (fuel only — Q_stage = 0)
%   Terminal cost: x'P x  (same P as above — valid CLF, see header)


Q_stage = zeros(6);                     % no state penalty at intermediate steps
R_eco = eye(3);
%R_eco   = diag([1e3, 1e3, 1e3]);        % equal fuel weight per axis, tune for fuel&performance trade-off!!

fprintf('Running Economic MPC...\n');

X_eco = zeros(6, n_steps);   U_eco = zeros(3, n_steps);
X_eco(:,1) = X0;

for k = 1:n_steps-1
    [u_opt, ~, ~]   = mpc_regulation(X_eco(:,k), Ad, Bd, Q_stage, R_eco, P, N, con);
    U_eco(:,k)      = u_opt;
    X_eco(:,k+1)    = Ad * X_eco(:,k) + Bd * u_opt;
end
U_eco(:,end) = mpc_regulation(X_eco(:,end), Ad, Bd, Q_stage, R_eco, P, N, con);

fprintf('Simulation complete.\n');


%%  POST-PROCESSING

rel_std  = vecnorm(X_std(1:3,:), 2, 1);
rel_eco  = vecnorm(X_eco(1:3,:), 2, 1);
dv_std   = cumsum(vecnorm(U_std, 2, 1) * dt);
dv_eco   = cumsum(vecnorm(U_eco, 2, 1) * dt);

V_std = arrayfun(@(k) X_std(:,k)'*P*X_std(:,k), 1:n_steps);
V_eco = arrayfun(@(k) X_eco(:,k)'*P*X_eco(:,k), 1:n_steps);

tol = 0.01;
conv_std = find(rel_std < tol, 1, 'first');
conv_eco = find(rel_eco < tol, 1, 'first');

fprintf('\n--- Standard MPC ---\n');
fprintf('  Total delta-v: %.4f m/s\n', dv_std(end));
if ~isempty(conv_std)
    fprintf('  Converged at:  %.1f min\n', t_min(conv_std));
else
    fprintf('  Did not converge within simulation time.\n');
end

fprintf('\n--- Economic MPC ---\n');
fprintf('  Total delta-v: %.4f m/s\n', dv_eco(end));
if ~isempty(conv_eco)
    fprintf('  Converged at:  %.1f min\n', t_min(conv_eco));
else
    fprintf('  Did not converge within simulation time.\n');
end

fprintf('\nDelta-v saving: %.2f%%\n', 100*(dv_std(end)-dv_eco(end))/dv_std(end));


%%  PLOTS

c_std = [0.00 0.45 0.70];   % blue  — standard MPC
c_eco = [0.85 0.33 0.10];   % orange — economic MPC
c_x   = [0.00 0.60 0.90];
c_y   = [0.90 0.40 0.00];
c_z   = [0.20 0.80 0.20];

fig1 = figure('Name', 'EMPC Regulation — Comparison');

%% Relative distance
ax1 = subplot(2,3,1);
hold(ax1,'on'); grid(ax1,'on');
plot(ax1, t_min, rel_std, 'Color', c_std, 'LineWidth', 1.4, 'DisplayName', 'Std MPC');
plot(ax1, t_min, rel_eco, 'Color', c_eco, 'LineWidth', 1.4, 'DisplayName', 'EMPC');
if ~isempty(conv_std)
    xline(ax1, t_min(conv_std), '--', 'Color', c_std, 'LineWidth', 0.8);
end
if ~isempty(conv_eco)
    xline(ax1, t_min(conv_eco), '--', 'Color', c_eco, 'LineWidth', 0.8);
end
xlabel(ax1,'Time [min]'); ylabel(ax1,'Distance [m]');
title(ax1,'Relative Distance to Origin');
legend(ax1,'Location','northeast');

%% Cumulative delta-v
ax2 = subplot(2,3,2);
hold(ax2,'on'); grid(ax2,'on');
plot(ax2, t_min, dv_std, 'Color', c_std, 'LineWidth', 1.4, 'DisplayName', 'Std MPC');
plot(ax2, t_min, dv_eco, 'Color', c_eco, 'LineWidth', 1.4, 'DisplayName', 'EMPC');
xlabel(ax2,'Time [min]');
ylabel(ax2,'$\Delta v$ [m/s]','Interpreter','latex');
title(ax2,'Cumulative $\Delta v$','Interpreter','latex');
legend(ax2,'Location','northwest');

%% Lyapunov function
ax3 = subplot(2,3,3);
hold(ax3,'on'); grid(ax3,'on');
semilogy(ax3, t_min, V_std, 'Color', c_std, 'LineWidth', 1.4, 'DisplayName', 'Std MPC');
semilogy(ax3, t_min, V_eco, 'Color', c_eco, 'LineWidth', 1.4, 'DisplayName', 'EMPC');
xlabel(ax3,'Time [min]');
ylabel(ax3,'$V = x^\top P x$','Interpreter','latex');
title(ax3,'Lyapunov Function (log scale)');
legend(ax3,'Location','northeast');

%% Hill frame trajectory
ax4 = subplot(2,3,4);
hold(ax4,'on'); grid(ax4,'on'); axis(ax4,'equal');
plot(ax4, X_std(2,:), X_std(1,:), 'Color', c_std, 'LineWidth', 1.4, 'DisplayName', 'Std MPC');
plot(ax4, X_eco(2,:), X_eco(1,:), 'Color', c_eco, 'LineWidth', 1.4, 'DisplayName', 'EMPC');
plot(ax4, X0(2), X0(1), 'o', 'Color', [0.5 0.5 0.5], ...
     'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 6, 'HandleVisibility','off');
plot(ax4, 0, 0, '+k', 'MarkerSize', 8, 'LineWidth', 1.5, 'HandleVisibility','off');
xlabel(ax4,'Along-track $y$ [m]','Interpreter','latex');
ylabel(ax4,'Radial $x$ [m]','Interpreter','latex');
title(ax4,'Hill Frame Trajectory');
legend(ax4,'Location','northeast');

%% Control inputs — standard MPC
ax5 = subplot(2,3,5);
hold(ax5,'on'); grid(ax5,'on');
plot(ax5, t_min, U_std(1,:)*1000, 'Color', c_x, 'LineWidth', 1.2, 'DisplayName', '$u_x$');
plot(ax5, t_min, U_std(2,:)*1000, 'Color', c_y, 'LineWidth', 1.2, 'DisplayName', '$u_y$');
plot(ax5, t_min, U_std(3,:)*1000, 'Color', c_z, 'LineWidth', 1.2, 'DisplayName', '$u_z$');
yline(ax5,  con.u_max*1000, '--k', 'LineWidth', 0.8);
yline(ax5, -con.u_max*1000, '--k', 'LineWidth', 0.8);
xlabel(ax5,'Time [min]');
ylabel(ax5,'Accel. [mm/s$^2$]','Interpreter','latex');
title(ax5,'Control Inputs — Std MPC');
legend(ax5,'Location','northeast','Interpreter','latex','FontSize',6);

%% Control inputs — EMPC
ax6 = subplot(2,3,6);
hold(ax6,'on'); grid(ax6,'on');
plot(ax6, t_min, U_eco(1,:)*1000, 'Color', c_x, 'LineWidth', 1.2, 'DisplayName', '$u_x$');
plot(ax6, t_min, U_eco(2,:)*1000, 'Color', c_y, 'LineWidth', 1.2, 'DisplayName', '$u_y$');
plot(ax6, t_min, U_eco(3,:)*1000, 'Color', c_z, 'LineWidth', 1.2, 'DisplayName', '$u_z$');
yline(ax6,  con.u_max*1000, '--k', 'LineWidth', 0.8);
yline(ax6, -con.u_max*1000, '--k', 'LineWidth', 0.8);
xlabel(ax6,'Time [min]');
ylabel(ax6,'Accel. [mm/s$^2$]','Interpreter','latex');
title(ax6,'Control Inputs — EMPC');
legend(ax6,'Location','northeast','Interpreter','latex','FontSize',6);


%%  EXPORT JSON

metadata = struct(...
    'a_km',           a / 1000, ...
    'e',              e, ...
    'i_deg',          rad2deg(i), ...
    'period_min',     T / 60, ...
    'altitude_km',    (a - R_earth) / 1000, ...
    'dt',             dt, ...
    'n_steps',        n_steps, ...
    'controller',     'EMPC_regulation', ...
    'dynamics_model', 'CWH', ...
    'horizon_N',      N, ...
    'u_max',          con.u_max, ...
    'Q_stage',        'zeros', ...
    'R_eco',          'eye3', ...
    'dv_std_ms',      dv_std(end), ...
    'dv_eco_ms',      dv_eco(end));

data = struct(...
    'metadata', metadata, ...
    't',        t, ...
    'X_std',    X_std', ...
    'X_eco',    X_eco', ...
    'U_std',    U_std', ...
    'U_eco',    U_eco');

json_str = jsonencode(data, 'PrettyPrint', true);
fid = fopen('exports/scenarios/sim_empc_regulation.json', 'w');
fprintf(fid, '%s', json_str);
fclose(fid);

fprintf('JSON exported.\n');
