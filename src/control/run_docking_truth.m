function res = run_docking_truth(cost_struct, con, params, plant)
% run_docking_truth.m
% Closed-loop docking run with a HIGH-FIDELITY TRUTH PLANT.
%
% Inputs (cost_struct, con, params) match run_docking.m. Additional:
%   plant - truth-plant configuration:
%             .a, .e, .i      target orbit (sma [m], ecc, incl [rad])
%             .RAAN,.argp,.nu target orientation/phase [rad] (default 0)
%             .p_t, .p_c      truth_propagator_at param structs for target and
%                             chaser (perturbation flags, ballistic coeffs).
%                             .mu is filled from constants if absent.
%
% Output: res struct
%   .X_hill    [6 x n_run]  TRUE relative state (Hill frame) at each step
%   .U         [3 x n_run]  applied control (Hill accel) per step
%   .ef        [1 x n_run]  true position-error norm [m]
%   .dv                     total Delta-v [m/s]
%   .term_dist              true terminal position distance to x_dock [m]
%   .conv_step              step index at convergence (NaN if none)
%   .t_vec     [1 x n_run]  time [s]
%   .n_ctrl                 n_run - 1
%   .los_viol_frac          fraction of steps violating the LoS cone (if active)
%   .kos_breach             true if TRUE trajectory entered r_KOS before conv.

%% Constants
constants;

%% Unpack params
x0      = params.x0;
x_dock  = params.x_dock;
dt      = params.dt;
N       = params.N;
u_max   = params.u_max;
t_final = params.t_final;

conv_tol = 1.0; if isfield(params, 'conv_tol'), conv_tol = params.conv_tol; end
r_KOS    = 5.0; if isfield(params, 'r_KOS'),    r_KOS    = params.r_KOS;    end
verbose  = true; if isfield(params, 'verbose'), verbose = params.verbose;   end

n_steps = round(t_final / dt);
if isfield(params, 'n_steps_override') && ~isempty(params.n_steps_override)
    n_steps = min(n_steps, params.n_steps_override);
end

%% Controller's internal CW model 
n_ctrl_mm = sqrt(mu / plant.a^3);
Ac = [0           0     0    1          0    0  ;
      0           0     0    0          1    0  ;
      0           0     0    0          0    1  ;
      3*n_ctrl_mm^2 0    0    0    2*n_ctrl_mm  0  ;
      0           0     0   -2*n_ctrl_mm 0    0  ;
      0           0  -n_ctrl_mm^2 0      0    0 ];
Bc = [zeros(3,3); eye(3)];
[Ad, Bd] = discretize(Ac, Bc, dt);

%% Constraint toggles consumed by constraints.m
con.u_max = u_max;
if ~isfield(con, 'y_min_active'),        con.y_min_active    = false; end
if ~isfield(con, 'los_cone'),            con.los_cone.active = false;
elseif ~isfield(con.los_cone, 'active'), con.los_cone.active = false; end
if con.los_cone.active
    if ~isfield(con.los_cone, 'half_angle'), con.los_cone.half_angle = pi/6; end
    if ~isfield(con.los_cone, 'n_faces'),    con.los_cone.n_faces    = 10;   end
end
if isfield(con, 'passive_safety'), con.passive_safety.active = false; end

%% Truth-plant setup: absolute ECI ICs
if ~isfield(plant, 'RAAN'), plant.RAAN = 0; end
if ~isfield(plant, 'argp'), plant.argp = 0; end
if ~isfield(plant, 'nu'),   plant.nu   = 0; end
p_t = plant.p_t;  if ~isfield(p_t, 'mu'), p_t.mu = mu; end
p_c = plant.p_c;  if ~isfield(p_c, 'mu'), p_c.mu = mu; end

X_t = oe2eci(plant.a, plant.e, plant.i, plant.RAAN, plant.argp, plant.nu, mu);
X_c = hill2eci(X_t, x0);                 % chaser ECI seeded from Hill IC

%% Allocate
U_warm = zeros(3*N, 1);
X_hill = zeros(6, n_steps + 1);
U      = zeros(3, n_steps + 1);
ef     = zeros(1, n_steps + 1);
los_viol = false(1, n_steps + 1);

X_hill(:,1) = x0;
e_k         = x0 - x_dock;
conv_step   = NaN;
kos_breach  = false;
n_run       = 0;
opts_ode    = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);

tan_half = tan(con.los_cone.half_angle);   % for LoS check on true state

for k = 1:n_steps
    ef(k) = norm(e_k(1:3));

    % LoS / KOS bookkeeping on the TRUE relative position (before conv check)
    p_rel = e_k(1:3);                      % relative to x_dock
    if con.los_cone.active
        % cone axis +y (along-track): require y>=0 and sqrt(x^2+z^2) <= y*tan
        los_viol(k) = (p_rel(2) < 0) || ...
                      (hypot(p_rel(1), p_rel(3)) > p_rel(2)*tan_half + 1e-9);
    end
    if norm(X_hill(1:3,k)) < r_KOS && isnan(conv_step)
        kos_breach = true;                 % entered keep-out sphere pre-convergence
    end

    if ef(k) < conv_tol
        conv_step = k;  n_run = k;  break;
    end

    %% MPC solve (CW internal model, true measured error state)
    [u_opt, U_opt, ~, ~, ~] = mpc_regulation( ...
        e_k, Ad, Bd, cost_struct.Q, cost_struct.R, cost_struct.P, N, con);
    U_warm = [U_opt(4:end); U_opt(end-2:end)];   %#ok<NASGU>
    U(:,k) = u_opt;

    %% Propagate TRUTH plant over dt with the ZOH Hill-frame control
    tk = (k-1) * dt;
    rhs = @(t, Y) plant_step_rhs(t, Y, p_t, p_c, u_opt, mu);
    [~, Ysol] = ode113(rhs, [tk, tk + dt], [X_t; X_c], opts_ode);
    X_t = Ysol(end, 1:6).';
    X_c = Ysol(end, 7:12).';

    %% Recover true relative state
    e_k = eci2hill(X_t, X_c) - x_dock;
    X_hill(:,k+1) = eci2hill(X_t, X_c);
    n_run = k + 1;
end

%% Trim
X_hill = X_hill(:, 1:n_run);
U      = U(:, 1:n_run);
ef     = ef(1:n_run);
los_viol = los_viol(1:n_run);
n_ctrl = n_run - 1;

%% Metrics
dv        = sum(vecnorm(U(:, 1:n_ctrl), 2, 1)) * dt;
term_dist = norm(X_hill(1:3, end) - x_dock(1:3));
t_vec     = (0:n_run-1) * dt;
los_frac  = mean(los_viol(1:max(n_ctrl,1)));

if verbose
    if isnan(conv_step)
        fprintf('  [truth] e=%.3f: NO CONVERGE in %d steps, term=%.2f m, dv=%.4f\n', ...
                plant.e, n_steps, term_dist, dv);
    else
        fprintf('  [truth] e=%.3f: conv @%d (%.1f min), term=%.3f m, dv=%.4f\n', ...
                plant.e, conv_step, t_vec(conv_step)/60, term_dist, dv);
    end
end

res.X_hill        = X_hill;
res.U             = U;
res.ef            = ef;
res.dv            = dv;
res.term_dist     = term_dist;
res.conv_step     = conv_step;
res.t_vec         = t_vec;
res.n_ctrl        = n_ctrl;
res.los_viol_frac = los_frac;
res.kos_breach    = kos_breach;

end

%local
function dY = plant_step_rhs(t, Y, p_t, p_c, u_hill, mu)
% Coupled target+chaser truth dynamics with a constant Hill-frame control
% acceleration applied to the chaser (rotated to ECI via the target basis).
X_t = Y(1:6);  X_c = Y(7:12);

dX_t = truth_propagator_at(t, X_t, p_t);
dX_c = truth_propagator_at(t, X_c, p_c);

% Instantaneous target LVLH basis (ECI -> Hill); Hill -> ECI is R.'
r_t = X_t(1:3);  v_t = X_t(4:6);
x_hat = r_t / norm(r_t);
h_vec = cross(r_t, v_t);
z_hat = h_vec / norm(h_vec);
y_hat = cross(z_hat, x_hat);
R = [x_hat.'; y_hat.'; z_hat.'];

dX_c(4:6) = dX_c(4:6) + R.' * u_hill;    % apply thrust in ECI

dY = [dX_t; dX_c];
end