function res = run_closing_truth(cost_struct, con, params, plant)
% run_closing_truth.m
% Closed-loop FAR-RANGE CLOSING run with a high-fidelity truth plant.
%
% Inputs match run_closing.m (cost_struct, con, params with .x_switch) plus
% the same `plant` struct as run_docking_truth.m (target orbit + p_t/p_c).
%
% Output: res struct
%   .X_hill .U .ef .dv .term_dist .conv_step .t_vec .n_ctrl .reached_ball

%% Constants
constants;

%% Unpack params
x0       = params.x0;
x_switch = params.x_switch;
dt       = params.dt;
N        = params.N;
u_max    = params.u_max;
t_final  = params.t_final;

conv_tol = 50;  if isfield(params, 'conv_tol'), conv_tol = params.conv_tol; end
verbose  = true; if isfield(params, 'verbose'), verbose = params.verbose;   end

n_steps = round(t_final / dt);
if isfield(params, 'n_steps_override') && ~isempty(params.n_steps_override)
    n_steps = min(n_steps, params.n_steps_override);
end

%% Controller CW model (mean motion at plant sma -> circular assumption)
n_mm = sqrt(mu / plant.a^3);
Ac = [0        0     0    1        0    0  ;
      0        0     0    0        1    0  ;
      0        0     0    0        0    1  ;
      3*n_mm^2 0     0    0    2*n_mm   0  ;
      0        0     0   -2*n_mm   0    0  ;
      0        0  -n_mm^2 0        0    0 ];
Bc = [zeros(3,3); eye(3)];
[Ad, Bd] = discretize(Ac, Bc, dt);

%% Constraint toggles
con.u_max = u_max;
if ~isfield(con, 'y_min_active'),        con.y_min_active    = false; end
if ~isfield(con, 'los_cone'),            con.los_cone.active = false;
elseif ~isfield(con.los_cone, 'active'), con.los_cone.active = false; end
if isfield(con, 'use_terminal_ball') && con.use_terminal_ball
    con.terminal_ball.active = true;
elseif isfield(con, 'terminal_ball')
    con.terminal_ball.active = false;
end
if isfield(con, 'passive_safety'), con.passive_safety.active = false; end

%% Truth-plant ICs
if ~isfield(plant, 'RAAN'), plant.RAAN = 0; end
if ~isfield(plant, 'argp'), plant.argp = 0; end
if ~isfield(plant, 'nu'),   plant.nu   = 0; end
p_t = plant.p_t;  if ~isfield(p_t, 'mu'), p_t.mu = mu; end
p_c = plant.p_c;  if ~isfield(p_c, 'mu'), p_c.mu = mu; end

X_t = oe2eci(plant.a, plant.e, plant.i, plant.RAAN, plant.argp, plant.nu, mu);
X_c = hill2eci(X_t, x0);

%% Allocate
X_hill = zeros(6, n_steps + 1);
U      = zeros(3, n_steps + 1);
ef     = zeros(1, n_steps + 1);

X_hill(:,1)  = x0;
e_k          = x0 - x_switch;
conv_step    = NaN;
n_run        = 0;
opts_ode     = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);

for k = 1:n_steps
    ef(k) = norm(e_k);                 % full-state handover norm

    if ef(k) < conv_tol
        conv_step = k;  n_run = k;  break;
    end

    [u_opt, ~, ~, ~, ~] = mpc_regulation( ...
        e_k, Ad, Bd, cost_struct.Q, cost_struct.R, cost_struct.P, N, con);
    U(:,k) = u_opt;

    tk  = (k-1) * dt;
    rhs = @(t, Y) closing_step_rhs(t, Y, p_t, p_c, u_opt);
    [~, Ysol] = ode113(rhs, [tk, tk + dt], [X_t; X_c], opts_ode);
    X_t = Ysol(end, 1:6).';
    X_c = Ysol(end, 7:12).';

    e_k          = eci2hill(X_t, X_c) - x_switch;
    X_hill(:,k+1) = eci2hill(X_t, X_c);
    n_run        = k + 1;
end

%% Trim + metrics
X_hill = X_hill(:, 1:n_run);
U      = U(:, 1:n_run);
ef     = ef(1:n_run);
n_ctrl = n_run - 1;

dv        = sum(vecnorm(U(:, 1:n_ctrl), 2, 1)) * dt;
term_dist = norm(X_hill(1:3, end) - x_switch(1:3));
t_vec     = (0:n_run-1) * dt;
reached   = ~isnan(conv_step);

if verbose
    if reached
        fprintf('  [closing] e=%.3f: handover @%d (%.2f h), term=%.1f m, dv=%.4f\n', ...
                plant.e, conv_step, t_vec(conv_step)/3600, term_dist, dv);
    else
        fprintf('  [closing] e=%.3f: NO handover in %d steps, term=%.1f m, dv=%.4f\n', ...
                plant.e, n_steps, term_dist, dv);
    end
end

res.X_hill    = X_hill;
res.U         = U;
res.ef        = ef;
res.dv        = dv;
res.term_dist = term_dist;
res.conv_step = conv_step;
res.t_vec     = t_vec;
res.n_ctrl    = n_ctrl;
res.reached_ball = reached;

end

% local
function dY = closing_step_rhs(t, Y, p_t, p_c, u_hill)
X_t = Y(1:6);  X_c = Y(7:12);
dX_t = truth_propagator_at(t, X_t, p_t);
dX_c = truth_propagator_at(t, X_c, p_c);

r_t = X_t(1:3);  v_t = X_t(4:6);
x_hat = r_t / norm(r_t);
h_vec = cross(r_t, v_t);
z_hat = h_vec / norm(h_vec);
y_hat = cross(z_hat, x_hat);
R = [x_hat.'; y_hat.'; z_hat.'];

dX_c(4:6) = dX_c(4:6) + R.' * u_hill;
dY = [dX_t; dX_c];
end