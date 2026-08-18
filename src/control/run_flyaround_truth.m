function res = run_flyaround_truth(run_cfg, matrices, params, plant)
% run_flyaround_truth.m
% Economic (pure-fuel) NMC fly-around run with a HIGH-FIDELITY TRUTH PLANT.
%
% Inputs
%   run_cfg  — as run_flyaround.m: .use_nmc_manifold .use_periodic_terminal
%              .use_band .x_star (6xP) .n_sim_override
%   matrices — CW condensed-QP matrices (.A_d .B_d .Sx .Su .Sx_pos_all
%              .Su_pos_all .H .A_u .b_u), built exactly as in S2
%   params   — .x0 .rho_min .rho_max .u_max .N .P .n .dt .n_sim [.verbose]
%   plant    — truth-plant config: .a .e .i [.RAAN .argp .nu] .p_t [.p_c]
%
% Output: res struct
%   .x_log (TRUE Hill states) .u_log .exitflag_log .range_log .dv_log .u_norm
%   .dv_total .dv_inject .dv_steady .track_err (vs Pi*) .n_sim .n_infeasible

%% Constants
constants;

%% Unpack matrices (controller CW model)
A_d        = matrices.A_d;   %#ok<NASGU>  (kept for interface parity)
Sx         = matrices.Sx;
Su         = matrices.Su;
Sx_pos_all = matrices.Sx_pos_all;
Su_pos_all = matrices.Su_pos_all;
H          = matrices.H;
A_u        = matrices.A_u;
b_u        = matrices.b_u;

%% Unpack params
x0      = params.x0;
rho_min = params.rho_min;
rho_max = params.rho_max;
N       = params.N;
P       = params.P;
n_orb   = params.n;
dt      = params.dt;
n_sim   = params.n_sim;
verbose = true; if isfield(params, 'verbose'), verbose = params.verbose; end

if isfield(run_cfg, 'n_sim_override') && ~isempty(run_cfg.n_sim_override)
    n_sim = min(n_sim, run_cfg.n_sim_override);
end

use_nmc      = run_cfg.use_nmc_manifold;
use_periodic = run_cfg.use_periodic_terminal;
use_band     = ~isfield(run_cfg, 'use_band') || run_cfg.use_band;
x_star = [];
if isfield(run_cfg, 'x_star') && ~isempty(run_cfg.x_star)
    x_star = run_cfg.x_star;
end

n_x = 6;

%% NMC manifold equality matrices (as run_flyaround.m)
Sx_term     = Sx((N-1)*n_x+1 : N*n_x, :);
Su_term     = Su((N-1)*n_x+1 : N*n_x, :);
Aeq_nmc     = [Su_term(4, :);  Su_term(5, :) + 2*n_orb * Su_term(1, :)];
Sx_term_nmc = [Sx_term(4, :);  Sx_term(5, :) + 2*n_orb * Sx_term(1, :)];

%% Truth-plant ICs
if ~isfield(plant, 'RAAN'), plant.RAAN = 0; end
if ~isfield(plant, 'argp'), plant.argp = 0; end
if ~isfield(plant, 'nu'),   plant.nu   = 0; end
p_t = plant.p_t;  if ~isfield(p_t, 'mu'), p_t.mu = mu; end
if isfield(plant, 'p_c'), p_c = plant.p_c; else, p_c = p_t; end
if ~isfield(p_c, 'mu'), p_c.mu = mu; end

X_t = oe2eci(plant.a, plant.e, plant.i, plant.RAAN, plant.argp, plant.nu, mu);
X_c = hill2eci(X_t, x0);

%% Allocate
x_log        = zeros(n_x, n_sim + 1);
u_log        = zeros(3,   n_sim);
exitflag_log = zeros(1,   n_sim);
track_err    = zeros(1,   n_sim);

x_log(:, 1) = x0;
U_warm      = zeros(3*N, 1);
qp_opts     = optimoptions('quadprog', 'Display', 'off');
opts_ode    = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);

%% Receding-horizon loop (controller CW, plant truth)
for k = 1:n_sim
    xk = x_log(:, k);           % TRUE measured Hill state

    % 1. SCA warm-start prediction (CW)
    x_pred = Sx * xk + Su * U_warm;
    r_nom  = reshape(x_pred, n_x, N);
    r_nom  = r_nom(1:3, :);

    % 2. Band inequality constraints (SCA-linearised)
    if use_band
        r_norms = vecnorm(r_nom, 2, 1);
        A_band  = zeros(N, 3*N);
        b_outer = zeros(N, 1);
        b_inner = zeros(N, 1);
        for ki = 1:N
            idx3 = (ki-1)*3+1 : ki*3;
            if r_norms(ki) < 1e-6, nk = [0;1;0]; else, nk = r_nom(:, ki)/r_norms(ki); end
            row           = nk' * Su_pos_all(idx3, :);
            b_sx          = nk' * Sx_pos_all(idx3, :) * xk;
            A_band(ki, :) = row;
            b_outer(ki)   = rho_max - b_sx;
            b_inner(ki)   = rho_min - b_sx;
        end
        A_ineq = [A_u;  A_band; -A_band];
        b_ineq = [b_u; b_outer; -b_inner];
    else
        A_ineq = A_u;  b_ineq = b_u;
    end

    % 3. Terminal equality (periodic Pi* pin, or NMC manifold, or none)
    if use_nmc && use_periodic && ~isempty(x_star)
        x_ref_col = x_star(:, mod(k + N - 1, P) + 1);
        Aeq = Su_term;
        beq = x_ref_col - Sx_term * xk;
    elseif use_nmc
        Aeq = Aeq_nmc;
        beq = -Sx_term_nmc * xk;
    else
        Aeq = [];  beq = [];
    end

    % 4. Solve QP (pure fuel: f = 0)
    [U_opt, ~, exitflag] = quadprog(H, zeros(3*N,1), A_ineq, b_ineq, Aeq, beq, [], [], U_warm, qp_opts);
    exitflag_log(k) = exitflag;

    if exitflag > 0
        u_opt  = U_opt(1:3);
        U_warm = [U_opt(4:end); zeros(3,1)];
    else
        u_opt  = zeros(3,1);
        U_warm = [U_warm(4:end); zeros(3,1)];
    end
    u_log(:, k) = u_opt;

    % Track error vs the CW periodic orbit at this phase (position)
    if ~isempty(x_star)
        xr = x_star(:, mod(k-1, P) + 1);
        track_err(k) = norm(xk(1:3) - xr(1:3));
    end

    % 5. Propagate TRUTH plant over dt with the ZOH Hill-frame control
    tk  = (k-1) * dt;
    rhs = @(t, Y) fly_step_rhs(t, Y, p_t, p_c, u_opt);
    [~, Ysol] = ode113(rhs, [tk, tk + dt], [X_t; X_c], opts_ode);
    X_t = Ysol(end, 1:6).';
    X_c = Ysol(end, 7:12).';

    x_log(:, k+1) = eci2hill(X_t, X_c);
end

%% Post-process
range_log = vecnorm(x_log(1:3, :));
dv_log    = vecnorm(u_log, 2, 1) * dt;
u_norm    = vecnorm(u_log, 2, 1);

rolling_mean = movmean(u_norm, max(1, round(P/4)));
inject_end   = find(rolling_mean < 1e-4, 1);
if isempty(inject_end), inject_end = min(P, n_sim); end   % fallback: 1 orbit

dv_total  = sum(dv_log);
dv_inject = sum(dv_log(1:inject_end));
dv_steady = sum(dv_log(inject_end+1:end));

n_infeasible = sum(exitflag_log <= 0);

if verbose
    fprintf('  [fly-truth] e=%.3f: dv_tot=%.4f (inj=%.4f, steady=%.4f), maxerr=%.1f m, infeas=%d\n', ...
        plant.e, dv_total, dv_inject, dv_steady, max(track_err), n_infeasible);
end

res.x_log        = x_log;
res.u_log        = u_log;
res.exitflag_log = exitflag_log;
res.range_log    = range_log;
res.dv_log       = dv_log;
res.u_norm       = u_norm;
res.dv_total     = dv_total;
res.dv_inject    = dv_inject;
res.dv_steady    = dv_steady;
res.track_err    = track_err;
res.n_sim        = n_sim;
res.n_infeasible = n_infeasible;

end

% local
function dY = fly_step_rhs(t, Y, p_t, p_c, u_hill)
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