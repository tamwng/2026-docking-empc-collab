function res = run_flyaround(run_cfg, matrices, params)
% run_flyaround.m — Single EMPC fly-around run (CWH, SCA band constraints).
%
% Inputs
%   run_cfg   — .use_nmc_manifold      (bool) add NMC manifold Aeq/beq
%               .use_periodic_terminal (bool) pin terminal state to Π*(k mod P)
%               .x_star                (6×P)  reference orbit, [] if unused
%               .n_sim_override        (scalar or []) cap simulation length
%   matrices  — .A_d .B_d .Sx .Su .Sx_pos_all .Su_pos_all .H .f_vec .A_u .b_u
%   params    — .x0 .rho_min .rho_max .u_max .N .P .n .dt .n_sim
%
% Output  res struct (all fields named per the caller spec):
%   .x_log .u_log .exitflag_log .range_log .dv_log .u_norm
%   .inject_end .b_est .ratio_2to1 .dv_total .dv_inject .dv_steady .n_sim

%% Unpack matrices
A_d        = matrices.A_d;
B_d        = matrices.B_d;
Sx         = matrices.Sx;
Su         = matrices.Su;
Sx_pos_all = matrices.Sx_pos_all;
Su_pos_all = matrices.Su_pos_all;
H          = matrices.H;
A_u        = matrices.A_u;
b_u        = matrices.b_u;

% F_mat: f_vec at step k = F_mat * x(k).  Zeros for fuel-only cost.
if isfield(matrices, 'F_mat')
    F_mat = matrices.F_mat;
else
    F_mat = zeros(size(H, 1), 6);
end

%% Unpack params
x0      = params.x0;
rho_min = params.rho_min;
rho_max = params.rho_max;
N       = params.N;
P       = params.P;
n_orb   = params.n;     % mean motion [rad/s]
dt      = params.dt;
n_sim   = params.n_sim;

if isfield(run_cfg, 'n_sim_override') && ~isempty(run_cfg.n_sim_override)
    n_sim = min(n_sim, run_cfg.n_sim_override);
end

use_nmc      = run_cfg.use_nmc_manifold;
use_periodic = run_cfg.use_periodic_terminal;
use_band     = ~isfield(run_cfg, 'use_band') || run_cfg.use_band;
x_star       = [];
if isfield(run_cfg, 'x_star') && ~isempty(run_cfg.x_star)
    x_star = run_cfg.x_star;   % 6 × P  (used by terminal equality and stage cost routes)
end

% f_ref_mat: reference-dependent linear term in f_vec.
% f_ref_k = f_ref_mat * X_ref_k  where X_ref_k is the stacked 6N horizon reference.
% Non-empty only for the periodic stage cost route (Run 5).
if isfield(matrices, 'f_ref_mat') && ~isempty(x_star)
    f_ref_mat    = matrices.f_ref_mat;   % 3N × 6N
    has_ref_cost = true;
else
    f_ref_mat    = [];
    has_ref_cost = false;
end

n_x = 6;

%% NMC manifold equality matrices — constant, computed once from Sx/Su
% Encodes: vx(N)=0  and  vy(N)+2*n*x(N)=0
Sx_term     = Sx((N-1)*n_x+1 : N*n_x, :);   % 6 × 6
Su_term     = Su((N-1)*n_x+1 : N*n_x, :);   % 6 × 3N
Aeq_nmc     = [Su_term(4, :);
               Su_term(5, :) + 2*n_orb * Su_term(1, :)];     % 2 × 3N
Sx_term_nmc = [Sx_term(4, :);
               Sx_term(5, :) + 2*n_orb * Sx_term(1, :)];     % 2 × 6

%% Allocate
x_log        = zeros(n_x, n_sim + 1);
u_log        = zeros(3,   n_sim);
exitflag_log = zeros(1,   n_sim);

x_log(:, 1) = x0;
U_warm      = zeros(3*N, 1);
opts        = optimoptions('quadprog', 'Display', 'off');

fprintf('  [run_flyaround] %d steps  nmc_manifold=%d  periodic_terminal=%d\n', ...
    n_sim, use_nmc, use_periodic);
t_start = tic;

%% Receding-horizon loop
for k = 1:n_sim

    % 1. Warm-start predicted positions (SCA linearisation point)
    x_pred = Sx * x_log(:, k) + Su * U_warm;   % 6N × 1
    r_nom  = reshape(x_pred, n_x, N);
    r_nom  = r_nom(1:3, :);                     % 3 × N

    % 2. Build inequality constraints
    if use_band
        r_norms = vecnorm(r_nom, 2, 1);
        A_band  = zeros(N, 3*N);
        b_outer = zeros(N, 1);
        b_inner = zeros(N, 1);
        for ki = 1:N
            idx3 = (ki-1)*3+1 : ki*3;
            if r_norms(ki) < 1e-6
                nk = [0; 1; 0];
            else
                nk = r_nom(:, ki) / r_norms(ki);
            end
            row           = nk' * Su_pos_all(idx3, :);
            b_sx          = nk' * Sx_pos_all(idx3, :) * x_log(:, k);
            A_band(ki, :) = row;
            b_outer(ki)   = rho_max - b_sx;
            b_inner(ki)   = rho_min - b_sx;
        end
        A_ineq = [A_u;  A_band; -A_band];
        b_ineq = [b_u; b_outer; -b_inner];
    else
        A_ineq = A_u;
        b_ineq = b_u;
    end

    % 3. Equality constraints depending on run configuration
    if use_nmc && use_periodic && ~isempty(x_star)
        % Run 2: pin terminal state to Π*(k+N-1 mod P) — Keerthi–Gilbert /
        % Theorem 2.24 (periodic terminal equality → closed-loop stability).
        % When N=P the phase mod(k+N-1,P) = mod(k-1,P), closing exactly one period.
        % The 6-row equality subsumes the 2-row NMC manifold constraint.
        x_ref_col = x_star(:, mod(k + N - 1, P) + 1);
        Aeq = Su_term;
        beq = x_ref_col - Sx_term * x_log(:, k);
    elseif use_nmc
        % Run 2: NMC manifold equality only — forces endpoint onto the manifold
        % without specifying which orbit
        Aeq = Aeq_nmc;
        beq = -Sx_term_nmc * x_log(:, k);
    else
        % Run 1: no terminal constraint — degenerate case, U*=0 is globally optimal
        Aeq = [];
        beq = [];
    end

    % 4. Solve QP
    f_vec_k = F_mat * x_log(:, k);
    if has_ref_cost
        cols    = mod(k - 1 + (1:N), P) + 1;
        X_ref_k = reshape(x_star(:, cols), [], 1);   % 6N × 1
        f_vec_k = f_vec_k + f_ref_mat * X_ref_k;
    end
    [U_opt, ~, exitflag] = quadprog(H, f_vec_k, A_ineq, b_ineq, Aeq, beq, [], [], U_warm, opts);
    exitflag_log(k) = exitflag;

    % 5. Apply first control; warm-start shift
    if exitflag > 0
        u_log(:, k) = U_opt(1:3);
        U_warm = [U_opt(4:end); zeros(3, 1)];
    else
        u_log(:, k) = zeros(3, 1);
        U_warm = [U_warm(4:end); zeros(3, 1)];
    end

    % 6. Propagate
    x_log(:, k+1) = A_d * x_log(:, k) + B_d * u_log(:, k);
end

fprintf('  [run_flyaround] done (%.1f s wall-clock)\n', toc(t_start));

%% Post-process
range_log    = vecnorm(x_log(1:3, :));            % 1 × (n_sim+1)
dv_log       = vecnorm(u_log, 2, 1) * dt;         % per-step ΔV [m/s]
u_norm       = vecnorm(u_log, 2, 1);              % 1 × n_sim

rolling_mean = movmean(u_norm, max(1, round(P/4)));
inject_end   = find(rolling_mean < 1e-4, 1);
if isempty(inject_end); inject_end = n_sim; end

dv_total  = sum(dv_log);
dv_inject = sum(dv_log(1:inject_end));
dv_steady = sum(dv_log(inject_end+1:end));

range_ss_start = max(1, n_sim - 2*P + 1);
range_ss       = range_log(range_ss_start : end);
b_est          = min(range_ss);
ratio_2to1     = max(range_ss) / max(b_est, 1e-6);

last_orbit_start = max(1, n_sim - P);
u_mean_last = mean(u_norm(last_orbit_start : end));
u_max_last  =  max(u_norm(last_orbit_start : end));

%% Summary (mirrors the existing sim_scenario2_flyaround summary block)
fprintf('Total ΔV:       %.4f m/s\n', dv_total);
fprintf('  Injection ΔV: %.4f m/s  (steps 1–%d)\n', dv_inject, inject_end);
fprintf('  Steady-state: %.4f m/s  (steps %d–%d)\n', dv_steady, inject_end+1, n_sim);
fprintf('Discovered NMT semi-axis b ≈ %.1f m\n', b_est);
fprintf('Range ratio max/min (expect ≈ 2.0): %.3f\n', ratio_2to1);
if use_band
    fprintf('Band violations: %d\n', sum(range_log < rho_min | range_log > rho_max));
end
fprintf('Infeasible QP steps: %d\n', sum(exitflag_log <= 0));
fprintf('Mean ‖u‖ last orbit: %.2e m/s²\n', u_mean_last);
fprintf('Max  ‖u‖ last orbit: %.2e m/s²\n', u_max_last);

%% Pack result
res.x_log        = x_log;
res.u_log        = u_log;
res.exitflag_log = exitflag_log;
res.range_log    = range_log;
res.dv_log       = dv_log;
res.u_norm       = u_norm;
res.inject_end   = inject_end;
res.b_est        = b_est;
res.ratio_2to1   = ratio_2to1;
res.dv_total     = dv_total;
res.dv_inject    = dv_inject;
res.dv_steady    = dv_steady;
res.n_sim        = n_sim;

end
