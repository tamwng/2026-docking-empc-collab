function res = run_docking(cost_struct, con, params)
% run_docking.m
% Single MPC/EMPC docking-phase run (CWH dynamics, condensed QP via mpc_regulation).
%
% Drives e_k = x_k − x_dock → 0 (x_dock is the docking port, typically origin).
% Convergence criterion: norm(e_k(1:3)) < conv_tol (position error [m], default 1 m)
%
% Inputs
%   cost_struct — .Q [6×6]  .R [3×3]  .P [6×6]
%                 When con.terminal_eq=true, P is overridden to zero internally
%                 by mpc_regulation (terminal cost replaced by equality).
%   con         — constraint config struct:
%                   .u_max              [required] symmetric thrust bound [m/s²]
%                   .y_min_active       [bool, default false]
%                   .los_cone.active    [bool, default false]
%                   .terminal_eq        [bool] hard x(N)=x_dock equality
%                   .use_passive_safety [bool]
%                   .passive_safety.*   (lazily populated if absent)
%   params      — run parameters:
%                   .x0       [6×1] initial state (Hill frame)
%                   .x_dock   [6×1] docking port target (typically zeros(6,1))
%                   .dt       [s] sampling time
%                   .N        prediction horizon
%                   .u_max    [m/s²] thrust bound
%                   .t_final  [s] max simulation time
%                   .conv_tol (optional, default 1.0 m)
%                   .r_KOS    (optional, default 5 m)
%                   .N_safe   (optional, default 554 steps ≈ 1 orbital period)
%                   .n_steps_override (optional, caps run length)
%                   .verbose  (optional, default true)
%
% Output: res struct
%   .X .U .ef .flags .is_safe  — [*×n_run] trimmed arrays
%   .dv .term_dist .ps_frac .conv_step    — scalars
%   .u_norm_seq .all_X_pred               — [1×n_ctrl] / cell
%   .turnpike_steps .turnpike_fraction    — vector / scalar
%   .tp_label                             — string
%   .t_vec                                — [1×n_run] time in seconds
%   .n_ctrl                               — n_run - 1
%   .n_ps_h                               — [1×n_ctrl] PS rows added per step

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
N_safe   = 554; if isfield(params, 'N_safe'),   N_safe   = params.N_safe;   end

n_steps = round(t_final / dt);
if isfield(params, 'n_steps_override') && ~isempty(params.n_steps_override)
    n_steps = min(n_steps, params.n_steps_override);
end

%% CWH state-space 
Ac = [0      0     0    1     0    0  ;
      0      0     0    0     1    0  ;
      0      0     0    0     0    1  ;
      3*n^2  0     0    0     2*n  0  ;
      0      0     0   -2*n   0    0  ;
      0      0    -n^2  0     0    0 ];
Bc = [zeros(3,3); eye(3)];
[Ad, Bd] = discretize(Ac, Bc, dt);

% Precompute Ad powers for post-hoc passive safety check
Ad_pow_safe = zeros(6, 6, N_safe);
Ad_cur = Ad;
for j = 1:N_safe
    Ad_pow_safe(:,:,j) = Ad_cur;
    Ad_cur = Ad_cur * Ad;
end

%% Resolve constraint toggles
use_passive_safety = isfield(con, 'use_passive_safety') && con.use_passive_safety;
use_terminal_eq    = isfield(con, 'terminal_eq')        && con.terminal_eq;

%% Set .active flags consumed by constraints.m
con.u_max = u_max;
if ~isfield(con, 'y_min_active'),          con.y_min_active    = false; end
if ~isfield(con, 'los_cone'),              con.los_cone.active = false;
elseif ~isfield(con.los_cone, 'active'),   con.los_cone.active = false; end
if con.los_cone.active
    if ~isfield(con.los_cone, 'half_angle'), con.los_cone.half_angle = pi/6; end  % 30 deg default
    if ~isfield(con.los_cone, 'n_faces'),    con.los_cone.n_faces    = 10;   end
end

if use_passive_safety
    con.passive_safety.active = true;
    if ~isfield(con.passive_safety, 'Phi_array')
        fprintf('  Precomputing CW STMs (%d matrices)... ', N_safe);
        Phi_ps = zeros(6, 6, N_safe);
        for i = 1:N_safe
            Phi_ps(:,:,i) = cw_stm(n, i * dt);
        end
        fprintf('done.\n');
        con.passive_safety.Phi_array = Phi_ps;
    end
    if ~isfield(con.passive_safety, 'C_pos'),         con.passive_safety.C_pos         = [eye(3), zeros(3,3)]; end
    if ~isfield(con.passive_safety, 'r_KOS'),         con.passive_safety.r_KOS         = r_KOS;                end
    if ~isfield(con.passive_safety, 'safety_margin'), con.passive_safety.safety_margin  = 0;                   end
    if ~isfield(con.passive_safety, 'x_switch'),      con.passive_safety.x_switch       = x_dock;              end
elseif isfield(con, 'passive_safety')
    con.passive_safety.active = false;
end

%% Allocate trajectories
U_warm     = zeros(3*N, 1);
X          = zeros(6, n_steps + 1);
U          = zeros(3, n_steps + 1);
ef         = zeros(1, n_steps + 1);
flags_arr  = zeros(1, n_steps + 1);
n_ps_h     = zeros(1, n_steps);
u_norm_seq = zeros(1, n_steps);
all_X_pred = cell(1, n_steps);

X(:,1)    = x0;
e_k       = x0 - x_dock;
conv_step = NaN;
n_run     = 0;

%% Verbose flag — suppress for sweep runs
verbose = true;
if isfield(params, 'verbose'), verbose = params.verbose; end

%% MPC / EMPC loop
if verbose
    fprintf('\n  %5s  %10s  %12s  %8s  %4s  %6s  %8s\n', ...
        'step', '||u||', '||e||', 'term_eq', 'n_ps', 'flag', 'time[s]');
end

for k = 1:n_steps
    ef(k) = norm(e_k(1:3));   % position error only [m]

    if ef(k) < conv_tol
        conv_step = k;
        n_run     = k;
        break;
    end

    if use_passive_safety
        con.passive_safety.U_warm  = U_warm;
        con.passive_safety.x_switch = x_dock;   % actual-coord offset
    end

    t0 = tic;
    [u_opt, U_opt, X_pred_k, exitflag, n_ps_step] = mpc_regulation( ...
        e_k, Ad, Bd, cost_struct.Q, cost_struct.R, cost_struct.P, N, con);
    solve_t       = toc(t0);
    flags_arr(k)  = exitflag;
    n_ps_h(k)     = n_ps_step;
    u_norm_seq(k) = norm(u_opt);
    all_X_pred{k} = X_pred_k + x_dock;   % error → Hill frame

    if verbose
        fprintf('  %5d  %10.3e  %12.3e  %8s  %4d  %6d  %8.4f\n', ...
            k, u_norm_seq(k), ef(k), mat2str(use_terminal_eq), ...
            n_ps_step, exitflag, solve_t);
    end

    U_warm   = [U_opt(4:end); U_opt(end-2:end)];
    U(:,k)   = u_opt;
    X(:,k+1) = Ad * X(:,k) + Bd * u_opt;
    e_k      = X(:,k+1) - x_dock;
    n_run    = k + 1;
end

%% Trim to actual run length
X          = X(:,        1:n_run);
U          = U(:,        1:n_run);
ef         = ef(         1:n_run);
flags      = flags_arr(  1:n_run);
n_ctrl     = n_run - 1;
n_ps_h     = n_ps_h(     1:n_ctrl);
u_norm_seq = u_norm_seq( 1:n_ctrl);
all_X_pred = all_X_pred( 1:n_ctrl);

%% Post-hoc passive safety check
% Propagates each logged state forward N_safe free-drift steps; flags safe
% iff the drift trajectory never enters the r_KOS sphere around x_dock.
is_safe   = check_passive_safety(X - x_dock, n_run, Ad_pow_safe, r_KOS, N_safe);
ps_frac   = mean(is_safe);
dv        = sum(vecnorm(U(:, 1:n_ctrl), 2, 1)) * dt;
term_dist = norm(X(:,end) - x_dock);
t_vec     = (0:n_run-1) * dt;

%% Console summary
if verbose
    fprintf('\n── Summary ──────────────────────────────────────────────────\n');
    if ~isnan(conv_step)
        fprintf('  Converged at step %d  (%.1f min)\n', conv_step, t_vec(conv_step)/60);
    else
        fprintf('  DID NOT CONVERGE in %d steps (%.1f min).\n', n_steps, n_steps*dt/60);
    end
    fprintf('  Total Δv:           %.4f m/s\n', dv);
    fprintf('  Terminal distance:  %.4f m\n',   term_dist);
    fprintf('  Passive safe:       %.1f%%\n',   ps_frac * 100);
    if use_passive_safety && any(n_ps_h > 0)
        fprintf('  PS rows — avg: %.1f / max: %d\n', mean(n_ps_h), max(n_ps_h));
    end
    bad = find(flags(1:n_ctrl) ~= 1 & flags(1:n_ctrl) ~= 0);
    if isempty(bad)
        fprintf('  All QP solves:      optimal.\n');
    else
        fprintf('  Non-optimal QPs at: %s\n', mat2str(bad));
    end
    fprintf('─────────────────────────────────────────────────────────────\n\n');
end

%% Turnpike detection (plateau in log(||u||) for near-zero-Q cost)
u_thresh          = 1e-4 * u_max;
tp_mask           = u_norm_seq < u_thresh;
turnpike_steps    = find(tp_mask);
turnpike_fraction = length(turnpike_steps) / max(n_ctrl, 1);
tp_label          = 'Turnpike (near-zero thrust)';

if isempty(turnpike_steps) && n_ctrl >= 5
    log_u   = log(max(u_norm_seq, 1e-20));
    dlog_u  = abs(diff(log_u));
    pl_mask = [dlog_u < 0.15, false];
    pl_pad  = [false, pl_mask, false];
    d_pl    = diff(pl_pad);
    pl_st   = find(d_pl == 1);
    pl_en   = find(d_pl == -1) - 1;
    if ~isempty(pl_st)
        [~, li]           = max(pl_en - pl_st);
        turnpike_steps    = pl_st(li):pl_en(li);
        turnpike_fraction = length(turnpike_steps) / n_ctrl;
        tp_label          = 'Turnpike (steady cruise)';
    end
end

%% Pack result struct
res.X                 = X;
res.U                 = U;
res.ef                = ef;
res.flags             = flags;
res.is_safe           = is_safe;
res.dv                = dv;
res.term_dist         = term_dist;
res.ps_frac           = ps_frac;
res.conv_step         = conv_step;
res.u_norm_seq        = u_norm_seq;
res.all_X_pred        = all_X_pred;
res.turnpike_steps    = turnpike_steps;
res.turnpike_fraction = turnpike_fraction;
res.tp_label          = tp_label;
res.t_vec             = t_vec;
res.n_ctrl            = n_ctrl;
res.n_ps_h            = n_ps_h;

end
