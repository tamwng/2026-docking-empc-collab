% plot_turnpike.m
% -------------------------------------------------------------------------
% OPEN-LOOP turnpike visualisation for the regularised EMPC.
%
% Unlike run_closing.m / run_flyaround.m (which are receding-horizon), this
% script solves the finite-horizon condensed QP ONCE per horizon N and plots
% the resulting open-loop optimal state trajectory x_0..x_N.  That is the
% object the turnpike property is actually a statement about.
%
%   Fig 1 (closing)     d_k = ||x_k - x_s||          vs k/N,  N in {10..160}
%   Fig 2 (rate)        theta^k overlay on the N=160 approach,
%                       theta = rho(Ad - Bd*K) from DARE(Ad,Bd,rho_c*I6,I3)
%   Fig 3 (fly-around)  d_k = ||x_k - x*_k||         vs k/N,  N in {46..368}
%
% Closing:    ell = ||u||^2 + rho_c*||x - x_s||^2, rho_c=1e-11, dt=120 s,
%             x0 = [0;10000;0;0;0;0], x_s = [0;300;0;0;0;0], u_max = 1e-2.
% Fly-around: ell = ||u||^2 + eps*||x - x*_k||^2,   eps=1e-4,   dt=T/92,
%             x0 = [0;300;0;0;0;0], Pi* = 2:1 NMC ellipse, b = 75 m.
% No terminal ingredients in either case.
%
% Reuses the (Ad,Bd) builder and condensed-QP structure of
% simulations/sim_scenario4_closing.m / sim_scenario2_flyaround.m.
% Does NOT modify any src/ file.
% -------------------------------------------------------------------------

clear; clc;

%% Path + constants
script_dir = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(script_dir, '..', '..', 'src')));
addpath(genpath(fullfile(script_dir, '..', '..', 'tools', 'matlab2tikz-master', 'src')));
constants;                                  % provides n, T

out_dir = fullfile(script_dir, '..', '..', 'results', 'figures', 'turnpike');
[~, ~]  = mkdir(out_dir);

% The condensed Hessian is dense and very badly conditioned at the long
% horizons (cond(H) ~ 1e13 at N=368, 3N=1104 variables); interior-point stalls
% there.  The thrust bound is a pure box constraint, so active-set solves all
% horizons to optimality and is used throughout for comparability.
qp_opts = optimoptions('quadprog', 'Display', 'off', ...
    'Algorithm', 'active-set', 'MaxIterations', 5000);

%% ========================================================================
%  PART A — CLOSING PHASE
%  ========================================================================
dt_c    = 120;                              % [s]
u_max   = 1e-2;                             % [m/s^2]
rho_c   = 1e-11;
x0_c    = [0; 10000; 0; 0; 0; 0];
xs_c    = [0;   300; 0; 0; 0; 0];           % V-bar hold = CWH equilibrium
R_eco   = eye(3);

[Ad_c, Bd_c] = cwh_discrete(n, dt_c);

% x_s is an equilibrium of the CWH dynamics, so the error e = x - x_s obeys
% the same LTI recursion and the QP can be written directly in e-coordinates.
e0_c = x0_c - xs_c;

N_list_c = [10, 20, 40, 80, 160];
res_c    = cell(1, numel(N_list_c));

fprintf('=================================================================\n');
fprintf(' OPEN-LOOP TURNPIKE — closing phase  (dt=%g s, rho_c=%.0e)\n', dt_c, rho_c);
fprintf('=================================================================\n');
fprintf('  %5s  %12s  %12s  %12s  %10s\n', ...
    'N', 'min d_k [m]', 'plateau [m]', 'd_N [m]', 'frac<1%');

for ii = 1:numel(N_list_c)
    N = N_list_c(ii);
    [Sx, Su] = condensed_matrices(Ad_c, Bd_c, N);

    % ell = ||u||^2 + rho_c*||e||^2, states k=1..N (e_0 is constant)
    H = 2 * (kron(eye(N), R_eco) + rho_c * (Su' * Su));
    H = (H + H') / 2;
    f = 2 * rho_c * (Su' * (Sx * e0_c));

    lb = -u_max * ones(3*N, 1);
    ub =  u_max * ones(3*N, 1);

    [U_opt, ~, flag] = quadprog(H, f, [], [], [], [], lb, ub, zeros(3*N,1), qp_opts);
    assert(flag > 0, 'Closing QP failed at N=%d (exitflag %d)', N, flag);

    E = reshape(Sx * e0_c + Su * U_opt, 6, N);      % e_1..e_N
    d = [norm(e0_c), vecnorm(E, 2, 1)];             % d_0..d_N

    res_c{ii} = struct('N', N, 'd', d, 'U', U_opt, 'E', [e0_c, E]);

    plateau = median(d(max(1,round(0.3*N)) : round(0.7*N) + 1));
    frac1   = mean(d <= 0.01 * d(1));
    fprintf('  %5d  %12.3e  %12.3e  %12.3e  %9.1f%%\n', ...
        N, min(d), plateau, d(end), 100*frac1);
end

%% Turnpike rate from the DARE
[K_c, ~]  = dlqr(Ad_c, Bd_c, rho_c * eye(6), R_eco);
theta_c   = max(abs(eig(Ad_c - Bd_c * K_c)));

ref        = res_c{end};                     % N = 160
d_ref      = ref.d;
N_ref      = ref.N;
frac_1pct  = mean(d_ref <= 0.01 * d_ref(1));
plateau_lv = min(d_ref);

% Empirical decay rate: regress log(d_k) on k over the clean exponential band,
% i.e. after the thrust-saturated entry arc and before the numerical floor.
[theta_emp, k_fit_lo, k_fit_hi] = fit_decay_rate(d_ref, 1e-1, 1e-5);

fprintf('-----------------------------------------------------------------\n');
fprintf('  theta = rho(Ad - Bd*K), DARE(Ad,Bd,rho_c*I6,I3) : %.6f\n', theta_c);
fprintf('  1 - theta                                       : %.4e\n', 1 - theta_c);
fprintf('  empirical rate on N=%d, steps %d..%d             : %.6f\n', ...
    N_ref, k_fit_lo, k_fit_hi, theta_emp);
fprintf('  ratio theta_emp / theta                         : %.4f\n', theta_emp/theta_c);
fprintf('  N=%d: fraction of horizon within 1%% of ||x0-xs||: %.1f%%  (%d of %d steps)\n', ...
    N_ref, 100*frac_1pct, round(frac_1pct*(N_ref+1)), N_ref+1);
fprintf('  N=%d: plateau level min_k d_k                   : %.4e m\n', N_ref, plateau_lv);
fprintf('  N=%d: d_0 = %.1f m,  d_N = %.4e m\n', N_ref, d_ref(1), d_ref(end));
fprintf('-----------------------------------------------------------------\n\n');

%% Figure 1 — the turnpike (closing)
cmap_c = turnpike_ramp(numel(N_list_c));

fig1 = figure('Name', 'Turnpike — closing (open-loop)');
ax1  = axes(fig1); hold(ax1, 'on'); grid(ax1, 'on');
set(ax1, 'YScale', 'log');
for ii = 1:numel(N_list_c)
    d = max(res_c{ii}.d, 1e-16);
    plot(ax1, (0:res_c{ii}.N) / res_c{ii}.N, d, '-', ...
        'Color', cmap_c(ii,:), 'LineWidth', 1.4, ...
        'DisplayName', sprintf('$N=%d$', res_c{ii}.N));
end
yline(ax1, 0.01 * d_ref(1), ':k', 'LineWidth', 0.9, ...
    'HandleVisibility', 'off');
text(ax1, 0.02, 0.01*d_ref(1), '$1\%$ of $\|x_0-x_s\|$', ...
    'Interpreter', 'latex', 'FontSize', 8, 'VerticalAlignment', 'bottom');
xlabel(ax1, 'Normalised time $k/N$', 'Interpreter', 'latex');
ylabel(ax1, '$d_k = \|x_k - x_s\|$ [m]', 'Interpreter', 'latex');
title(ax1, sprintf(['Open-loop turnpike --- closing, $\\ell=\\|u\\|^2+' ...
    '\\rho_c\\|x-x_s\\|^2$, $\\rho_c=10^{%d}$'], round(log10(rho_c))), ...
    'Interpreter', 'latex');
legend(ax1, 'Interpreter', 'latex', 'Location', 'southwest');

%% Figure 2 — the rate
fig2 = figure('Name', 'Turnpike rate — theta^k overlay');
ax2  = axes(fig2); hold(ax2, 'on'); grid(ax2, 'on');
set(ax2, 'YScale', 'log');

k_ref = 0:N_ref;
plot(ax2, k_ref, max(d_ref, 1e-16), '-', 'Color', [0.12 0.47 0.71], ...
    'LineWidth', 1.6, 'DisplayName', sprintf('$d_k$, $N=%d$', N_ref));

% Reference decay at the DARE turnpike rate, anchored at the start of the
% exponential band (the first ~20 steps run thrust-saturated, so anchoring at
% k=0 would compare against an arc the rate does not describe).
k_ov  = k_fit_lo : N_ref;
d_anc = d_ref(k_fit_lo + 1);
plot(ax2, k_ov, d_anc * theta_c.^(k_ov - k_fit_lo), '--', ...
    'Color', [0.85 0.33 0.10], 'LineWidth', 1.4, ...
    'DisplayName', sprintf('$\\theta^k$ (DARE), $\\theta=%.4f$', theta_c));
plot(ax2, k_ov, d_anc * theta_emp.^(k_ov - k_fit_lo), '-.', ...
    'Color', [0.20 0.63 0.17], 'LineWidth', 1.2, ...
    'DisplayName', sprintf('fitted, $\\hat\\theta=%.4f$', theta_emp));
xline(ax2, k_fit_lo, ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
xline(ax2, k_fit_hi, ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');

yline(ax2, plateau_lv, ':', 'LineWidth', 0.9, 'Color', [0.35 0.35 0.35], ...
    'HandleVisibility', 'off');
text(ax2, 0.55*N_ref, plateau_lv, ...
    sprintf('plateau $\\approx %.1e$ m', plateau_lv), ...
    'Interpreter', 'latex', 'FontSize', 8, 'VerticalAlignment', 'bottom');
xlabel(ax2, 'Step $k$', 'Interpreter', 'latex');
ylabel(ax2, '$d_k = \|x_k - x_s\|$ [m]', 'Interpreter', 'latex');
title(ax2, sprintf(['Approach rate vs.\\ DARE turnpike rate ($N=%d$, ' ...
    '$%.0f\\%%$ of horizon on the turnpike)'], N_ref, 100*frac_1pct), ...
    'Interpreter', 'latex');
legend(ax2, 'Interpreter', 'latex', 'Location', 'southwest');

%% ========================================================================
%  PART B — FLY-AROUND (NMC orbit Pi*)
%  ========================================================================
P_orb  = 92;
dt_f   = T / P_orb;
eps_f  = 1e-4;
b_nmc  = 75;                                % [m] radial semi-axis
x0_f   = [0; 300; 0; 0; 0; 0];              % V-bar hold

[Ad_f, Bd_f] = cwh_discrete(n, dt_f);

closure_res = norm(Ad_f^P_orb * [100;0;0;0;-2*n*100;0] - [100;0;0;0;-2*n*100;0]);
assert(closure_res < 1e-4, 'Discrete NMT does not close — check dt and P');

% Prescribed periodic orbit Pi*: propagate the on-manifold IC one full period
x_star      = zeros(6, P_orb);
x_star(:,1) = [b_nmc; 0; 0; 0; -2*n*b_nmc; 0];
for k = 1:P_orb-1
    x_star(:, k+1) = Ad_f * x_star(:, k);
end
assert(norm(Ad_f * x_star(:,end) - x_star(:,1)) < 1e-6, 'Pi* does not close');

% x*_k is a trajectory of the free dynamics, so e_k = x_k - x*_k again obeys
% e_{k+1} = Ad e_k + Bd u_k — identical condensed structure to Part A.
e0_f = x0_f - x_star(:,1);

N_list_f = [46, 92, 184, 368];
res_f    = cell(1, numel(N_list_f));

fprintf('=================================================================\n');
fprintf(' OPEN-LOOP TURNPIKE — fly-around  (dt=T/%d=%.1f s, eps=%.0e, b=%d m)\n', ...
    P_orb, dt_f, eps_f, b_nmc);
fprintf('=================================================================\n');
fprintf('  %5s  %8s  %12s  %12s  %12s  %10s\n', ...
    'N', 'N/P', 'min d_k [m]', 'plateau [m]', 'd_N [m]', 'frac<1%');

for ii = 1:numel(N_list_f)
    N = N_list_f(ii);
    [Sx, Su] = condensed_matrices(Ad_f, Bd_f, N);

    H = 2 * (kron(eye(N), R_eco) + eps_f * (Su' * Su));
    H = (H + H') / 2;
    f = 2 * eps_f * (Su' * (Sx * e0_f));

    lb = -u_max * ones(3*N, 1);
    ub =  u_max * ones(3*N, 1);

    [U_opt, ~, flag] = quadprog(H, f, [], [], [], [], lb, ub, zeros(3*N,1), qp_opts);
    assert(flag > 0, 'Fly-around QP failed at N=%d (exitflag %d)', N, flag);

    E = reshape(Sx * e0_f + Su * U_opt, 6, N);
    d = [norm(e0_f), vecnorm(E, 2, 1)];

    res_f{ii} = struct('N', N, 'd', d, 'U', U_opt);

    plateau = median(d(max(1,round(0.3*N)) : round(0.7*N) + 1));
    frac1   = mean(d <= 0.01 * d(1));
    fprintf('  %5d  %8.1f  %12.3e  %12.3e  %12.3e  %9.1f%%\n', ...
        N, N/P_orb, min(d), plateau, d(end), 100*frac1);
end

[K_f, ~] = dlqr(Ad_f, Bd_f, eps_f * eye(6), R_eco);
theta_f  = max(abs(eig(Ad_f - Bd_f * K_f)));
ref_f    = res_f{end};
frac_f   = mean(ref_f.d <= 0.01 * ref_f.d(1));

fprintf('-----------------------------------------------------------------\n');
% Rate is fitted on the N=P curve: at N=2P and 4P the approach reaches the
% QP solver's accuracy floor (~1e-5 m) less than a third into the horizon, so
% the rest of those curves is numerical noise, not dynamics.
i_rate      = find(N_list_f == P_orb, 1);
d_rate      = res_f{i_rate}.d;
[theta_emp_f, kf_lo, kf_hi] = fit_decay_rate(d_rate, 1e0, 1e-7);
fprintf('  theta_flyaround = rho(Ad - Bd*K), DARE(Ad,Bd,eps*I6,I3) : %.6f\n', theta_f);
fprintf('  empirical rate on N=%d, steps %d..%d                    : %.6f\n', ...
    P_orb, kf_lo, kf_hi, theta_emp_f);
fprintf('  ratio theta_emp / theta                                 : %.4f\n', ...
    theta_emp_f / theta_f);
fprintf('  NOTE: the N=%d and N=%d tails sit on the QP accuracy floor (~1e-5 m);\n', ...
    2*P_orb, 4*P_orb);
fprintf('        their plateau level is solver noise, not a dynamic limit.\n');
fprintf('  N=%d: fraction within 1%% of ||x0-x*_0||                : %.1f%%\n', ...
    ref_f.N, 100*frac_f);
fprintf('  N=%d: plateau level min_k d_k                          : %.4e m\n', ...
    ref_f.N, min(ref_f.d));
fprintf('-----------------------------------------------------------------\n\n');

%% Figure 3 — the turnpike (fly-around)
cmap_f = turnpike_ramp(numel(N_list_f));

fig3 = figure('Name', 'Turnpike — fly-around (open-loop)');
ax3  = axes(fig3); hold(ax3, 'on'); grid(ax3, 'on');
set(ax3, 'YScale', 'log');
for ii = 1:numel(N_list_f)
    d = max(res_f{ii}.d, 1e-16);
    plot(ax3, (0:res_f{ii}.N) / res_f{ii}.N, d, '-', ...
        'Color', cmap_f(ii,:), 'LineWidth', 1.4, ...
        'DisplayName', sprintf('$N=%d$ ($%.1f\\times P$)', ...
            res_f{ii}.N, res_f{ii}.N / P_orb));
end
ylim(ax3, [1e-8, 1e4]);
yline(ax3, 0.01 * ref_f.d(1), ':k', 'LineWidth', 0.9, 'HandleVisibility', 'off');
text(ax3, 0.02, 0.01*ref_f.d(1), '$1\%$ of $\|x_0-x^*_0\|$', ...
    'Interpreter', 'latex', 'FontSize', 8, 'VerticalAlignment', 'bottom');
xlabel(ax3, 'Normalised time $k/N$', 'Interpreter', 'latex');
ylabel(ax3, '$d_k = \|x_k - x^*_k\|_{\Pi^*}$ [m]', 'Interpreter', 'latex');
title(ax3, sprintf(['Open-loop turnpike --- fly-around onto $\\Pi^*$, ' ...
    '$\\ell=\\|u\\|^2+\\varepsilon\\|x-x^*_k\\|^2$, $\\varepsilon=10^{%d}$'], ...
    round(log10(eps_f))), 'Interpreter', 'latex');
legend(ax3, 'Interpreter', 'latex', 'Location', 'southwest');

%% Export
print(fig1, fullfile(out_dir, 'turnpike_closing'),        '-dpng', '-r150');
print(fig2, fullfile(out_dir, 'turnpike_rate_closing'),   '-dpng', '-r150');
print(fig3, fullfile(out_dir, 'turnpike_flyaround'),      '-dpng', '-r150');
if exist('matlab2tikz', 'file')
    matlab2tikz(fullfile(out_dir, 'turnpike_closing.tikz'), ...
        'figurehandle', fig1, 'showInfo', false);
    matlab2tikz(fullfile(out_dir, 'turnpike_rate_closing.tikz'), ...
        'figurehandle', fig2, 'showInfo', false);
    matlab2tikz(fullfile(out_dir, 'turnpike_flyaround.tikz'), ...
        'figurehandle', fig3, 'showInfo', false);
    fprintf('TikZ + PNG figures written to %s\n', out_dir);
else
    fprintf('PNG figures written to %s  (TikZ skipped — run setup.m first)\n', out_dir);
end

%% ---- local helpers ------------------------------------------------------
function [Ad, Bd] = cwh_discrete(n, dt)
% Same CWH continuous pair + ZOH discretisation used across the simulations.
    Ac = [0      0     0    1     0    0  ;
          0      0     0    0     1    0  ;
          0      0     0    0     0    1  ;
          3*n^2  0     0    0     2*n  0  ;
          0      0     0   -2*n   0    0  ;
          0      0    -n^2  0     0    0 ];
    Bc = [zeros(3,3); eye(3)];
    [Ad, Bd] = discretize(Ac, Bc, dt);
end

function [Sx, Su] = condensed_matrices(Ad, Bd, N)
% Stacked prediction e_1..e_N = Sx*e_0 + Su*U   (Sx: 6N x 6, Su: 6N x 3N).
    Sx = zeros(6*N, 6);
    Su = zeros(6*N, 3*N);
    Ad_pow = eye(6);
    for i = 1:N
        Ad_pow = Ad_pow * Ad;
        Sx((i-1)*6+1 : i*6, :) = Ad_pow;
        for j = 1:i
            Su((i-1)*6+1 : i*6, (j-1)*3+1 : j*3) = Ad^(i-j) * Bd;
        end
    end
end

function [theta_hat, k_lo, k_hi] = fit_decay_rate(d, d_hi, d_lo)
% Least-squares rate on the clean exponential band d_lo <= d_k <= d_hi:
% log d_k ~ a + k*log(theta_hat).  Band edges returned as step indices.
    idx = find(d <= d_hi & d >= d_lo);
    if numel(idx) < 5
        theta_hat = NaN; k_lo = NaN; k_hi = NaN; return;
    end
    k_lo = idx(1) - 1;            % d is indexed d(1) = d_0
    k_hi = idx(end) - 1;
    kk   = (k_lo:k_hi)';
    p    = polyfit(kk, log(d(k_lo+1 : k_hi+1))', 1);
    theta_hat = exp(p(1));
end

function c = turnpike_ramp(m)
% Blue (short horizon) -> dark red (long horizon).
    c = zeros(m, 3);
    for i = 1:m
        t = (i - 1) / max(m - 1, 1);
        c(i,:) = (1 - t) * [0.12 0.47 0.71] + t * [0.65 0.08 0.08];
    end
end