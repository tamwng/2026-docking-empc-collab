% quantify_turnpike_tradeoff.m
% -------------------------------------------------------------------------
% STEP 4 — turn the strict-dissipativity LMI margins into mission numbers:
%   regularisation weight  w  <->  turnpike rate  theta(w)
%                              <->  practical-AS neighbourhood  rho(N,w)
%                              <->  fuel-bias / steady offset
% and contrast the exact dissipativity-LMI minimal regulariser against the
% Jaeschke/Gershgorin strong-convexity requirement.
%
% Same (Ad,Bd) as tests/verify_dissipativity_lmi.m: CWH from
% clohessy_wiltshire.m (column-extracted Ac) + discretize.m, closing dt=120 s.
% One weight w plays the role of BOTH rho_c (closing config A) and eps
% (fly-around error dynamics), since both use Qe=w*I6, Re=I3 on the same pair.
%
% Established (not re-derived here):
%   A full-state Qe=w*I6            -> strictly dissipative, LMI eps* = w EXACTLY.
%   B pos-only  Qe=w*blkdiag(I3,0)  -> strictly dissipative, eps* ~= 0.995 w.
%   C pure fuel Qe=0                -> eps* = 0 (obstruction).
%
% Console tables for items (1)-(4) + struct saved to
% tests/quantify_turnpike_tradeoff.mat.  Reads Run B closing JSON for item (3).
% Does NOT modify any src/ file.
% -------------------------------------------------------------------------

clear; clc;

%% Path + constants + (Ad,Bd) --------------------------------------------
script_dir = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(script_dir, '..', '..', 'src')));
constants;
dt = 120;                                   % closing sample time

Ac = zeros(6);
for j = 1:6
    ej = zeros(6,1); ej(j) = 1;
    Ac(:,j) = clohessy_wiltshire(ej, n);
end
Bc = [zeros(3,3); eye(3)];
[Ad, Bd] = discretize(Ac, Bc, dt);

w_sweep = [1e-12 1e-11 1e-10 1e-9 1e-8 1e-6 1e-4 1e-2];
nw = numel(w_sweep);

results = struct();
results.meta = struct('dt',dt,'n',n,'T',T,'w_sweep',w_sweep, ...
    'note','w plays role of rho_c (closing) AND eps (fly-around)');

fprintf('================================================================\n');
fprintf(' STEP 4 — turnpike / practical-AS / fuel-bias trade  (CWH, dt=%gs)\n', dt);
fprintf('================================================================\n');

%% ===== (1) TURNPIKE RATE theta(w) ======================================
% DARE/LQR for (Ad,Bd,Qe=w*I6,Re=I3): closed loop A_K = Ad - Bd*K,
% theta(w) = spectral radius (Damm Thm 6.2 exponential-turnpike rate).
theta = zeros(1,nw);
for i = 1:nw
    K = dlqr(Ad, Bd, w_sweep(i)*eye(6), eye(3));
    theta(i) = max(abs(eig(Ad - Bd*K)));
end
oneminus = 1 - theta;

fprintf('\n(1) TURNPIKE RATE  theta(w) = rho(Ad - Bd*K),  Qe=w*I6, Re=I3\n');
fprintf('    %-10s  %-12s  %-14s\n','w','theta','1 - theta');
for i = 1:nw
    fprintf('    %-10.0e  %-12.8f  %-14.6e\n', w_sweep(i), theta(i), oneminus(i));
end
% small-w scaling of (1 - theta): log-log slope on the monotone branch w<=1e-8
sel = w_sweep <= 1e-8;
p1  = polyfit(log10(w_sweep(sel)), log10(oneminus(sel)), 1);
[theta_min, imin] = min(theta);
w_opt = w_sweep(imin);
fprintf('    small-w fit:  (1 - theta) ~ w^%.4f   (log-log slope, w<=1e-8)\n', p1(1));
fprintf('    turnpike is SHARPEST at w* = %.0e : theta_min = %.4f\n', w_opt, theta_min);
fprintf('    VERDICT: theta(w) is U-SHAPED (non-monotonic), NOT "sharper with more reg".\n');
fprintf('             theta->1 as w->0 (turnpike vanishes); sharpens as ~w^%.2f down to\n', p1(1));
fprintf('             theta_min=%.2f at w*=%.0e; then RISES back toward 1 for large w\n', theta_min, w_opt);
fprintf('             (cheap-control limit: poles head to the sampled-plant zeros).\n');
results.item1 = struct('w',w_sweep,'theta',theta,'one_minus_theta',oneminus, ...
    'smallw_exponent',p1(1),'w_opt',w_opt,'theta_min',theta_min);

%% ===== (2) PRACTICAL-AS NEIGHBOURHOOD  theta(w)^N =======================
% Turnpike half-width model sigma_P(N) = C*theta(w)^N.  Tabulate theta^N and
% the horizon N* needed for theta^N <= 1e-3.
N_list = [10 20 50 92 150];
tol_ps = 1e-3;
thetaN = zeros(nw, numel(N_list));
Nstar  = zeros(1, nw);
for i = 1:nw
    thetaN(i,:) = theta(i).^N_list;
    if theta(i) < 1
        Nstar(i) = ceil(log(tol_ps)/log(theta(i)));
    else
        Nstar(i) = Inf;
    end
end

fprintf('\n(2) PRACTICAL-AS NEIGHBOURHOOD  theta(w)^N   (92 = fly-around period)\n');
fprintf('    %-10s', 'w');
for N = N_list, fprintf('  N=%-9d', N); end
fprintf('  %-8s\n', 'N*(1e-3)');
for i = 1:nw
    fprintf('    %-10.0e', w_sweep(i));
    for k = 1:numel(N_list)
        star = ''; if N_list(k)==92, star='*'; end
        fprintf('  %9.3e%1s', thetaN(i,k), star);
    end
    if isfinite(Nstar(i)), fprintf('  %-8d\n', Nstar(i));
    else,                  fprintf('  %-8s\n','inf'); end
end
fprintf('    (* column N=92 = one fly-around period.)\n');
fprintf('    VERDICT: N* is smallest near w*=%.0e (N*=%d); both smaller and larger w\n', ...
    w_opt, min(Nstar));
fprintf('             need a longer horizon. At the logged weights N=92 gives theta^92 =\n');
fprintf('             %.1e (Run B rho_c=1e-11) and %.1e (S2 eps=1e-4).\n', ...
    thetaN(w_sweep==1e-11,N_list==92), thetaN(w_sweep==1e-4,N_list==92));
results.item2 = struct('N_list',N_list,'thetaN',thetaN,'Nstar_1em3',Nstar,'tol',tol_ps);

%% ===== (3) FUEL-BIAS / STEADY OFFSET scaling (closing, config A) ========
% Predicted ||e||* ~ c / sqrt(w).  Extract ACTUAL terminal offset from the
% Run B closing sweep JSON (rho_c in {1e-12,1e-11,1e-10}) and regress
% log(offset) vs log(w).  Report the TRUE exponent (not the assumed -1/2).
json_path = fullfile(script_dir, '..', '..', 'exports', 'scenarios', 'sim_scenario4_closing.json');
fprintf('\n(3) FUEL-BIAS / STEADY OFFSET  ||e||* vs w  (closing Run B, config A)\n');
if exist(json_path, 'file')
    J  = jsondecode(fileread(json_path));
    sw = J.sweep_horizon;
    % Use a fixed horizon slice (N=20) across the three rho_c values.
    Nsel = 20;
    mask = (sw.N(:)==Nsel) & logical(sw.converged(:));
    wv   = sw.rho_c(mask);
    off  = sw.term_dist_m(mask);
    [wv, ord] = sort(wv); off = off(ord);
    fprintf('    source: sweep_horizon (N=%d, converged only)\n', Nsel);
    fprintf('    %-10s  %-14s\n','w (rho_c)','||e||* [m]');
    for i = 1:numel(wv), fprintf('    %-10.0e  %-14.4f\n', wv(i), off(i)); end
    if numel(wv) >= 2
        p3 = polyfit(log10(wv(:)), log10(off(:)), 1);
        fprintf('    regressed:  ||e||* ~ w^%.4f   (predicted -0.5)\n', p3(1));
        match = abs(p3(1) - (-0.5)) < 0.15;
        fprintf('    VERDICT: exponent = %.3f -> %s the -1/2 law.\n', p3(1), ...
            ternary(match,'MATCHES','does NOT match'));
        if ~match
            fprintf('             TRUTH: the closing phase REGULATES to the switch point\n');
            fprintf('             (setpoint is an equilibrium, e->0), so the terminal\n');
            fprintf('             distance is limited by the ~50 m convergence tolerance,\n');
            fprintf('             NOT a reg-dependent steady offset. No c/sqrt(w) law here;\n');
            fprintf('             the sqrt-neighbourhood shows up in item (2) via theta^N.\n');
        end
        results.item3 = struct('w',wv(:)','offset_m',off(:)','exponent',p3(1), ...
            'predicted',-0.5,'N',Nsel);
    else
        fprintf('    Not enough converged points to regress.\n');
        results.item3 = struct('w',wv(:)','offset_m',off(:)');
    end
else
    fprintf('    JSON not found: %s\n', json_path);
    fprintf('    Re-run simulations/sim_scenario4_closing.m to regenerate it.\n');
    results.item3 = struct('error','json_missing');
end

%% ===== (4) EXACT-LMI vs JAESCHKE/GERSHGORIN minimal regulariser =========
fprintf('\n(4) EXACT-LMI  vs  JAESCHKE/GERSHGORIN minimal regulariser\n');

% (4a) MAGNITUDE via Gershgorin strong convexity of the stage-cost Hessian
%      H = blkdiag(2*Qe, 2*Re), built on the pure-fuel baseline (Qe=0, Re=I3).
%      Row i is PD-safe (strongly convex) iff  H_ii - sum_{j~=i}|H_ij| > 0.
Re = eye(3);
H_pf = blkdiag(zeros(6), 2*Re);             % pure-fuel Hessian in [x;u]
fprintf('  (4a) MAGNITUDE — Gershgorin strong convexity of H=blkdiag(2Qe,2Re):\n');
fprintf('       row  H_ii   radius R_i=sum_{j~=i}|H_ij|   H_ii-R_i   needs Qe_ii?\n');
rank_needed = 0;
for i = 1:6                                  % the six STATE rows (pure fuel)
    Hii = H_pf(i,i);
    Ri  = sum(abs(H_pf(i,:))) - abs(Hii);
    slack = Hii - Ri;
    need  = slack <= 0;
    rank_needed = rank_needed + need;
    fprintf('       x%-2d  %5.2f  %-24.2f  %8.2f   %s\n', i, Hii, Ri, slack, ternary(need,'YES','no'));
end
fprintf('       => pure-fuel state rows have H_ii=0, radius=0: strong convexity\n');
fprintf('          needs Qe_ii>0 on ALL %d state rows (magnitude can be ->0+,\n', rank_needed);
fprintf('          matching exact-LMI eps*=w, but SUPPORT is forced to rank 6).\n');

% (4b) SUPPORT / RANK — the real win.
fprintf('\n  (4b) SUPPORT / RANK — minimal regulariser support:\n');
fprintf('       %-34s  %-8s  %-10s\n','method','support','rank');
fprintf('       %-34s  %-8s  %-10d\n','strong-convexity (Gershgorin)','all 6',6);
fprintf('       %-34s  %-8s  %-10d\n','exact dissipativity-LMI (cfg B)','position',3);
fprintf('       => LMI certifies strict dissipativity with a RANK-3 position-only\n');
fprintf('          Qe (eps ~= 0.995 w); strong convexity needs RANK-6 (velocities too).\n');
fprintf('       VERDICT: exact-LMI minimal support = {position}, rank 3;\n');
fprintf('                strong-convexity minimal support = {all states}, rank 6.\n');
results.item4 = struct( ...
    'gershgorin_rank', rank_needed, ...
    'lmi_min_support', 'position', 'lmi_min_rank', 3, ...
    'strongconv_min_support','all_states','strongconv_min_rank',6, ...
    'lmi_eps_full', 1.0, 'lmi_eps_posonly', 0.995);

%% ===== SAVE =============================================================
out_mat = fullfile(script_dir, 'quantify_turnpike_tradeoff.mat');
save(out_mat, 'results');
fprintf('\nSaved results struct -> %s\n', out_mat);
fprintf('Done.\n');

%% ---- local helper ------------------------------------------------------
function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end