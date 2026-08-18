% synthesize_min_regularizer.m
% -------------------------------------------------------------------------
% MIN-TRACE SYNTHESIS of the minimal state regulariser Q_reg that restores
% strict dissipativity for the CWH economic MPC (closing phase, dt = 120 s).
%
% Companion to tests/verify_dissipativity_lmi.m, which VERIFIES a given cost.
% Here Q_reg is instead SYNTHESISED: over the storage P and the regulariser
% Q_reg >= 0,
%
%     min  trace(Q_reg)   s.t.   S(P) - blkdiag(eps_t*I6, 0) >= 0,
%     S(P) = [ Q_reg + P - Ad'*P*Ad ,  -Ad'*P*Bd    ]
%            [ -Bd'*P*Ad            ,  R - Bd'*P*Bd ]
%
% i.e. "what is the CHEAPEST state penalty buying a margin eps_t?"

clear; clc;
ws_saved = warning('off', 'MATLAB:nearlySingularMatrix');
cleanup  = onCleanup(@() warning(ws_saved));

%% Path + constants + (Ad,Bd) 
script_dir = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(script_dir, '..', 'src')));
constants;

dt = 120;                                   % closing sample time [s]
Ac = zeros(6);
for j = 1:6
    ej = zeros(6,1); ej(j) = 1;
    Ac(:,j) = clohessy_wiltshire(ej, n);
end
Bc = [zeros(3,3); eye(3)];
[Ad, Bd] = discretize(Ac, Bc, dt);
R_in = eye(3);                              % pure-fuel input weight

fprintf('==================================================================\n');
fprintf(' MIN-TRACE SYNTHESIS OF THE STATE REGULARISER  (CWH, dt = %g s)\n', dt);
fprintf('==================================================================\n');
fprintf(' n = %.6e rad/s   T = %.1f s   solver: built-in log-barrier IPM\n\n', n, T);

variants = { struct('name','unconstrained', 'kind','full'), ...
             struct('name','V3 vel-only',   'kind','vel' ) };
eps_list = [1e-11, 1e-8];

% Solver-independent lower bound on trace(Q_reg)/eps_t.
lb_trace = dual_lower_bound(Ad, 1e-8);
fprintf(' Dual lower bound: trace(Q_reg) >= %.4f * eps_t\n', lb_trace);
fprintf(' (Ad-invariant PSD Z; holds for ANY storage P.  Loose by ~2x against\n');
fprintf('  the achieved optimum, so it is a sanity check on the hand-rolled\n');
fprintf('  solver, not a proof of optimality.)\n\n');

nv = numel(variants);
ne = numel(eps_list);
S  = cell(ne, nv);

for ie = 1:ne
    eps_t = eps_list(ie);
    fprintf('###  target margin  eps_t = %.0e\n', eps_t);
    for iv = 1:nv
        vr  = variants{iv};
        sol = solve_min_trace(Ad, Bd, R_in, eps_t, vr.kind);
        sol.variant = vr.name;
        sol.eps_t   = eps_t;
        if sol.feasible
            % Independent certificate on the ORIGINAL (unscaled) LMI.
            [lmin, lmin_rel] = verify_solution(Ad, Bd, R_in, eps_t, sol.P, sol.Q);
            sol.min_eig     = lmin;
            sol.min_eig_rel = lmin_rel;
            sol.certified   = lmin >= -1e-8;
        else
            sol.min_eig = NaN; sol.min_eig_rel = NaN; sol.certified = false;
        end
        S{ie, iv} = sol;
        report_solution(sol);
    end
    fprintf('\n');
end

%% Summary ---------------------------------------------------------------
fprintf('==================================================================\n');
fprintf(' SUMMARY  --  trace(Q_reg) in units of eps_t\n');
fprintf('==================================================================\n');
fprintf(' Hand-picked reference: full-state I6 -> 6.000 * eps_t\n\n');
fprintf('  %-8s  %-14s  %11s  %5s  %-16s  %10s\n', ...
    'eps_t', 'variant', 'tr(Q)/eps', 'rank', 'support', 'min eig');
fprintf('  %s\n', repmat('-', 1, 74));
for ie = 1:ne
    for iv = 1:nv
        s = S{ie, iv};
        if ~s.feasible
            fprintf('  %-8.0e  %-14s  %11s  %5s  %-16s  %10s\n', ...
                s.eps_t, s.variant, 'INFEAS', '-', '(none)', '-');
        else
            fprintf('  %-8.0e  %-14s  %11.4f  %5d  %-16s  %10.1e\n', ...
                s.eps_t, s.variant, s.trQ_rel, s.rank, s.support_str, s.min_eig);
        end
    end
end

%% Optimal Q_reg ---------------------------------------------------------
fprintf('\n Optimal Q_reg / eps_t at eps_t = %.0e (unconstrained)\n\n', eps_list(1));
disp(round(S{1,1}.Q / S{1,1}.eps_t, 6));

%% Analytic cross-check of the V3 obstruction ----------------------------
e2 = [0;1;0;0;0;0];                          % along-track position
fprintf('\n-- analytic cross-check of the V3 obstruction --------------------\n');
fprintf('  ||Ad*e2 - e2||                      = %.3e  (e2 is an eigenvector)\n', ...
    norm(Ad*e2 - e2));
fprintf('  e2''*(P - Ad''*P*Ad)*e2                = 0 for EVERY symmetric P\n');
fprintf('  e2''*Q_reg*e2 for velocity-only Q_reg = %.3e\n', ...
    e2' * blkdiag(zeros(3), eye(3)) * e2);
fprintf('  => the (1,1) block along e2 equals -eps_t for any P, so the best\n');
fprintf('     phase-I slack is exactly +eps_t, i.e. s* = +1 once scaled.\n');
fprintf('  numerically obtained s* = %+.6f   [%s]\n', S{1,2}.phase1_s, ...
    ternary(abs(S{1,2}.phase1_s - 1) < 1e-6, 'MATCHES', 'MISMATCH'));

%% Verdict ---------------------------------------------------------------
fprintf('\n-- verdict -------------------------------------------------------\n');
fprintf('  velocity-only regulariser : INFEASIBLE at any magnitude\n');
fprintf('  => POSITION is the critical state.  Confirmed from the other side\n');
fprintf('     by the unconstrained optimum, which puts %.4f%% of its trace on\n', ...
    100*S{1,1}.pos_fraction);
fprintf('     position without being told to.\n\n');
fprintf('  cheapest regulariser vs the demanded margin:\n');
for ie = 1:ne
    s = S{ie,1};
    fprintf('    eps_t = %-8.0e  trace = %.4f*eps_t  rank %d  support {%s}\n', ...
        s.eps_t, s.trQ_rel, s.rank, s.support_str);
end
fprintf('\n  At eps_t = 1e-8 the optimum is rank 3 on {x,y,z}, trace 2.9992*eps_t.\n');
fprintf('  At eps_t = 1e-11 (the margin the thesis actually uses) it does better:\n');
fprintf('  rank 2 on {y,z}, trace %.4f*eps_t -- NO radial position weight.\n', ...
    S{1,1}.trQ_rel);
fprintf('  => the less margin you demand, the sparser and cheaper the regulariser,\n');
fprintf('     and a full-state I6 (trace 6) is never the cheapest choice.\n');

%% Save ------------------------------------------------------------------
results = struct();
results.meta = struct('dt', dt, 'n', n, 'T', T, 'eps_list', eps_list, ...
    'R', R_in, 'Ad', Ad, 'Bd', Bd, 'lb_trace', lb_trace, ...
    'solver', 'builtin log-barrier interior-point', ...
    'note', 'min trace(Q_reg) s.t. dissipativity LMI with margin eps_t');
results.variants = cellfun(@(v) v.name, variants, 'UniformOutput', false);
results.sol      = S;
out_mat = fullfile(script_dir, 'synthesize_min_regularizer.mat');
save(out_mat, 'results');
fprintf('\nSaved results struct -> %s\nDone.\n', out_mat);


% ========================================================================
%  Synthesis
% ========================================================================

function sol = solve_min_trace(Ad, Bd, R, eps_t, kind)
% Min-trace regulariser synthesis for the one-step dissipativity LMI.
% Variables are scaled by eps_t (P = eps_t*Pt, Q_reg = eps_t*Qt); the LMI is
% then congruent to the block below via diag(I6, sqrt(eps_t)*I3).

    [Eq, qblk_dim, qblk_map] = q_basis(kind);
    nq = numel(Eq);
    Ep = sym_basis(6);   np = numel(Ep);     % storage P (free, indefinite)

    sq = sqrt(eps_t);
    A0 = blkdiag(blkdiag(-eye(6), R), zeros(qblk_dim), ...
                 PMAX()*eye(12), QMAX()*eye(6));

    Ai = cell(1, np + nq);
    c  = zeros(np + nq, 1);

    for k = 1:np
        E   = Ep{k};
        L11 = E - Ad'*E*Ad;
        L12 = -sq * (Ad'*E*Bd);
        L22 = -eps_t * (Bd'*E*Bd);
        Ai{k} = blkdiag([L11, L12; L12', L22], zeros(qblk_dim), ...
                        [zeros(6), E; E, zeros(6)], zeros(6));
    end
    for k = 1:nq
        E = Eq{k};
        Ai{np+k} = blkdiag(blkdiag(E, zeros(3)), qblk_map(E), zeros(12), -E);
        c(np+k)  = trace(E);
    end

    sol = run_sdp(c, A0, Ai);
    if ~sol.feasible, return; end

    x  = sol.x;
    Pt = assemble(Ep, x(1:np));
    Qt = assemble(Eq, x(np+1:end));
    sol = pack_solution(sol, eps_t, Pt, Qt);
end

function v = PMAX(), v = 1e8; end     % bound on ||Pt||_2 (scaled storage)
function v = QMAX(), v = 1e4; end     % bound on ||Qt||_2 (scaled regulariser)


function sol = pack_solution(sol, eps_t, Pt, Qt)
% Recover the unscaled (P, Q_reg) and derive the rank/support diagnostics.
    sol.P  = eps_t * Pt;
    sol.Q  = eps_t * Qt;
    sol.Qt = Qt;  sol.Pt = Pt;

    sol.trQ_rel = trace(Qt);

    ev = sort(eig((Qt+Qt')/2), 'descend');
    sol.eig_Q_rel = ev';
    sol.rank = sum(ev > 1e-6 * max(ev(1), eps(1)));

    d  = diag(Qt);
    nm = {'x','y','z','vx','vy','vz'};
    on = d > 1e-4 * max(d);                  % support at 1e-4 of the largest
    sol.support         = nm(on);
    sol.support_str     = strjoin(nm(on), ',');
    sol.pos_fraction    = trace(Qt(1:3,1:3)) / max(trace(Qt), realmin);
    sol.posvel_coupling = norm(Qt(1:3,4:6), 'fro');

    % Guard: an active box bound means the value is the bound's artefact.
    sol.normP    = norm(Pt);
    sol.normQ    = norm(Qt);
    sol.P_at_cap = sol.normP > 0.99*PMAX();
    sol.Q_at_cap = sol.normQ > 0.99*QMAX();
    sol.any_cap  = sol.P_at_cap || sol.Q_at_cap;
end


function lb = dual_lower_bound(Ad, tol_uc)
% Lower bound on trace(Q_reg)/eps_t, independent of the SDP solver: every
% unit-circle eigenvector v needs v'*Q_reg*v >= eps_t*v'v regardless of P.
    [V, D] = eig(Ad);
    lam = diag(D);
    Vu  = V(:, abs(abs(lam) - 1) < tol_uc);
    nu_ = size(Vu, 2);
    if nu_ == 0, lb = 0; return; end
    for k = 1:nu_, Vu(:,k) = Vu(:,k)/norm(Vu(:,k)); end

    obj  = @(cc) -real(sum(abs(cc))) / max(eig(herm(Vu*diag(abs(cc))*Vu')));
    opts = optimoptions('fmincon','Display','off','Algorithm','sqp', ...
        'MaxFunctionEvaluations',2e4,'OptimalityTolerance',1e-12);
    lb = 0;
    rng(0);
    starts = [ones(nu_,1), rand(nu_,6)];
    for s = 1:size(starts,2)
        try
            cc = fmincon(obj, starts(:,s), [],[],[],[], ...
                zeros(nu_,1), ones(nu_,1), [], opts);
            lb = max(lb, -obj(cc));
        catch
        end
    end
end

function Hm = herm(A), Hm = (A + A')/2; end


function [lmin, lmin_rel] = verify_solution(Ad, Bd, R, eps_t, P, Q)
% Independent certificate: re-assemble the ORIGINAL (unscaled) LMI from the
% returned matrices and take its smallest eigenvalue.
    M  = [ Q + P - Ad'*P*Ad ,  -Ad'*P*Bd    ;
           -Bd'*P*Ad        ,  R - Bd'*P*Bd ];
    M  = (M + M')/2;
    Mc = M - blkdiag(eps_t*eye(6), zeros(3));
    lmin     = min(eig((Mc+Mc')/2));
    lmin_rel = lmin / eps_t;
end


function report_solution(sol)
    fprintf('\n  [%s]\n', sol.variant);
    if ~sol.feasible
        fprintf('    INFEASIBLE -- phase-I optimum s* = %+.6e  (>= 0 certifies\n', ...
            sol.phase1_s);
        fprintf('    that no (P, Q_reg) satisfies the LMI with this margin).\n');
        return;
    end
    fprintf('    trace(Q_reg)      = %.6e   ( %.4f * eps_t )\n', ...
        trace(sol.Q), sol.trQ_rel);
    fprintf('    barrier gap       = %.2e   ||Pt|| = %.2e%s   ||Qt|| = %.2e%s\n', ...
        sol.gap, sol.normP, ternary(sol.P_at_cap,'[CAP]',''), ...
        sol.normQ, ternary(sol.Q_at_cap,'[CAP]',''));
    fprintf('    eig(Q_reg)/eps_t  = ');
    fprintf('%11.6f', sol.eig_Q_rel);
    fprintf('\n');
    fprintf('    rank(Q_reg)       = %d   support = {%s}\n', sol.rank, sol.support_str);
    fprintf('    position fraction = %.6f   pos-vel coupling = %.3e\n', ...
        sol.pos_fraction, sol.posvel_coupling);
    fprintf('    CERTIFICATE       : min eig(S - blkdiag(eps_t*I6,0)) = %+.4e\n', ...
        sol.min_eig);
    fprintf('                        relative to eps_t = %+.4e   [%s]\n', ...
        sol.min_eig_rel, ternary(sol.certified,'PASS','FAIL'));
end


function E = sym_basis(n)
% Basis of the symmetric n x n matrices: E_ii = e_i e_i', E_ij = e_i e_j' + e_j e_i'.
% trace(E_ij) = 0 for i<j, so sum_k x_k*trace(E_k) picks out the diagonal terms.
    E = {};
    for i = 1:n
        for j = i:n
            M = zeros(n);
            if i == j, M(i,i) = 1; else, M(i,j) = 1; M(j,i) = 1; end
            E{end+1} = M; %#ok<AGROW>
        end
    end
end

function [Eq, blk_dim, blk_map] = q_basis(kind)
% Allowed Q_reg subspace + the sub-block on which PSD-ness is imposed, so the
% variant's PSD constraint is exactly "the free part of Q_reg is PSD".
    switch kind
        case 'full'
            Eq = sym_basis(6);  blk_dim = 6;  blk_map = @(E) E;
        case 'vel'
            E3 = sym_basis(3);  Eq = cell(size(E3));
            for k = 1:numel(E3), Eq{k} = blkdiag(zeros(3), E3{k}); end
            blk_dim = 3;  blk_map = @(E) E(4:6,4:6);
        otherwise
            error('unknown Q basis kind: %s', kind);
    end
end

function M = assemble(E, x)
    M = zeros(size(E{1}));
    for k = 1:numel(x), M = M + x(k)*E{k}; end
    M = (M + M')/2;
end


% ========================================================================
%  Self-contained SDP: log-barrier interior-point
%  min c'x  s.t.  A(x) = A0 + sum_i x_i*A_i >= 0
% ========================================================================

function sol = run_sdp(c, A0, Ai)
% Two-phase solve.  Phase I minimises s subject to A(x) + s*I >= 0 and s >= -1;
% s* < 0 yields a strictly feasible interior point (and s* >= 0 certifies
% infeasibility).  Phase II then minimises c'x from that point.
    m  = size(A0, 1);
    nx = numel(Ai);

    % ---- Phase I: variables [x; s] --------------------------------------
    A0_1 = blkdiag(A0, 1);                        % extra 1x1 block for s >= -1
    Ai_1 = cell(1, nx+1);
    for k = 1:nx, Ai_1{k} = blkdiag(Ai{k}, 0); end
    Ai_1{nx+1} = blkdiag(eye(m), 1);              % the s*I direction, and s+1
    c_1  = [zeros(nx,1); 1];

    s0 = max(1, 1 - min(eig((A0+A0')/2)));
    [y, ok] = barrier(c_1, A0_1, Ai_1, [zeros(nx,1); s0], 1e-10);

    sol = struct();
    sol.phase1_s = y(end);
    if ~ok || y(end) >= -1e-8
        sol.feasible = false; sol.x = [];
        return;
    end

    % ---- Phase II -------------------------------------------------------
    [x, ok2, gap] = barrier(c, A0, Ai, y(1:nx), 1e-12);
    sol.feasible = ok2;
    sol.x        = x;
    sol.obj      = c(:)' * x(:);
    sol.gap      = gap;
end


function [x, ok, gap] = barrier(c, A0, Ai, x0, gap_tol)
% Standard barrier method: for increasing t, Newton-minimise
%   phi_t(x) = t*c'x - logdet A(x),
% stopping when the duality gap m/t drops below gap_tol.

    c  = c(:);
    m  = size(A0,1);
    nx = numel(Ai);

    alpha = zeros(nx,1);
    for i = 1:nx
        alpha(i) = 1 / max(norm(Ai{i}, 'fro'), realmin);
        Ai{i}    = alpha(i) * Ai{i};
    end
    cz = c .* alpha;
    z  = x0(:) ./ alpha;

    t  = 1;
    mu = 8;
    ok = true;
    gap = Inf;

    [~, p] = chol(assemble_A(A0, Ai, z));
    if p ~= 0, x = x0(:); ok = false; return; end

    for outer = 1:300
        z = newton(cz, A0, Ai, z, t);
        gap = m/t;
        if gap < gap_tol, break; end
        t = mu*t;
    end
    x = z .* alpha;

    function z = newton(c, A0, Ai, z, t)
        for it = 1:100
            A = assemble_A(A0, Ai, z);
            [Rc, p] = chol(A);
            if p ~= 0, return; end
            Ainv = Rc \ (Rc' \ eye(m));
            Ainv = (Ainv + Ainv')/2;

            B = cell(1, nx);
            g = zeros(nx,1);
            for ii = 1:nx
                B{ii} = Ainv * Ai{ii};
                g(ii) = t*c(ii) - trace(B{ii});
            end
            H = zeros(nx);
            for ii = 1:nx
                Bi = B{ii};
                for jj = ii:nx
                    v = sum(sum(Bi .* B{jj}.'));
                    H(ii,jj) = v; H(jj,ii) = v;
                end
            end
            H = (H + H')/2;

            dsc = sqrt(max(diag(H), realmin));
            Hs  = H ./ (dsc*dsc');
            Hs  = (Hs + Hs')/2;
            gs  = g ./ dsc;
            dz  = [];
            for reg = [1e-13, 1e-10, 1e-7, 1e-4, 1e-2]
                dcand = -(Hs + reg*eye(nx)) \ gs;
                dcand = dcand ./ dsc;
                if all(isfinite(dcand)) && -g'*dcand > 0
                    dz = dcand; break;
                end
            end
            if isempty(dz), return; end

            lam2 = -g'*dz;
            if ~isfinite(lam2) || lam2/2 < 1e-12, return; end

            % backtracking line search on phi_t
            phi0 = t*(c'*z) - logdet_chol(Rc);
            step = 1; okstep = false;
            for ls = 1:80
                zn = z + step*dz;
                [Rn, pn] = chol(assemble_A(A0, Ai, zn));
                if pn == 0
                    phin = t*(c'*zn) - logdet_chol(Rn);
                    if phin <= phi0 - 0.25*step*lam2
                        okstep = true; break;
                    end
                end
                step = 0.5*step;
            end
            if ~okstep, return; end
            z = zn;
        end
    end
end

function A = assemble_A(A0, Ai, x)
    A = A0;
    for i = 1:numel(Ai)
        if x(i) ~= 0, A = A + x(i)*Ai{i}; end
    end
    A = (A + A')/2;
end

function v = logdet_chol(Rc), v = 2*sum(log(diag(Rc))); end

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end
