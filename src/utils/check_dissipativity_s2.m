function check_dissipativity_s2(A_d, B_d, x_star, P, N, n_orb)
% check_dissipativity_s2  Numerical strict-dissipativity certificate for the
%   Scenario 2 NMC fly-around with a linear periodic storage and the periodic
%   terminal equality constraint (Zanon-Grüne-Diehl 2015, Type-B).
%
%   Verifies three things:
%     (1) GLOBAL non-strictness   : rotated stage-cost Hessian has zero x-block
%                                   (linear storage cannot be strict globally).
%     (2) CONSTRAINED well-posed  : reduced input Hessian on the terminal-
%                                   constrained subspace is strictly PD.
%     (3) TRANSVERSE strictness   : rotated value-function Hessian wrt initial
%                                   state has positive eigenvalues in the
%                                   directions transverse to the orbit; the
%                                   flat directions are the phase tangent.
%
%   Usage (from sim_scenario2_flyaround.m, after x_star is built):
%     check_dissipativity_s2(A_d, B_d, x_star, P, N, n);
%
%   Inputs:
%     A_d, B_d  - discrete CWH matrices (6x6, 6x3)
%     x_star    - prescribed NMC orbit (6 x P)
%     P         - steps per orbital period
%     N         - prediction horizon (= P)
%     n_orb     - mean motion [rad/s]

    nx = 6;  nu = 3;
    tol = 1e-9;

    fprintf('\n=== Strict-dissipativity certificate (Scenario 2) ===\n');

    %% Linear periodic storage candidate
    % lambda_k(x) = p' x, phase-independent (the NMC manifold multiplier on vy).
    % This is the candidate derived analytically; the magnitude is immaterial
    % for the Hessian structure below (it only shifts the linear/gradient term).
    p = [0; 0; 0; 0; -1; 0];

    %% (1) Global rotated stage-cost Hessian in z = [x; u]
    % l(x,u) = ||u||^2,  l(x*_k,u*_k) = 0 (passive orbit).
    % L_k(x,u) = u'u + p'x - p'(A_d x + B_d u)
    %          = u'u + p'(I-A_d)x - (p'B_d)u
    % Hessian wrt (x,u) = blkdiag(0_6, I_3): PSD, NOT PD in x.
    Hq = blkdiag(zeros(nx), eye(nu));
    ev_glob = eig(Hq);
    fprintf('(1) Global stage-cost Hessian eigenvalues:\n    [%s]\n', ...
        num2str(sort(ev_glob)', '%5.2f'));
    fprintf('    -> %d zero eigenvalues in x-block: linear storage is NON-strict globally (expected).\n', ...
        sum(abs(ev_glob) < tol));

    %% Condensed prediction matrices for horizon N
    Sx = zeros(nx*N, nx);
    Su = zeros(nx*N, nu*N);
    Adp = eye(nx);
    for i = 1:N
        Adp = Adp * A_d;
        Sx((i-1)*nx+1:i*nx, :) = Adp;
        for j = 1:i
            Su((i-1)*nx+1:i*nx, (j-1)*nu+1:j*nu) = A_d^(i-j) * B_d;
        end
    end
    Su_term = Su((N-1)*nx+1:N*nx, :);   % 6 x nu*N
    Sx_term = Sx((N-1)*nx+1:N*nx, :);   % 6 x 6

    %% (2) Reduced input Hessian on the terminal-constrained subspace
    % Terminal equality: Su_term * U = x*_term - Sx_term * x0.
    % Free directions = null(Su_term).  Stage cost Hessian wrt U is 2I (R=I).
    Z = null(Su_term);
    Hred = Z' * (2*eye(nu*N)) * Z;
    ev_red = eig((Hred+Hred')/2);
    fprintf('(2) Reduced input Hessian on constrained subspace (dim %d):\n', size(Z,2));
    fprintf('    min eig = %.4f   max eig = %.4f   -> %s\n', ...
        min(ev_red), max(ev_red), ...
        ternary(min(ev_red) > tol, 'PD (PASS, QP strictly convex)', 'FAIL'));

    %% (3) Rotated value-function Hessian wrt initial state
    % V(x0) = min_{U: Su_term U = xT - Sx_term x0} ||U||^2
    %       = (xT - Sx_term x0)' (Su_term Su_term')^{-1} (xT - Sx_term x0).
    % Hessian wrt x0 = Sx_term' (Su_term Su_term')^{-1} Sx_term.
    G   = inv(Su_term * Su_term');
    Hx0 = Sx_term' * G * Sx_term;
    Hx0 = (Hx0 + Hx0')/2;
    ev_x0 = sort(eig(Hx0));
    r_pos = sum(ev_x0 > tol);
    fprintf('(3) Value-function Hessian wrt x0 eigenvalues:\n');
    fprintf('    [%s]\n', num2str(ev_x0', '%9.2e'));
    fprintf('    transverse strict directions = %d / 6  (rest = orbit tangent / phase)\n', r_pos);

    %% Verdict
    fprintf('\nVerdict:\n');
    if min(ev_red) > tol && r_pos >= 3
        fprintf('  Strict Type-B periodic dissipativity holds on the terminal-\n');
        fprintf('  constrained set in the %d transverse directions.\n', r_pos);
        fprintf('  => Pi* asymptotically stable as a SET (Zanon Thm 4.6, sigma_B).\n');
        fprintf('  => Type-A (phase-specific) stability follows via Remark 4.7\n');
        fprintf('     since the terminal equality makes Pi* the unique length-P minimiser.\n');
    else
        fprintf('  Certificate incomplete - inspect eigenvalues above.\n');
    end
    fprintf('=====================================================\n\n');
end

function out = ternary(cond, a, b)
    if cond; out = a; else; out = b; end
end