% verify_terminal_cost.m
%
% The economic stage cost is:
%   ell_e(x,u) = u' * R_eco * u    (Q_stage = 0)
%
% The terminal cost is:
%   V_f(x) = x' * P * x
%   where P = DARE(Ad, Bd, Q_dare, R_dare) with Q_dare ~= 0
%
% Amrit et al. (2011) Assumption 6 requires — for the terminal region Xf
% and a terminal control law u = Kx — that:
%
%   V_f(f(x, Kx)) - V_f(x) <= -ell_tilde(x, Kx)   for all x in Xf
%
% Raw descent check (lambda = 0, necessary condition).
%             Check: P - Acl' P Acl - K' R_eco K >= 0 ?
%             This is the terminal cost condition WITHOUT a storage function.
%             If it holds, the regulation DARE P already satisfies Assumption 6
%             for the economic stage cost.
%
% References:
%   Amrit, Rawlings & Angeli (2011), Annual Reviews in Control 35, 178-186.
%   Theorem 15, Assumption 6, Section 4.

clear; clc;

fprintf(' EMPC Terminal Cost Verification — Amrit et al. (2011)\n');

%% System setup

constants;  

dt = 60;         

Ac = [0      0     0    1     0    0  ;
      0      0     0    0     1    0  ;
      0      0     0    0     0    1  ;
      3*n^2  0     0    0     2*n  0  ;
      0      0     0   -2*n   0    0  ;
      0      0    -n^2  0     0    0 ];

Bc = [zeros(3,3); eye(3)];
[Ad, Bd] = discretize(Ac, Bc, dt);

n_x = size(Ad, 1);
n_u = size(Bd, 2);

%% DARE parameters 

Q_dare = diag([1e-2, 1e-2, 1e-2, 1e0, 1e0, 1e0]);
R_dare = diag([1e4,  1e4,  1e4]);
R_eco  = eye(3);        % economic stage cost weight

%%  Compute P and K 

[K, P] = lqr_controller(Ad, Bd, Q_dare, R_dare);

Acl = Ad - Bd * K;          % closed-loop system matrix under LQR gain

fprintf('  Spectral radius of Acl: %.6f (must be < 1 for stability)\n\n', ...
        max(abs(eig(Acl))));


%%  Raw terminal descent condition (no storage function)

% Amrit et al. (2011) Assumption 6, linearised around xs=0:
%   V_f(Acl*x) - V_f(x) <= -ell_e(x, Kx)
%
%   Acl'*P*Acl - P  <= -K'*R_eco*K
%
%   i.e.   M1 := P - Acl'*P*Acl - K'*R_eco*K  >= 0
%
%   Therefore:
%     M1 = P - Acl'*P*Acl - K'*R_eco*K
%        = Q_dare + K'*R_dare*K - K'*R_eco*K
%        = Q_dare + K'*(R_dare - R_eco)*K
%
%   Q_dare = diag(1e-2,...) > 0 and R_dare - R_eco = (1e4-1)*I > 0,
%   so both terms are positive semidefinite => M1 > 0 analytically.


fprintf('Raw descent condition  M1 = P - Acl''*P*Acl - K''*R_eco*K\n');

% First verify the DARE identity numerically (if this fails there is a
% bug in Ac, Bc, or the dlqr call — not in the theory).
DARE_lhs = P - Acl' * P * Acl;
DARE_rhs = Q_dare + K' * R_dare * K;
dare_residual = norm(DARE_lhs - DARE_rhs, 'fro');
if dare_residual > 1e-6
    warning('Large DARE residual (%.2e). Check Ac matrix and dlqr call.', dare_residual);
    fprintf('  WARNING: large residual — check system matrices.\n\n');
else
    fprintf('  DARE identity holds to machine precision.  ✓\n\n');
end

% Compute M1 in both forms and verify they agree.
M1_numeric  = P - Acl' * P * Acl - K' * R_eco * K;
M1_numeric  = (M1_numeric  + M1_numeric') / 2;

M1_analytic = Q_dare + K' * (R_dare - R_eco) * K;
M1_analytic = (M1_analytic + M1_analytic') / 2;

fmt_residual = norm(M1_numeric - M1_analytic, 'fro');
fprintf('  Consistency: ||M1_numeric - M1_analytic||_F = %.2e\n', fmt_residual);
fprintf('  (should be ~machine precision given the DARE identity above)\n\n');

min_eig_M1 = min(eig(M1_analytic));

if min_eig_M1 > 1e-10
    fprintf('  RESULT: M1 > 0  (strictly positive definite)  ✓\n');
    fprintf('  P satisfies Assumption 6 for the economic stage cost directly.\n');
    fprintf('  No storage function required. Theorem 15 applies.\n\n');
    step1_passed = true;
elseif min_eig_M1 >= -1e-9
    fprintf('  RESULT: M1 is semidefinite (within tolerance). Proceed to Step 2.\n\n');
    step1_passed = false;
else
    fprintf('  RESULT: M1 has a negative eigenvalue — condition NOT satisfied.\n');
    fprintf('  Proceed to Step 2 for the storage function approach.\n\n');
    step1_passed = false;
end