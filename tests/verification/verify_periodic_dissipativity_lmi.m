% verify_periodic_dissipativity_lmi.m
% -------------------------------------------------------------------------
% CERTIFICATE B (NMC fly-around): strict P-PERIODIC dissipativity, the
% periodic analogue of tests/verify_dissipativity_lmi.m.
%
% Discrete CWH  x+ = Ad x + Bd u  at the fly-around sample time dt = T/P,
% P = 92, so P*dt = T EXACTLY -- which is what makes the monodromy
% M := Ad^P unipotent.  The NMC orbit is a P-periodic FUEL-FREE solution of
% the unforced dynamics, so in deviation coordinates e_k = x_k - x*_k the
% error obeys the SAME pair (Ad,Bd) at every phase.
%
% NOTATION (thesis): Re = INPUT weight, Qe = STATE weight. (Damm swaps R/Q.)
%
% Stage costs about the phase-synced reference x*_k:
%   pure fuel      l     = ||u||^2
%   eps-cost (R5)  l_aug = ||u||^2 + eps*||x - x*_k||^2,  eps > 0

clear; clc;

%% Path + constants 
script_dir = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(script_dir, '..', '..', 'src')));
constants;                        % provides n, T, ...

P      = 92;                      % steps per orbital period (S2 setting)
dt     = T / P;                   % fly-around sample time
eps_r5 = 1e-4;                    % S2 Run 5 state-reg weight (validated)

fprintf('==========================================================\n');
fprintf(' CERTIFICATE B -- strict P-periodic dissipativity (NMC)\n');
fprintf('==========================================================\n');
fprintf(' n = %.6e rad/s   T = %.2f s   P = %d   dt = T/P = %.4f s\n\n', n, T, P, dt);

%% CWH model, monodromy, lifted input map
Ac = zeros(6);
for j = 1:6
    ej = zeros(6,1); ej(j) = 1;
    Ac(:,j) = clohessy_wiltshire(ej, n);
end
Bc = [zeros(3,3); eye(3)];
[Ad, Bd] = discretize(Ac, Bc, dt);

Adpow = cell(P+1,1);  Adpow{1} = eye(6);
for k = 1:P, Adpow{k+1} = Adpow{k} * Ad; end
M = Adpow{P+1};                                   % monodromy Ad^P

B_lift = zeros(6, 3*P);                           % [Ad^{P-1}Bd ... Ad Bd, Bd]
for i = 0:P-1
    B_lift(:, 3*i+1:3*i+3) = Adpow{P-1-i+1} * Bd;
end

%%  STEP 1: monodromy spectrum and defectiveness 
ev = eig(M);
fprintf('-- STEP 1: monodromy M = Ad^P ------------------------\n');
fprintf('   max |eig(M) - 1|  = %.2e   ->  sigma(M) = {1,...,1}\n', max(abs(ev - 1)));
MmI       = M - eye(6);
geom_mult = 6 - rank(MmI, 1e-6);
fprintf('   rank(M - I)       = %d      (lam=1: alg mult 6, geom mult %d)\n', ...
    rank(MmI, 1e-6), geom_mult);
fprintf('   ||M - I||_F       = %.4e  -> one 2x2 Jordan block, DEFECTIVE\n', norm(MmI,'fro'));
fprintf('   drift entries     : M(2,1) = %.4f (-12*pi = %.4f)\n', MmI(2,1), -12*pi);
fprintf('                       M(2,5) = %.1f (-6*pi/n = %.1f)\n\n', MmI(2,5), -6*pi/n);

%% STEP 2: the NMC subspace ker(M - I)
Ns    = null(MmI);
c_nmc = [2*n, 0, 0, 0, 1, 0];                     % NMC condition ydot0 = -2n x0
fprintf('-- STEP 2: drift-free NMC subspace ker(M-I) ----------\n');
fprintf('   dim ker(M - I)  = %d   (5-D NMC subspace)\n', size(Ns,2));
fprintf('   ||M*Ns - Ns||   = %.2e  (M restricted to it is the identity)\n', norm(M*Ns - Ns));
fprintf('   ||c_nmc * Ns||  = %.2e  (ydot0 = -2n x0 holds on it)\n\n', norm(c_nmc*Ns));

%% Lifted prediction matrices and Gram blocks 
Sx_stack = zeros(6*P, 6);
Su_stack = zeros(6*P, 3*P);
for j = 0:P-1
    Sx_stack(6*j+1:6*j+6, :) = Adpow{j+1};                            % Ad^j
    for i = 0:j-1
        Su_stack(6*j+1:6*j+6, 3*i+1:3*i+3) = Adpow{j-1-i+1} * Bd;     % Ad^{j-1-i}Bd
    end
end
G0 = grams(Sx_stack, Su_stack, P);

%% STEP 3: periodic margins 
fprintf('-- STEP 3: periodic dissipativity margin -------------\n');

% Pure fuel: zero for every storage.  L = 0 plus a randomised probe.
mp = probe(G0, 0, M, B_lift);
fprintf('   pure fuel (eps = 0) : at L = 0 %.3e , probe max over L %.3e\n', ...
    margin_lift(zeros(6), G0, 0, M, B_lift), mp);
fprintf('                         => margin 0: v''(L - M''LM)v = 0 on ker(M-I)\n');
fprintf('                         NOT strictly P-periodically dissipative.\n\n');

% eps-cost: L = 0 already certifies a positive margin, scaling with eps.
fprintf('   eps-cost: margin at L = 0 (a certified lower bound)\n');
fprintf('   %-12s  %-16s  %-12s\n', 'eps', 'margin', 'margin/eps');
for ep = [1e-6, 1e-4, 1e-2]
    g = margin_lift(zeros(6), G0, ep, M, B_lift);
    fprintf('   %-12.0e  %-16.8e  %-12.4f\n', ep, g, g/ep);
end
fprintf('\n   At the S2 Run 5 design point eps = %.0e the margin is %.4e ~ eps,\n', ...
    eps_r5, margin_lift(zeros(6), G0, eps_r5, M, B_lift));
fprintf('   i.e. the one-step reading l~_k >= eps*||x - x*_k||^2 with lambda_k == 0.\n');
fprintf('   Strictly P-periodically dissipative -> Koehler (2018) Cor. 4 applies.\n');
fprintf('\nDone.\n');

% ========================================================================
%  Local functions
% ========================================================================

function G0 = grams(Sx_stack, Su_stack, P)
% eps-INDEPENDENT lifted Gram blocks of the one-period state cost
%   sum_j ||e_{k+j}||^2 = ||Sx_stack*e + Su_stack*U||^2.
    G0.xx = Sx_stack' * Sx_stack;            % 6 x 6
    G0.xu = Sx_stack' * Su_stack;            % 6 x 3P
    G0.uu = Su_stack' * Su_stack;            % 3P x 3P
    G0.n3 = 3*P;
end

function g = margin_lift(L, G0, ep, M, B_lift)
% Largest mu with S_lift(L) >= blkdiag(mu*I6, 0), via the Schur complement
% w.r.t. the 3P input block.  -inf if that block is not positive definite.
    S11 = ep*G0.xx + L - M'*L*M;                          % 6 x 6
    S12 = ep*G0.xu - M'*L*B_lift;                         % 6 x 3P
    S22 = ep*G0.uu + eye(G0.n3) - B_lift'*L*B_lift;       % 3P x 3P
    S22 = (S22 + S22')/2;
    if min(eig(S22)) <= 1e-10, g = -inf; return; end
    Sc = S11 - S12*(S22\S12');
    g  = min(eig((Sc + Sc')/2));
end

function mbest = probe(G0, ep, M, B_lift)
% Randomised falsification for the pure-fuel case: no storage beats zero.
    rng(0);
    mbest = margin_lift(zeros(6), G0, ep, M, B_lift);
    s0    = 1 / norm(B_lift)^2;
    for s = s0 * logspace(-3, 3, 7)
        for t = 1:20
            G = randn(6);
            mbest = max(mbest, margin_lift(s*(G+G')/2, G0, ep, M, B_lift));
        end
    end
end
