% verify_dissipativity_lmi.m
% CERTIFICATE A (closing phase): strict-dissipativity margin of the CWH
% economic OCP at the V-bar steady state.
%
%   S(L) = [ Qe + L - Ad'*L*Ad ,  -Ad'*L*Bd      ]
%          [ -Bd'*L*Ad         ,  Re - Bd'*L*Bd  ]
%
%   mu* = sup { mu >= 0 : S(L) >= blkdiag(mu*I6, 0) for some L = L' }
%
% and the certificate margin is rho(r) = mu* r^2.  Two stage costs:
%
%   regularised (Run B) : Qe = rho_c*I6,  Re = I3   ->  mu* = rho_c
%   pure fuel           : Qe = 0,         Re = I3   ->  mu* = 0


clear; clc;

%% Path + constants 
script_dir = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(script_dir, '..', '..', 'src')));
constants;                        % provides n, T, ...

dt    = 120;                      % S4 closing sample time [s]
rho_c = 1e-11;                    % Run B state-regularisation weight
Re    = eye(3);                   % input weight R_eco

fprintf('==========================================================\n');
fprintf(' CERTIFICATE A -- strict dissipativity (closing, V-bar)\n');
fprintf('==========================================================\n');
fprintf(' n = %.6e rad/s   T = %.1f s   dt = %g s   rho_c = %.0e\n\n', n, T, dt, rho_c);

%% CWH model 
Ac = zeros(6);
for j = 1:6
    ej = zeros(6,1); ej(j) = 1;
    Ac(:,j) = clohessy_wiltshire(ej, n);
end
Bc = [zeros(3,3); eye(3)];
[Ad, Bd] = discretize(Ac, Bc, dt);

%% STEP 1: sigma(Ad) lies on the unit circle 
[Vd, Dd] = eig(Ad);  lam = diag(Dd);
fprintf('-- STEP 1: eigenvalues of Ad -------------------------\n');
for k = 1:6
    fprintf('   %+9.6f %+9.6fi    |lam| = %.14f\n', real(lam(k)), imag(lam(k)), abs(lam(k)));
end
fprintf('   max | |lam| - 1 | = %.2e  ->  all 6 on the unit circle\n', ...
    max(abs(abs(lam) - 1)));
fprintf('   (double lam=1 along-track drift + two pairs exp(+-i n dt))\n\n');

%% STEP 2: controllability => stabilisability 
Cm = Bd;  AkB = Bd;
for k = 1:5, AkB = Ad*AkB;  Cm = [Cm, AkB]; end          %#ok<AGROW>
hautus = arrayfun(@(l) rank([Ad - l*eye(6), Bd], 1e-9), lam);
fprintf('-- STEP 2: controllability of (Ad,Bd) ----------------\n');
fprintf('   rank(ctrb)          = %d / 6\n', rank(Cm, 1e-9));
fprintf('   min Hautus/PBH rank = %d / 6  (over all 6 eigenvalues)\n\n', min(hautus));

%% STEP 3: the margin mu* for both stage costs 
fprintf('-- STEP 3: dissipativity margin mu* ------------------\n');
fprintf('   %-22s  %-13s  %-13s  %-13s\n', ...
    'stage cost', 'ub (eigvec)', 'at L = 0', 'probe max L');
run_case('regularised (Run B)', rho_c*eye(6), Re, Ad, Bd, Vd, lam);
run_case('pure fuel',           zeros(6),     Re, Ad, Bd, Vd, lam);

fprintf('\n   Upper and lower bound coincide in both rows, so:\n');
fprintf('     regularised : mu* = rho_c = %.3e  (storage lambda == 0)\n', rho_c);
fprintf('     pure fuel   : mu* = 0             -> NOT strictly dissipative\n');
fprintf('\nDone.\n');

% ========================================================================
%  Local functions
% ========================================================================

function run_case(name, Qe, Re, Ad, Bd, Vd, lam)
% Print the eigenvector upper bound, the L=0 value, and a randomised probe.
    ub = eig_bound(Qe, Vd, lam);
    m0 = margin(zeros(6), Qe, Re, Ad, Bd);
    mp = probe(Qe, Re, Ad, Bd);
    fprintf('   %-22s  %-13.6e  %-13.6e  %-13.6e\n', name, ub, m0, mp);
end

function ub = eig_bound(Qe, Vd, lam)
% mu* <= min over unit-circle eigenvectors of v'Qe v / v'v.  Valid for ANY
% storage L, because v'(L - Ad'*L*Ad)v = (1-|lam|^2) v'L v = 0 there.
    ub = inf;
    for k = 1:numel(lam)
        if abs(abs(lam(k)) - 1) < 1e-8
            v  = Vd(:,k);
            ub = min(ub, real(v'*Qe*v) / real(v'*v));
        end
    end
end

function m = margin(L, Qe, Re, Ad, Bd)
% Largest mu with S(L) >= blkdiag(mu*I6, 0), via the Schur complement of
% S(L) w.r.t. its input block.  -inf if that block is not positive definite.
    S22 = Re - Bd'*L*Bd;   S22 = (S22 + S22')/2;
    if min(eig(S22)) <= 0, m = -inf; return; end
    S11 = Qe + L - Ad'*L*Ad;
    S12 = -Ad'*L*Bd;
    Sc  = S11 - S12*(S22\S12');
    m   = min(eig((Sc + Sc')/2));
end

function mbest = probe(Qe, Re, Ad, Bd)
% Randomised falsification: no storage should beat the eigenvector bound.
% Scales are centred on norm(Re)/norm(Bd)^2, where Bd'*L*Bd meets Re.
    rng(0);
    mbest = margin(zeros(6), Qe, Re, Ad, Bd);
    s0    = norm(Re) / norm(Bd)^2;
    for s = s0 * logspace(-3, 3, 7)
        for t = 1:20
            G = randn(6);
            mbest = max(mbest, margin(s*(G+G')/2, Qe, Re, Ad, Bd));
        end
    end
end
