function [u_opt, U_opt, X_pred] = mpc_tracking(x0, R_preview, Ad, Bd, Q, R, P, N, con)
% mpc_tracking.m
% Solves finite-horizon tracking MPC with input constraints at each time step.
%
%   min   Σᵢ₌₁ᴺ⁻¹ eᵢ'Qeᵢ + eₙ'Peₙ + Σᵢ₌₀ᴺ⁻¹ (uᵢ-u_ss)'R(uᵢ-u_ss)
%    U
%   s.t.  x(i+1) = Ad*x(i) + Bd*u(i)    (dynamics)
%         -u_max ≤ u(i) ≤ u_max          (thrust limits)
%         x(0)   = x0                    (initial condition)
%
% where eᵢ = xᵢ - r(k+i) uses the full N-step reference preview, and u_ss
% is the steady-state feedforward for the terminal reference.
%
% Input:
%   x0         - current state [6x1]
%   R_preview  - reference trajectory over horizon [6 x N]
%   Ad         - discrete state matrix [6x6]
%   Bd         - discrete input matrix [6x3]
%   Q          - state cost matrix [6x6]
%   R          - control cost matrix [3x3]
%   P          - terminal cost matrix [6x6]
%   N          - prediction horizon
%   con        - constraint config struct (see constraints.m)
%
% Output:
%   u_opt  - optimal first control input [3x1]
%   U_opt  - full optimal input sequence [3N x 1]
%   X_pred - predicted state trajectory [6 x N+1]

n_x = size(Ad, 1);   % 6 states
n_u = size(Bd, 2);   % 3 inputs

%% Stack future references into a single column vector [6N x 1]
R_vec = R_preview(:);

%% Prediction matrices Sx and Su
Sx = zeros(n_x * N, n_x);
Su = zeros(n_x * N, n_u * N);

Ad_pow = eye(n_x);
for i = 1:N
    Ad_pow = Ad_pow * Ad;
    Sx((i-1)*n_x+1 : i*n_x, :) = Ad_pow;
    for j = 1:i
        Su((i-1)*n_x+1 : i*n_x, (j-1)*n_u+1 : j*n_u) = ...
            Ad^(i-j) * Bd;
    end
end

%% Build block diagonal cost matrices
Qbar = blkdiag(kron(eye(N-1), Q), P);
Rbar = kron(eye(N), R);

%% Steady-state feedforward based on terminal reference
% u_ss satisfies r_terminal = Ad*r_terminal + Bd*u_ss in discrete time.
% Valid for fixed equilibrium setpoints (hold points).
% For zero-input references (e.g. NMC orbit), pass con.u_ss = zeros(n_u,1)
% to bypass this computation and avoid a spurious feedforward bias.
if isfield(con, 'u_ss')
    u_ss = con.u_ss;
else
    r_terminal = R_preview(:, end);
    u_ss = pinv(Bd) * (eye(n_x) - Ad) * r_terminal;
end
U_ss = repmat(u_ss, N, 1);

%% QP matrices
H = 2 * (Su' * Qbar * Su + Rbar);
f = 2 * (Su' * Qbar * (Sx * x0 - R_vec) - Rbar * U_ss);
H = (H + H') / 2;

%% Constraints
[A_ineq, b_ineq, lb, ub] = constraints(x0, Sx, Su, N, n_u, con);

%% Solve QP
opts = optimoptions('quadprog', 'Display', 'off');
[U_opt, ~, exitflag] = quadprog(H, f, A_ineq, b_ineq, [], [], lb, ub, [], opts);

if exitflag ~= 1
    warning('MPC QP did not solve to optimality. exitflag = %d', exitflag);
end

%% First control input
u_opt = U_opt(1:n_u);

%% Predicted trajectory
X_pred = zeros(n_x, N+1);
X_pred(:, 1) = x0;
for i = 1:N
    u_i = U_opt((i-1)*n_u+1 : i*n_u);
    X_pred(:, i+1) = Ad * X_pred(:, i) + Bd * u_i;
end

end
