function [u_opt, U_opt, X_pred, exitflag, n_ps_added] = mpc_regulation(x0, Ad, Bd, Q, R, P, N, con)
% mpc_regulation.m
% Solves the finite-horizon constrained tracking MPC (regulation to origin
%  x = 0) problem at each timestep.
% Equivalent to LQR with finite horizon and hard input constraints.
%
%   min   Σᵢ₌₀ᴺ⁻¹ (xᵢ'Qxᵢ + uᵢ'Ruᵢ) + xₙ'Pxₙ
%    U
%   s.t.  x(i+1) = Ad*x(i) + Bd*u(i)    (dynamics)
%         -u_max ≤ u(i) ≤ u_max          (thrust limits)
%         x(0)   = x0                    (initial condition)

%   x0    - current state [6x1]
%   Ad    - discrete state matrix [6x6]
%   Bd    - discrete input matrix [6x3]
%   Q     - state cost matrix [6x6]
%   R     - control cost matrix [3x3]
%   P     - terminal cost matrix [6x6]
%   N     - prediction horizon
%   con   - constraint config struct with fields:
%             .u_max          [required] symmetric thrust bound [m/s^2]
%             .y_min_active   [optional] enforce y >= 0 over horizon (bool)
%
% Output:
%   u_opt  - optimal first control input [3x1]
%   U_opt  - full optimal input sequence [3N x 1]
%   X_pred - predicted state trajectory [6 x N+1]
%
n_x = size(Ad, 1);
n_u = size(Bd, 2);

%% Prediciton matrices (for condensed QP, only optimize over u)
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

%% Terminal constraint selection
% terminal_eq: x(N) = 0     (regulation to origin)
% periodic:    x(N) = x(0)  (orbit closes on itself)
% default:     terminal cost P*x(N)'*x(N)
use_terminal_eq = isfield(con, 'terminal_eq') && con.terminal_eq;
use_periodic    = isfield(con, 'periodic')    && con.periodic;

if use_terminal_eq
    P_cost = zeros(n_x);
    Aeq = Su(end-n_x+1:end, :);
    if isfield(con, 'x_target')
        beq = con.x_target - Sx(end-n_x+1:end, :) * x0;
    else
        beq = -Sx(end-n_x+1:end, :) * x0;
    end
elseif use_periodic
    P_cost = zeros(n_x);
    Aeq = Su(end-n_x+1:end, :);
    beq = (eye(n_x) - Sx(end-n_x+1:end, :)) * x0;
else
    P_cost = P;
    Aeq = [];
    beq = [];
end

%% Cost matrices
% Q at intermediate stages, P_cost at terminal stage
Qbar = blkdiag(kron(eye(N-1), Q), P_cost);
Rbar = kron(eye(N), R);

%% Form QP matrices
H = 2 * (Su' * Qbar * Su + Rbar);
f = 2 * Su' * Qbar * Sx * x0;
H = (H + H') / 2;   % symmetrize for numerical stability

%% Constraints (input bounds + optional state constraints)
[A_ineq, b_ineq, lb, ub, n_ps_added] = constraints(x0, Sx, Su, N, n_u, con);

%% Solve QP
opts = optimoptions('quadprog', 'Display', 'off');
[U_opt, ~, exitflag] = quadprog(H, f, A_ineq, b_ineq, Aeq, beq, lb, ub, [], opts);

if exitflag ~= 1
    warning('MPC QP did not solve to optimality. exitflag = %d', exitflag);
end

%% first control input (1:3 vector for MPC)
u_opt = U_opt(1:n_u);

%% ── Reconstruct predicted trajectory ─────────────────────────────────────
X_pred = zeros(n_x, N+1);
X_pred(:, 1) = x0;
for i = 1:N
    u_i = U_opt((i-1)*n_u+1 : i*n_u); % 1:3, then 4:6, 7:9 ....
    X_pred(:, i+1) = Ad * X_pred(:, i) + Bd * u_i;
end

end