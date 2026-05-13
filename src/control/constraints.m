function [A_ineq, b_ineq, lb, ub] = constraints(x0, Sx, Su, N, n_u, con)
% constraints.m
% Builds QP constraint matrices for MPC from a constraint config struct.
%
% Input box constraints are returned as lb/ub (quadprog bound vectors).
% State constraints are returned as A_ineq*U <= b_ineq inequality matrices.
%
% Input:
%   x0   - current state [n_x x 1]
%   Sx   - prediction matrix mapping x0 to free response [n_x*N x n_x]
%   Su   - prediction matrix mapping U to forced response [n_x*N x n_u*N]
%   N    - prediction horizon
%   n_u  - number of inputs
%   con  - constraint config struct with fields:
%            .u_max          [required] symmetric input bound [m/s^2]
%            .y_min_active   [optional] enforce y >= 0 over horizon (bool)
%
% Output:
%   A_ineq - inequality LHS [n_con x n_u*N]  ([] if no state constraints)
%   b_ineq - inequality RHS [n_con x 1]       ([] if no state constraints)
%   lb     - input lower bound [n_u*N x 1]
%   ub     - input upper bound [n_u*N x 1]

%% Input box constraints
lb = -con.u_max * ones(n_u * N, 1);
ub =  con.u_max * ones(n_u * N, 1);

%% State constraints  A_ineq * U <= b_ineq
A_list = {};
b_list = {};

% Half-space: y >= 0 (chaser stays behind target on V-bar)
% Selects along-track position (state index 2) from each horizon step.
% C_y*(Sx*x0 + Su*U) >= 0  =>  -C_y*Su*U <= C_y*Sx*x0
if isfield(con, 'y_min_active') && con.y_min_active
    C_y = kron(eye(N), [0 1 0 0 0 0]);   % [N x n_x*N]
    A_list{end+1} = -C_y * Su;
    b_list{end+1} =  C_y * Sx * x0;
end

%% Assemble
if isempty(A_list)
    A_ineq = [];
    b_ineq = [];
else
    A_ineq = vertcat(A_list{:});
    b_ineq = vertcat(b_list{:});
end

end
