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
%            .u_max              [required] symmetric input bound [m/s^2]
%            .y_min_active       [optional] enforce y >= 0 over horizon (bool)
%            .los_cone.active     [optional] enforce V-bar LoS cone (bool)
%            .los_cone.half_angle [optional] cone half-angle [rad]
%            .los_cone.n_faces    [optional] faces in inner polygon (e.g. 10)
%            .los_cone.apex       [optional] cone apex [x;y;z], default [0;0;0]
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

% LoS cone: inner polyhedral approximation of the V-bar approach corridor.
%
% Physical setup: docking port at origin, docking axis along +y (chaser
% approaches from positive y). The cone is:
%   sqrt(x^2 + z^2) <= y * tan(alpha)   for y >= 0
%
% Inner M-face polyhedral approximation — face j gives one linear row:
%   cos(2*pi*j/M)*x - tan(alpha)*cos(pi/M)*y + sin(2*pi*j/M)*z <= 0
%
% The cos(pi/M) factor is the inradius/circumradius ratio of the inscribed
% polygon, making the approximation strictly inner (conservative).
%
% con.los_cone fields:
%   .active     - bool
%   .half_angle - cone half-angle [rad]
%   .n_faces    - number of polyhedral faces (10 matches Weiss 2015)
if isfield(con, 'los_cone') && con.los_cone.active
    alpha  = con.los_cone.half_angle;
    M      = con.los_cone.n_faces;
    n_x    = size(Sx, 2);

    % Cone apex = docking port position [x_a; y_a; z_a]; defaults to origin.
    if isfield(con.los_cone, 'apex')
        apex = con.los_cone.apex;
    else
        apex = zeros(3, 1);
    end

    alpha_eff = tan(alpha) * cos(pi / M);   % effective radial slope (inradius factor)

    Acone_local = zeros(M, n_x);
    bcone_local = zeros(M, 1);
    for j = 0:M-1
        theta = 2*pi*j / M;
        Acone_local(j+1, :) = [cos(theta), -alpha_eff, sin(theta), 0, 0, 0];
        bcone_local(j+1)    =  cos(theta)*apex(1) + sin(theta)*apex(3) - alpha_eff*apex(2);
    end

    Acone = kron(eye(N), Acone_local);      % [M*N x n_x*N]
    bcone = repmat(bcone_local, N, 1);      % same RHS at every horizon step

    A_list{end+1} = Acone * Su;
    b_list{end+1} = bcone - Acone * Sx * x0;
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
