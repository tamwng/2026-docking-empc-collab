function [A_ineq, b_ineq, lb, ub, n_ps_added] = constraints(x0, Sx, Su, N, n_u, con)
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
%            .approach_corridor.active [optional] rectangular tube along V-bar (bool)
%            .approach_corridor.d_max  [optional] half-width [m] (|x|,|z| <= d_max)
%
% Output:
%   A_ineq - inequality LHS [n_con x n_u*N]  ([] if no state constraints)
%   b_ineq - inequality RHS [n_con x 1]       ([] if no state constraints)
%   lb     - input lower bound [n_u*N x 1]
%   ub     - input upper bound [n_u*N x 1]

n_x        = size(Sx, 2);    % state dimension
n_ps_added = 0;              % passive safety rows added this call

%% Input box constraints
lb = -con.u_max * ones(n_u * N, 1);
ub =  con.u_max * ones(n_u * N, 1);

%% State constraints  A_ineq * U <= b_ineq
A_list = {};
b_list = {};

%% Half-space: y >= 0 (chaser stays behind target on V-bar)
% Selects along-track position (state index 2) from each horizon step.
% C_y*(Sx*x0 + Su*U) >= 0  =>  -C_y*Su*U <= C_y*Sx*x0
if isfield(con, 'y_min_active') && con.y_min_active
    C_y = kron(eye(N), [0 1 0 0 0 0]);   % [N x n_x*N]
    A_list{end+1} = -C_y * Su;
    b_list{end+1} =  C_y * Sx * x0;
end

%% LoS cone: inner polyhedral approximation of the V-bar approach corridor.
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

%% Approach corridor: rectangular tube along the V-bar axis
% Box constraint on radial (x) and cross-track (z) at every horizon step;
% y (along-track) is unconstrained.  Adds 4N rows (±x, ±z at each step).
if isfield(con, 'approach_corridor') && con.approach_corridor.active
    d_max  = con.approach_corridor.d_max;
    C_x    = kron(eye(N), [1,0,0,0,0,0]);   % [N × n_x*N]
    C_z    = kron(eye(N), [0,0,1,0,0,0]);   % [N × n_x*N]
    free_x = C_x * Sx * x0;                 % free-response x component [N×1]
    free_z = C_z * Sx * x0;                 % free-response z component [N×1]
    A_list{end+1} = [ C_x * Su; -C_x * Su; C_z * Su; -C_z * Su];
    b_list{end+1} = [ d_max*ones(N,1) - free_x;
                      d_max*ones(N,1) + free_x;
                      d_max*ones(N,1) - free_z;
                      d_max*ones(N,1) + free_z];
end

%% Terminal ball: position-only L-inf box  -r <= e_pos(N) <= r
% e_pos(N) = C_pos3*(Sx_N*x0 + Su_N*U)  [3×1, position part of terminal error]
% Velocity rows removed: with V_f=0 (no terminal cost) there is no stabilising
% mechanism for terminal velocity, so boxing all 6 states over-constrains the QP
% without providing any additional passive-safety guarantee beyond position
% containment.  Exact L2 terminal ball would require SOCP; this L-inf box is
% a conservative inner approximation on position only.
if isfield(con, 'terminal_ball') && con.terminal_ball.active
    r_sw    = con.terminal_ball.r_switch;
    n_x_l   = size(Sx, 2);                     % = n_x = 6
    C_pos3  = [eye(3), zeros(3,3)];            % 3×6 position selector
    Su_N    = Su(end-n_x_l+1:end, :);          % 6 × n_u*N  (full terminal block)
    Sx_N    = Sx(end-n_x_l+1:end, :);          % 6 × n_x
    Su_Np   = C_pos3 * Su_N;                   % 3 × n_u*N
    free_Np = C_pos3 * Sx_N * x0;             % 3 × 1
    A_list{end+1} = [ Su_Np; -Su_Np];
    b_list{end+1} = [ r_sw*ones(3,1) - free_Np;
                      r_sw*ones(3,1) + free_Np];
end

%% Passive safety: SCA-linearised KOS half-space per prediction step ─────────
%
% For each prediction step k, free-drift the reference predicted actual state
% forward over N_safe steps and find the worst-case approach to the KOS.
% A single linearised half-space n'*C*Phi*x_k >= r_KOS is added only when the
% reference comes within (r_KOS + safety_margin) of the KOS boundary.
%
% con.passive_safety fields:
%   .active         — bool
%   .Phi_array      — [6×6×N_safe] precomputed CW STMs
%   .C_pos          — [3×6] position extractor [eye(3), zeros(3,3)]
%   .r_KOS          — keep-out sphere radius [m]
%   .safety_margin  — activation margin beyond r_KOS [m]
%   .x_switch       — [6×1] switch point (actual state offset for error coords)
%   .U_warm         — [3N×1] SCA warm-start input sequence
if isfield(con, 'passive_safety') && con.passive_safety.active
    ps     = con.passive_safety;
    C_p    = ps.C_pos;
    r_k    = ps.r_KOS;
    s_mg   = ps.safety_margin;
    Ph_arr = ps.Phi_array;
    N_s    = size(Ph_arr, 3);
    U_w    = ps.U_warm;
    x_sw   = ps.x_switch;

    % Reference predicted states in error coords, then shift to actual coords
    x_pred_ref_e = reshape(Sx * x0 + Su * U_w, n_x, N);   % n_x × N (error)

    for k = 1:N
        x_k_ref = x_pred_ref_e(:, k) + x_sw;   % actual coordinates

        % Find worst-case free-drift time (exhaustive search over all steps)
        min_d  = inf;
        best_i = 1;
        for i = 1:N_s
            p_i = C_p * Ph_arr(:,:,i) * x_k_ref;
            d_i = norm(p_i);
            if d_i < min_d
                min_d  = d_i;
                best_i = i;
            end
        end

        if min_d < r_k + s_mg
            Phi_star = Ph_arr(:,:,best_i);
            p_ref    = C_p * Phi_star * x_k_ref;
            if norm(p_ref) < 1e-10; continue; end   % degenerate — skip
            n_vec = p_ref / norm(p_ref);

            Su_k  = Su((k-1)*n_x+1 : k*n_x, :);
            Sx_k  = Sx((k-1)*n_x+1 : k*n_x, :);
            % n_vec' * C_p * Phi_star * (Sx_k*x0 + Su_k*U + x_sw) >= r_k
            A_row = -(n_vec' * C_p * Phi_star * Su_k);
            b_row = -(r_k - n_vec' * C_p * Phi_star * (Sx_k * x0 + x_sw));
            A_list{end+1} = A_row;  %#ok<AGROW>
            b_list{end+1} = b_row;  %#ok<AGROW>
            n_ps_added = n_ps_added + 1;
        end
    end
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
