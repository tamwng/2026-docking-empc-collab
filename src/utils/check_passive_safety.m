function is_safe = check_passive_safety(X_log, n_logged, Ad_pow_safe, r_KOS, N_safe)
% check_passive_safety  Post-hoc passive safety check.
%
% For each state X_log(:,k), free-drifts the state forward N_safe steps
% using precomputed CWH STM powers and flags it safe iff the position never
% enters the keep-out sphere of radius r_KOS.
%
% Inputs:
%   X_log      [6 x n_logged]  logged state trajectory (Hill/LVLH frame)
%   n_logged   scalar          number of states to check
%   Ad_pow_safe[6 x 6 x N_safe] Ad^1 ... Ad^N_safe
%   r_KOS      scalar [m]      keep-out sphere radius
%   N_safe     scalar          free-drift check horizon (steps)
%
% Output:
%   is_safe    [1 x n_logged] logical

is_safe = false(1, n_logged);
for k = 1:n_logged
    xk   = X_log(:, k);
    safe = true;
    for j = 1:N_safe
        x_free = Ad_pow_safe(:,:,j) * xk;
        if norm(x_free(1:3)) < r_KOS
            safe = false;
            break;
        end
    end
    is_safe(k) = safe;
end
end
