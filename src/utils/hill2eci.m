function X_c = hill2eci(X_t, Rho)
% hill2eci.m
% Inverse of eci2hill.m: reconstructs the absolute chaser ECI state from the
% target ECI state and the relative Hill/LVLH state.
%
% Input:
%   X_t - target ECI state       [rx; ry; rz; vx; vy; vz] [m, m/s]
%   Rho - relative Hill state    [x; y; z; x_dot; y_dot; z_dot] [m, m/s]
%
% Output:
%   X_c - chaser ECI state       [rx; ry; rz; vx; vy; vz] [m, m/s]

r_t = X_t(1:3);  v_t = X_t(4:6);

% Hill basis expressed in ECI (same as eci2hill.m)
x_hat = r_t / norm(r_t);
h_vec = cross(r_t, v_t);
z_hat = h_vec / norm(h_vec);
y_hat = cross(z_hat, x_hat);

R = [x_hat.'; y_hat.'; z_hat.'];    % ECI -> Hill;  R.' is Hill -> ECI

omega = h_vec / norm(r_t)^2;        % frame angular velocity (in ECI)

rho     = Rho(1:3);
rho_dot = Rho(4:6);

% Invert the transport theorem (see eci2hill.m):
%   rho_dot = R*(dv - omega x dr)  =>  dv = R.'*rho_dot + omega x dr
dr = R.' * rho;
dv = R.' * rho_dot + cross(omega, dr);

X_c = [r_t + dr; v_t + dv];

end
