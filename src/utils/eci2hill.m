function Rho = eci2hill(X_t, X_c)
% eci2hill.m
% Converts an absolute chaser ECI state into the RELATIVE state in the
% target-centred Hill / LVLH frame.
%
% FRAME (Schaub Hill frame, matches src/dynamics/ convention):
%   x_hat = r_t/|r_t|            radial      (R-bar)
%   z_hat = h_t/|h_t|            orbit normal (cross-track)
%   y_hat = z_hat x x_hat        along-track (V-bar), right-handed
%
% ANGULAR VELOCITY:
%   omega = (r_t x v_t)/|r_t|^2  = h_vec/r_t^2
%
% Input:
%   X_t - target ECI state [rx; ry; rz; vx; vy; vz] [m, m/s]
%   X_c - chaser ECI state [rx; ry; rz; vx; vy; vz] [m, m/s]
%
% Output:
%   Rho - relative Hill state [x; y; z; x_dot; y_dot; z_dot] [m, m/s]
%
% Inverse: hill2eci.m

r_t = X_t(1:3);  v_t = X_t(4:6);
r_c = X_c(1:3);  v_c = X_c(4:6);

% Hill basis expressed in ECI
x_hat = r_t / norm(r_t);
h_vec = cross(r_t, v_t);
z_hat = h_vec / norm(h_vec);
y_hat = cross(z_hat, x_hat);

R = [x_hat.'; y_hat.'; z_hat.'];    % rotation ECI -> Hill

% Frame angular velocity (in ECI)
omega = h_vec / norm(r_t)^2;

% Relative state via the transport theorem:
%   rho_dot_rot = rho_dot_inertial - omega x rho     (all in ECI, then rotate)
dr = r_c - r_t;
dv = v_c - v_t;

rho     = R * dr;
rho_dot = R * (dv - cross(omega, dr));

Rho = [rho; rho_dot];

end
