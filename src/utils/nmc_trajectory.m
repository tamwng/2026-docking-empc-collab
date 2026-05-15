function X_nmc = nmc_trajectory(rho, phi_0, n, dt, n_steps, y_c)
% nmc_trajectory.m
% Closed-form Natural Motion Circumnavigation (NMC) orbit in the Hill frame.
%
% The NMC is the drift-free solution to the unforced CWH equations obtained
% by setting ẋ₀ = 0 and ẏ₀ = -2n·x₀.  The result is a 2:1 ellipse:
%
%   x(t)  =  rho · cos(n·t + phi)
%   y(t)  =  y_c − 2·rho · sin(n·t + phi)
%   ẋ(t)  = −rho·n · sin(n·t + phi)
%   ẏ(t)  = −2·rho·n · cos(n·t + phi)
%
% Semi-axes: radial (x) = rho  [semi-minor],  along-track (y) = 2·rho  [semi-major].
% Period: T = 2π/n  (= orbital period).
%
% Input:
%   rho     - radial semi-axis [m]
%   phi_0   - phase at step k=0 [rad]
%   n       - mean motion [rad/s]
%   dt      - sample time [s]
%   n_steps - number of steps to generate
%   y_c     - along-track drift centre [m]  (default 0)
%
% Output:
%   X_nmc   - NMC state trajectory [6 x n_steps]

if nargin < 6
    y_c = 0;
end

phi = phi_0 + n * (0:n_steps-1) * dt;   % phase at each step

X_nmc = [
     rho   *  cos(phi);
     y_c   -  2*rho * sin(phi);
     zeros(1, n_steps);
    -rho*n *  sin(phi);
    -2*rho*n * cos(phi);
     zeros(1, n_steps)
];

end
