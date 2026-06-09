function Phi = cw_stm(n_orb, tau)
% cw_stm  Clohessy-Wiltshire state transition matrix (closed form).
%   Phi = cw_stm(n_orb, tau) returns the 6x6 CW STM for mean motion
%   n_orb [rad/s] and free-drift time tau [s].
%   State order: [x; y; z; vx; vy; vz] (radial, along-track, cross-track).
%
%   Reference: Schaub & Junkins, "Analytical Mechanics of Space Systems",
%   Chapter 14 (Hill / CW equations, Eq. 14.19).
c  = cos(n_orb * tau);
s  = sin(n_orb * tau);
nt = n_orb * tau;
Phi = [4 - 3*c,           0,  0,   s / n_orb,         2*(1-c) / n_orb,     0       ;
       6*(s - nt),        1,  0,  -2*(1-c) / n_orb,   (4*s - 3*nt) / n_orb, 0       ;
       0,                 0,  c,   0,                  0,                   s / n_orb;
       3*n_orb*s,         0,  0,   c,                  2*s,                 0       ;
      -6*n_orb*(1 - c),   0,  0,  -2*s,                4*c - 3,             0       ;
       0,                 0, -n_orb*s, 0,               0,                   c       ];
end
