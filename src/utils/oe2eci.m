function X = oe2eci(a, e, i, RAAN, argp, nu, mu)
% oe2eci.m
% Classical (Keplerian) orbital elements -> ECI state vector.
% Used to seed target initial conditions in the model-validity sweeps
% directly from the swept parameters (a, e, i).
%
% Input (angles in RADIANS):
%   a    - semi-major axis [m]
%   e    - eccentricity [-]
%   i    - inclination [rad]
%   RAAN - right ascension of ascending node, Omega [rad]
%   argp - argument of perigee, omega [rad]
%   nu   - true anomaly [rad]
%   mu   - gravitational parameter [m^3/s^2]
%
% Output:
%   X    - ECI state [rx; ry; rz; vx; vy; vz] [m, m/s]
%
% Reference: Vallado, "Fundamentals of Astrodynamics and Applications",
% Algorithm 10 (COE2RV).

p = a * (1 - e^2);                 % semi-latus rectum [m]
r = p / (1 + e*cos(nu));           % orbital radius [m]

% State in perifocal (PQW) frame
r_pf = [ r*cos(nu);  r*sin(nu);  0 ];
v_pf = sqrt(mu/p) * [ -sin(nu);  e + cos(nu);  0 ];

% Perifocal -> ECI rotation (3-1-3: RAAN, incl, argp)
cO = cos(RAAN); sO = sin(RAAN);
ci = cos(i);    si = sin(i);
cw = cos(argp); sw = sin(argp);

Q = [ cO*cw - sO*sw*ci,  -cO*sw - sO*cw*ci,   sO*si;
      sO*cw + cO*sw*ci,  -sO*sw + cO*cw*ci,  -cO*si;
      sw*si,              cw*si,              ci    ];

X = [ Q * r_pf;  Q * v_pf ];

end
