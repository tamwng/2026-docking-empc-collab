function dX = truth_propagator_at(t, X, p)
% truth_propagator_at.m
% High-fidelity ECI equations of motion for a SINGLE spacecraft, built on
% the Aerospace Toolbox.
%
% Input:
%   t  - time since epoch [s] (scalar; needed for Earth rotation / atmosphere)
%   X  - spacecraft ECI state [rx; ry; rz; vx; vy; vz] [m, m/s]
%   p  - parameter struct with fields:
%          .mu           gravitational parameter [m^3/s^2]
%          .epoch        UTC epoch as datetime OR [Y M D h m s] vector
%          .use_zonal    logical, enable J2-J4 (default false)
%          .zonal_degree 2, 3, or 4 (default 4)
%          .use_drag     logical, enable atmospheric drag (default false)
%          .Cd           drag coefficient [-]        (drag only)
%          .A            reference cross-section [m^2](drag only)
%          .m            spacecraft mass [kg]        (drag only)
%          .omega_earth  Earth rotation rate [rad/s] (drag only)
%          .reduction    dcmeci2ecef reduction, default 'IAU-2000/2006'
%
% Output:
%   dX - time derivative [vx; vy; vz; ax; ay; az]
%
% References:
%   Vallado, "Fundamentals of Astrodynamics and Applications".
%   Aerospace Toolbox: gravityzonal, dcmeci2ecef, ecef2lla, atmosnrlmsise00.

% defaults
if ~isfield(p, 'use_zonal'),    p.use_zonal = false;              end
if ~isfield(p, 'use_drag'),     p.use_drag  = false;              end
if ~isfield(p, 'zonal_degree'), p.zonal_degree = 4;               end
if ~isfield(p, 'reduction'),    p.reduction = 'IAU-2000/2006';    end

r = X(1:3);
v = X(4:6);
r_norm = norm(r);

% central-body gravity
if p.use_zonal
    utc     = epoch_to_utc(p.epoch, t);
    dcm     = dcmeci2ecef(p.reduction, utc);   % 3x3 ECI -> ECEF
    r_ecef  = dcm * r;                          % position in ECEF [m]

    [gx, gy, gz] = gravityzonal(r_ecef.', p.zonal_degree);  % ECEF accel [m/s^2]
    a_grav_ecef  = [gx; gy; gz];
    a_grav       = dcm.' * a_grav_ecef;         % rotate back to ECI
else
    % Pure point-mass gravity
    a_grav = -p.mu * r / r_norm^3;
end

a = a_grav;

% drag
if p.use_drag
    if ~p.use_zonal
        utc    = epoch_to_utc(p.epoch, t);
        dcm    = dcmeci2ecef(p.reduction, utc);
        r_ecef = dcm * r;
    end

    lla = ecef2lla(r_ecef.');            % [lat(deg) lon(deg) alt(m)]
    lat = lla(1); lon = lla(2); alt = lla(3);

    % Atmospheric density from NRLMSISE-00 (total mass density = column 6).
    [year, doy, ut_sec] = utc_to_doy(utc);
    [~, rho_vec] = atmosnrlmsise00(alt, lat, lon, year, doy, ut_sec);
    rho = rho_vec(6);                    % total mass density [kg/m^3]

    % Velocity relative to the co-rotating atmosphere, in ECI.
    omega_vec = [0; 0; p.omega_earth];
    v_rel     = v - cross(omega_vec, r);
    v_rel_n   = norm(v_rel);

    % Drag acceleration a = -1/2 * rho * (Cd*A/m) * |v_rel| * v_rel
    a_drag = -0.5 * rho * (p.Cd * p.A / p.m) * v_rel_n * v_rel;
    a = a + a_drag;
end

dX = [v; a];

end

% util
function utc = epoch_to_utc(epoch, t)
% Return a [Y M D h m s] vector t seconds after the given epoch.
if isdatetime(epoch)
    dt  = epoch + seconds(t);
    utc = [dt.Year, dt.Month, dt.Day, dt.Hour, dt.Minute, dt.Second];
else
    % epoch already a [Y M D h m s] vector -> add t via datetime arithmetic
    dt  = datetime(epoch(1), epoch(2), epoch(3), ...
                   epoch(4), epoch(5), epoch(6)) + seconds(t);
    utc = [dt.Year, dt.Month, dt.Day, dt.Hour, dt.Minute, dt.Second];
end
end

function [year, doy, ut_sec] = utc_to_doy(utc)
% Split a [Y M D h m s] vector into year, day-of-year, and UT seconds.
year   = utc(1);
doy    = day(datetime(utc(1), utc(2), utc(3)), 'dayofyear');
ut_sec = utc(4)*3600 + utc(5)*60 + utc(6);
end
