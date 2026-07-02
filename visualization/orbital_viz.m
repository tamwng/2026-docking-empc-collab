function orbital_viz()
% orbital_viz.m — trajectory visualizer for thesis
% Imports a JSON file and replays the recorded ECI trajectory
% JSON format: { "metadata": {...}, "t": [N floats], "X": [Nx6 floats] }
% X columns: [rx, ry, rz, vx, vy, vz] in SI units (m, m/s)
% Reference: "Analytical Mechanics of Space Systems" - H.Schaub, J.L.Junkins

%%  CONSTANTS 
RE = 6378.137;      % Earth mean equatorial radius [km]

%%  PLAYBACK STATE
pb_t    = [];       % time vector [s], Nx1
pb_Xt    = [];      % target ECI state
pb_Rho   = [];      % chaser LVLH state
pb_idx  = 1;        % current frame index
pb_meta = [];       % metadata struct from JSON
running = false;    % is animation currently playing?

pb_tick_acc = 0;        % tick accumulator for real-time pacing
prox_zoom   = false;    % true = center ECI view on target, auto-scale to range
r_max_orbit = 1;        % [km] full-orbit axis half-range, set on JSON load

%% FIGURE 
fig = figure('Name', 'RPO Trajectory Visualizer', ...
    'Color', [0.02 0.05 0.09], ...
    'Position', [50 50 1400 700], ...
    'NumberTitle', 'off');

%% 3D AXES 
ax3d = axes('Parent', fig, ...
    'Position', [0.03 0.12 0.40 0.82], ...
    'Color',    [0.02 0.05 0.09], ...
    'XColor',   [0.29 0.42 0.53], ...
    'YColor',   [0.29 0.42 0.53], ...
    'ZColor',   [0.29 0.42 0.53], ...
    'FontSize', 7, 'FontName', 'Courier New');
hold(ax3d, 'on'); grid(ax3d, 'on'); axis(ax3d, 'equal');
view(ax3d, 35, 25);
xlabel(ax3d, 'Radial x [m]',      'FontSize', 7, 'Color', [0.69 0.80 0.91]);
ylabel(ax3d, 'Along-track y [m]', 'FontSize', 7, 'Color', [0.69 0.80 0.91]);
zlabel(ax3d, 'Cross-track z [m]', 'FontSize', 7, 'Color', [0.69 0.80 0.91]);
title(ax3d, 'HILL FRAME — 3D RELATIVE TRAJECTORY', ...
    'Color', [0 0.78 1], 'FontName', 'Courier New', 'FontSize', 9);
 
%%  2D PLOT 1: Radial vs Along-track (Hill frame)
ax_xy = axes('Parent', fig, ...
    'Position', [0.46 0.55 0.24 0.38], ...
    'Color',    [0.02 0.05 0.09], ...
    'XColor',   [0.29 0.42 0.53], ...
    'YColor',   [0.29 0.42 0.53], ...
    'FontSize', 7, 'FontName', 'Courier New');
hold(ax_xy, 'on'); grid(ax_xy, 'on'); axis(ax_xy, 'equal');
xlabel(ax_xy, 'Along-track y [m]', 'FontSize', 7, 'Color', [0.69 0.80 0.91]);
ylabel(ax_xy, 'Radial x [m]',      'FontSize', 7, 'Color', [0.69 0.80 0.91]);
title(ax_xy, 'HILL FRAME — RADIAL vs ALONG-TRACK', ...
    'Color', [0 0.78 1], 'FontName', 'Courier New', 'FontSize', 8);
 
%% 2D PLOT 2: Relative distance vs time 
ax_dt = axes('Parent', fig, ...
    'Position', [0.46 0.12 0.24 0.34], ...
    'Color',    [0.02 0.05 0.09], ...
    'XColor',   [0.29 0.42 0.53], ...
    'YColor',   [0.29 0.42 0.53], ...
    'FontSize', 7, 'FontName', 'Courier New');
hold(ax_dt, 'on'); grid(ax_dt, 'on');
xlabel(ax_dt, 'Time [s]',             'FontSize', 7, 'Color', [0.69 0.80 0.91]);
ylabel(ax_dt, 'Relative distance [m]', 'FontSize', 7, 'Color', [0.69 0.80 0.91]);
title(ax_dt, 'RELATIVE DISTANCE vs TIME', ...
    'Color', [0 0.78 1], 'FontName', 'Courier New', 'FontSize', 8);

% Origin marker — target is always at (0,0,0) in Hill frame
plot3(ax3d, 0, 0, 0, '+', 'Color', [0 0.78 1], 'MarkerSize', 12, 'LineWidth', 1.5);

%%  PLOT ARTISTS  (Hill-frame 3D)
% Full mission relative trajectory — very dim static trail drawn on load
h_trail_t = plot3(ax3d, nan, nan, nan, 'Color', [0.07 0.14 0.20], ...
                  'LineWidth', 0.4, 'LineStyle', '-');   % unused; kept for handle parity
h_sc_t    = plot3(ax3d, nan, nan, nan, 'o', 'Color', [0 0.78 1], ...
                  'MarkerFaceColor', [0 0.78 1], 'MarkerSize', 6);  % target (origin)
h_rad_t   = plot3(ax3d, nan, nan, nan, 'Color', [0.05 0.10 0.15], 'LineWidth', 0.3);

% Chaser: dim full-path trail, bright animated dot
h_trail_c = plot3(ax3d, nan, nan, nan, 'Color', [0.65 0.35 0.06], ...
                  'LineWidth', 0.8, 'LineStyle', '-');
h_sc_c    = plot3(ax3d, nan, nan, nan, 'o', 'Color', [1 0.55 0], ...
                  'MarkerFaceColor', [1 0.55 0], 'MarkerSize', 8);
h_rad_c   = plot3(ax3d, nan, nan, nan, 'Color', [0.10 0.06 0.02], 'LineWidth', 0.3);

% Range line: origin → chaser
h_los = plot3(ax3d, nan, nan, nan, 'Color', [0.20 0.28 0.35], ...
              'LineWidth', 0.7, 'LineStyle', ':');

legend(ax3d, [h_sc_t, h_sc_c], {'Target', 'Chaser'}, ...
    'TextColor', [0.69 0.80 0.91], 'Color', [0.04 0.07 0.13], ...
    'EdgeColor', [0.10 0.18 0.28], 'FontSize', 7, 'Location', 'northwest');

%% PLOT ARTISTS — 2D Hill frame 
% Full trajectory trail (static, drawn on load)
h_xy_trail = plot(ax_xy, nan, nan, 'Color', [0 0.78 1]*0.5, ...
                  'LineWidth', 0.8, 'LineStyle', '-');
% Live position dot
h_xy_dot   = plot(ax_xy, nan, nan, 'o', 'Color', [1 0.55 0], ...
                  'MarkerFaceColor', [1 0.55 0], 'MarkerSize', 7);
% Origin cross (target position)
plot(ax_xy, 0, 0, '+', 'Color', [0 0.78 1], 'MarkerSize', 10, 'LineWidth', 1.5);
 
%% PLOT ARTISTS — 2D distance vs time
% Full distance history (static, drawn on load)
h_dt_trail = plot(ax_dt, nan, nan, 'Color', [0 0.78 1]*0.5, ...
                  'LineWidth', 0.8);
% Live progress marker
h_dt_dot   = plot(ax_dt, nan, nan, 'o', 'Color', [1 0.55 0], ...
                  'MarkerFaceColor', [1 0.55 0], 'MarkerSize', 7);

%%  TELEMETRY PANEL 
bg  = [0.02 0.05 0.09];
dim = [0.29 0.42 0.53];
blu = [0.00 0.78 1.00];
org = [1.00 0.55 0.00];

% Target telemetry (left column)
tnames_t = {'TARGET ALT', 'TARGET VEL', 'PERIOD', 'ELAPSED TIME'};
tx_t = [0.74 0.74 0.74 0.74];
ty_t = [0.84 0.72 0.60 0.48];
h_tval_t = gobjects(1, 4);
for k = 1:4
    annotation(fig, 'textbox', [tx_t(k) ty_t(k)+0.04 0.12 0.04], ...
        'String', tnames_t{k}, 'Color', dim, 'FontName', 'Courier New', ...
        'FontSize', 7, 'EdgeColor', 'none', 'BackgroundColor', bg);
    h_tval_t(k) = annotation(fig, 'textbox', [tx_t(k) ty_t(k)-0.01 0.12 0.05], ...
        'String', '-', 'Color', blu, 'FontName', 'Courier New', ...
        'FontSize', 11, 'EdgeColor', 'none', 'BackgroundColor', bg);
end
 
% Chaser / relative telemetry (right column) 
tnames_c = {'CHASER ALT', 'CHASER VEL', 'REL DISTANCE', 'REL SPEED'};
tx_c = [0.88 0.88 0.88 0.88];
ty_c = [0.84 0.72 0.60 0.48];
h_tval_c = gobjects(1, 4);
for k = 1:4
    annotation(fig, 'textbox', [tx_c(k) ty_c(k)+0.04 0.12 0.04], ...
        'String', tnames_c{k}, 'Color', dim, 'FontName', 'Courier New', ...
        'FontSize', 7, 'EdgeColor', 'none', 'BackgroundColor', bg);
    h_tval_c(k) = annotation(fig, 'textbox', [tx_c(k) ty_c(k)-0.01 0.12 0.05], ...
        'String', '-', 'Color', org, 'FontName', 'Courier New', ...
        'FontSize', 11, 'EdgeColor', 'none', 'BackgroundColor', bg);
end

% h_file_label = shows the loaded filename
h_file_label = annotation(fig, 'textbox', [0.74 0.91 0.24 0.05], ...
    'String', 'No file loaded', 'Color', dim, ...
    'FontName', 'Courier New', 'FontSize', 8, ...
    'EdgeColor', 'none', 'BackgroundColor', bg);

%% CONTROLS 
sc_bg = [0.02 0.05 0.09];

% btn_load  = open file browser to import JSON
btn_load = uicontrol(fig, 'Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', [0.74 0.025 0.055 0.045], 'String', 'LOAD JSON', ...
    'BackgroundColor', [0.04 0.08 0.14], 'ForegroundColor', [0 0.78 1], ...
    'FontName', 'Courier New', 'FontSize', 8);

% btn_play  = start / resume playback
btn_play = uicontrol(fig, 'Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', [0.80 0.025 0.055 0.045], 'String', 'PLAY', ...
    'BackgroundColor', [0.04 0.08 0.14], 'ForegroundColor', [0.22 1 0.08], ...
    'FontName', 'Courier New', 'FontSize', 8);

% btn_pause = pause playback
btn_pause = uicontrol(fig, 'Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', [0.86 0.025 0.055 0.045], 'String', 'PAUSE', ...
    'BackgroundColor', [0.04 0.08 0.14], 'ForegroundColor', [1 0.55 0], ...
    'FontName', 'Courier New', 'FontSize', 8);

% btn_reset = rewind to first frame
btn_reset = uicontrol(fig, 'Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', [0.92 0.025 0.055 0.045], 'String', 'RESET', ...
    'BackgroundColor', [0.04 0.08 0.14], 'ForegroundColor', [0.69 0.80 0.91], ...
    'FontName', 'Courier New', 'FontSize', 8);

% btn_zoom = toggle between full-orbit ECI view and proximity view centred on target
btn_zoom = uicontrol(fig, 'Style', 'togglebutton', 'Units', 'normalized', ...
    'Position', [0.03 0.025 0.09 0.045], 'String', 'PROX ZOOM', ...
    'BackgroundColor', [0.04 0.08 0.14], 'ForegroundColor', [0.29 0.42 0.53], ...
    'FontName', 'Courier New', 'FontSize', 8);

% sl_spd  = playback speed slider (1x to 200x)
% vh_spd  = label showing current speed
uicontrol(fig, 'Style', 'text', 'Units', 'normalized', ...
    'Position', [0.74 0.110 0.05 0.025], 'String', 'Speed', ...
    'HorizontalAlignment', 'left', 'BackgroundColor', sc_bg, ...
    'ForegroundColor', [0.69 0.80 0.91], 'FontName', 'Courier New', 'FontSize', 8);
sl_spd = uicontrol(fig, 'Style', 'slider', 'Units', 'normalized', ...
    'Position', [0.74 0.085 0.20 0.022], ...
    'Min', 1, 'Max', 200, 'Value', 1, ...
    'BackgroundColor', [0.06 0.13 0.22]);
vh_spd = uicontrol(fig, 'Style', 'text', 'Units', 'normalized', ...
    'Position', [0.93 0.085 0.05 0.022], 'String', '1x', ...
    'HorizontalAlignment', 'left', 'BackgroundColor', sc_bg, ...
    'ForegroundColor', [0.22 1 0.08], 'FontName', 'Courier New', 'FontSize', 8);

set(sl_spd, 'Callback', @(~,~) update_speed());
    function update_speed()
        vh_spd.String = sprintf('%dx', round(sl_spd.Value));
    end


%% LVLH TO ECI ROTATION 
% lvlh_to_eci: converts a position vector from LVLH frame to ECI frame
% The LVLH (Local Vertical Local Horizontal) frame is defined by:
%   x = radial     — points from Earth center outward through target
%   y = along-track — in direction of orbital motion
%   z = cross-track — normal to orbit plane (= orbit angular momentum direction)
%
% pos_t_m = target ECI position [m], 3x1
% vel_t_m = target ECI velocity [m/s], 3x1
% rho_m   = relative position in LVLH [m], 3x1
% returns chaser ECI position [m], 3x1
 
    function pos_chaser = lvlh_to_eci(pos_t_m, vel_t_m, rho_m)
        r_hat = pos_t_m / norm(pos_t_m);               % radial unit vector
        h_vec = cross(pos_t_m, vel_t_m);               % orbit angular momentum vector
        h_hat = h_vec / norm(h_vec);                   % cross-track unit vector
        y_hat = cross(h_hat, r_hat);                   % along-track unit vector
        R = [r_hat, y_hat, h_hat];                     % rotation matrix LVLH -> ECI (3x3)
        pos_chaser = pos_t_m + R * rho_m;              % chaser ECI position [m]
    end    

%%  PLAYBACK UPDATE 
% pb_update: reads current frame from pb_X, moves spacecraft, updates telemetry


    function pb_update()
        if isempty(pb_Xt), return; end
 
        % loop back to start
        if pb_idx > size(pb_Xt, 1)
            pb_idx = 1;
        end
 
        %  Target state 
        row_t    = pb_Xt(pb_idx, :);
        pos_t_m  = row_t(1:3)';         % target ECI position [m], column vector
        vel_t_m  = row_t(4:6)';         % target ECI velocity [m/s], column vector
        pos_t_km = pos_t_m / 1000;      % convert to km
        vel_t_kms = vel_t_m / 1000;
 
        %  Chaser state 
        row_rho  = pb_Rho(pb_idx, :);
        rho_m    = row_rho(1:3)';       % relative position in LVLH [m], column vector
        rho_dot_m = row_rho(4:6)';      % relative velocity in LVLH [m/s]
 
        % convert chaser to ECI
        pos_c_m  = lvlh_to_eci(pos_t_m, vel_t_m, rho_m);
        pos_c_km = pos_c_m / 1000;      % chaser ECI position [km]
 
        % Update 3D — Hill frame: target at origin, chaser at rho_m [m]
        set(h_sc_t,  'XData', 0,        'YData', 0,        'ZData', 0);
        set(h_rad_t, 'XData', nan,      'YData', nan,      'ZData', nan);
        set(h_sc_c,  'XData', rho_m(1), 'YData', rho_m(2), 'ZData', rho_m(3));
        set(h_rad_c, 'XData', nan,      'YData', nan,      'ZData', nan);
        set(h_los,   'XData', [0 rho_m(1)], ...
                     'YData', [0 rho_m(2)], ...
                     'ZData', [0 rho_m(3)]);
 
        % Update Plot 1: live dot on x-y phase plot
        set(h_xy_dot, 'XData', rho_m(2), 'YData', rho_m(1));
 
        % Update Plot 2: live dot on distance-time plot
        t_now = pb_t(pb_idx);
        rel_d = norm(rho_m);
        set(h_dt_dot, 'XData', t_now, 'YData', rel_d);
 
 
        % Telemetry 
        t_alt  = norm(pos_t_km) - RE;       % target altitude [km]
        t_spd  = norm(vel_t_kms);           % target speed [km/s]
        c_alt  = norm(pos_c_km) - RE;       % chaser altitude [km]
        rel_d  = norm(rho_m);               % relative distance [m]
        rel_v  = norm(rho_dot_m);           % relative speed [m/s]
        t_now  = pb_t(pb_idx);              % elapsed time [s]
 
        T_min = nan;
        if ~isempty(pb_meta) && isstruct(pb_meta)
            T_min = pb_meta.period_min;
        end
 
        % target column (cyan)
        vals_t = {sprintf('%.0f km',  t_alt), ...
                  sprintf('%.3f km/s', t_spd), ...
                  sprintf('%.2f min',  T_min), ...
                  sprintf('%.1f s',    t_now)};
        for k = 1:4
            h_tval_t(k).String = vals_t{k};
        end
 
        % chaser / relative column (orange)
        vals_c = {sprintf('%.0f km',  c_alt), ...
                  sprintf('%.1f m/s', norm(rho_dot_m + vel_t_m/1000)*1000), ...
                  sprintf('%.1f m',   rel_d), ...
                  sprintf('%.2f m/s', rel_v)};
        for k = 1:4
            h_tval_c(k).String = vals_c{k};
        end
 
        % Proximity zoom: auto-scale Hill-frame view to 3× current range [m]
        if prox_zoom
            d_m = max(norm(rho_m) * 3, 50);   % 3× range in m, min 50 m
            axis(ax3d, [-d_m d_m -d_m d_m -d_m d_m]);
        end

        % Frame advance
        sim_dt       = pb_t(2) - pb_t(1);
        ticks_needed = sim_dt / (0.05 * sl_spd.Value);
        pb_tick_acc  = pb_tick_acc + 1;
        if pb_tick_acc >= ticks_needed
            pb_idx      = pb_idx + max(1, round(pb_tick_acc / ticks_needed));
            pb_tick_acc = 0;
        end

        drawnow limitrate;
    end

%% JSON LOAD 

    function load_json(src, evt)
        [fname, fpath] = uigetfile('*.json', 'Select RPO trajectory JSON file');
        if isequal(fname, 0), return; end
 
        fid = fopen(fullfile(fpath, fname), 'r');
        raw = fread(fid, inf, 'uint8=>char')';
        fclose(fid);
        data = jsondecode(raw);
 
        % validate required fields
        if ~isfield(data, 't') || ~isfield(data, 'X_t') || ~isfield(data, 'Rho')
            errordlg('JSON must contain fields "t", "X_t" and "Rho"', 'Invalid file');
            return;
        end
 
        pb_t        = data.t;
        pb_Xt       = data.X_t;
        pb_Rho      = data.Rho;
        pb_idx      = 1;
        pb_tick_acc = 0;
        pb_meta     = [];
        if isfield(data, 'metadata')
            pb_meta = data.metadata;
        end
 
        % Draw full Hill-frame relative trajectory as dim static trail [m]
        rho_xyz = pb_Rho(:, 1:3);     % N × 3 [m]
        set(h_trail_t, 'XData', nan, 'YData', nan, 'ZData', nan);  % unused
        set(h_trail_c, 'XData', rho_xyz(:,1)', ...
                       'YData', rho_xyz(:,2)', ...
                       'ZData', rho_xyz(:,3)');

        % Fit axis to relative trajectory extent [m]
        r_hill      = max(vecnorm(rho_xyz, 2, 2)) * 1.2;
        r_max_orbit = r_hill;
        axis(ax3d, [-r_hill r_hill -r_hill r_hill -r_hill r_hill]);

        % Draw full Hill frame trail on Plot 1 (x vs y)
        rho_y = pb_Rho(:, 2);   % along-track
        rho_x = pb_Rho(:, 1);   % radial
        set(h_xy_trail, 'XData', rho_y', 'YData', rho_x');
        axis(ax_xy, 'auto');
 
        % Draw full distance history on Plot 2
        rel_dist = vecnorm(pb_Rho(:, 1:3), 2, 2);   % Nx1
        set(h_dt_trail, 'XData', pb_t', 'YData', rel_dist');
        axis(ax_dt, 'auto');
 
        h_file_label.String = fname;
        h_file_label.Color  = [0 0.78 1];
 
        running = true;
    end

%% TIMER 

    function tick(src, evt)
        if running
            pb_update();
        end
    end

    function on_close(src, evt)
        stop(t); delete(t); delete(fig);
    end

%% ── BUTTON CALLBACKS ─────────────────────────────────────────────────────
set(btn_load,  'Callback', @load_json);
set(btn_play,  'Callback', @(~,~) set_running(true));
set(btn_pause, 'Callback', @(~,~) set_running(false));
set(btn_reset, 'Callback', @(~,~) do_reset());
set(btn_zoom,  'Callback', @toggle_zoom);
set(fig,       'CloseRequestFcn', @on_close);

    function toggle_zoom(src, ~)
        prox_zoom = logical(src.Value);
        if prox_zoom
            src.ForegroundColor = [1 0.55 0];   % orange = proximity mode active
        else
            src.ForegroundColor = [0.29 0.42 0.53];
            if r_max_orbit > 0
                axis(ax3d, [-r_max_orbit r_max_orbit -r_max_orbit r_max_orbit ...
                            -r_max_orbit r_max_orbit]);
            end
        end
    end

    function set_running(val)
        running = val;
    end

    function do_reset()
        pb_idx  = 1;
        pb_tick_acc = 0;
        running = false;
    end

%% ── START ────────────────────────────────────────────────────────────────
t = timer('ExecutionMode', 'fixedRate', 'Period', 0.05, 'TimerFcn', @tick);
start(t);

end