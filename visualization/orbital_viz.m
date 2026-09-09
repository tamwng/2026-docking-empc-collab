function orbital_viz()
% orbital_viz.m — trajectory visualizer for thesis
% Imports a JSON file and replays the recorded ECI trajectory
% JSON format: { "metadata": {...}, "t": [N floats], "X": [Nx6 floats] }
% X columns: [rx, ry, rz, vx, vy, vz] in SI units (m, m/s)
% Reference: "Analytical Mechanics of Space Systems" - H.Schaub, J.L.Junkins
%
% ECI panel offset exaggeration: the EXAG slider inflates only the drawn
% chaser offset from the target (trail + live dot + LOS line) around the
% true ECI orbits — it never touches the JSON data, the Hill-frame panel,
% or the telemetry readouts, which always show the true, unscaled
% separation. Default is 1x (true scale) precisely so the true separation
% is what's shown unless deliberately exaggerated.

%%  CONSTANTS
RE = 6378.137;      % Earth mean equatorial radius [km]

%%  PLAYBACK STATE
pb_t          = [];       % time vector [s], Nx1
pb_Xt         = [];       % target ECI state, Nx6
pb_Rho        = [];       % chaser LVLH state, Nx6
pb_pos_t_eci  = [];       % target ECI position, Nx3 [km]
pb_offset_eci = [];       % true chaser-minus-target ECI-aligned offset, Nx3 [km]
pb_idx        = 1;        % current frame index
pb_meta       = [];       % metadata struct from JSON
running       = false;    % is animation currently playing?

pb_tick_acc = 0;        % tick accumulator for real-time pacing
prox_zoom   = false;    % true = center ECI view on target, auto-scale to range

%% FIGURE
fig = figure('Name', 'RPO Trajectory Visualizer', ...
    'Color', [0.02 0.05 0.09], ...
    'Position', [50 50 1400 700], ...
    'NumberTitle', 'off');

%% 3D AXES — ECI FRAME (ABSOLUTE ORBITS, WITH EARTH)
ax_eci = axes('Parent', fig, ...
    'Position', [0.03 0.12 0.42 0.82], ...
    'Color',    [0.02 0.05 0.09], ...
    'XColor',   [0.29 0.42 0.53], ...
    'YColor',   [0.29 0.42 0.53], ...
    'ZColor',   [0.29 0.42 0.53], ...
    'FontSize', 7, 'FontName', 'Courier New');
hold(ax_eci, 'on'); grid(ax_eci, 'on'); axis(ax_eci, 'equal');
view(ax_eci, 35, 25);
xlabel(ax_eci, 'X (km)', 'FontSize', 7, 'Color', [0.69 0.80 0.91]);
ylabel(ax_eci, 'Y (km)', 'FontSize', 7, 'Color', [0.69 0.80 0.91]);
zlabel(ax_eci, 'Z (km)', 'FontSize', 7, 'Color', [0.69 0.80 0.91]);
title(ax_eci, 'ECI FRAME — ABSOLUTE ORBITS', ...
    'Color', [0 0.78 1], 'FontName', 'Courier New', 'FontSize', 9);

% Earth sphere
[ue, ve] = meshgrid(linspace(0, 2*pi, 36), linspace(0, pi, 18));
mesh(ax_eci, RE*sin(ve).*cos(ue), RE*sin(ve).*sin(ue), RE*cos(ve), ...
    'EdgeColor', [0.10 0.23 0.35], 'FaceColor', [0.04 0.08 0.14], ...
    'EdgeAlpha', 0.4, 'FaceAlpha', 0.6);

%% PLOT ARTISTS — ECI FRAME (target + chaser absolute orbits around Earth)
h_trail_t_eci = plot3(ax_eci, nan, nan, nan, 'Color', [0 0.78 1], ...
                  'LineWidth', 1.0, 'LineStyle', '-');
h_sc_t_eci    = plot3(ax_eci, nan, nan, nan, 'o', 'Color', [0 0.78 1], ...
                  'MarkerFaceColor', [0 0.78 1], 'MarkerSize', 8);
h_rad_t_eci   = plot3(ax_eci, nan, nan, nan, 'Color', [0.10 0.25 0.35], 'LineWidth', 0.6);

h_trail_c_eci = plot3(ax_eci, nan, nan, nan, 'Color', [1 0.55 0], ...
                  'LineWidth', 1.0, 'LineStyle', '--');
h_sc_c_eci    = plot3(ax_eci, nan, nan, nan, 'o', 'Color', [1 0.55 0], ...
                  'MarkerFaceColor', [1 0.55 0], 'MarkerSize', 8);
h_rad_c_eci   = plot3(ax_eci, nan, nan, nan, 'Color', [0.30 0.18 0.05], 'LineWidth', 0.6);

h_los_eci = plot3(ax_eci, nan, nan, nan, 'Color', [0.8 0.8 0.8], ...
              'LineWidth', 0.8, 'LineStyle', ':');

legend(ax_eci, [h_sc_t_eci, h_sc_c_eci], {'Target', 'Chaser'}, ...
    'TextColor', [0.69 0.80 0.91], 'Color', [0.04 0.07 0.13], ...
    'EdgeColor', [0.10 0.18 0.28], 'FontSize', 7, 'Location', 'northwest');

%%  2D PLOT: Radial vs Along-track (Hill frame)
ax_xy = axes('Parent', fig, ...
    'Position', [0.50 0.12 0.22 0.82], ...
    'Color',    [0.02 0.05 0.09], ...
    'XColor',   [0.29 0.42 0.53], ...
    'YColor',   [0.29 0.42 0.53], ...
    'FontSize', 7, 'FontName', 'Courier New');
hold(ax_xy, 'on'); grid(ax_xy, 'on'); axis(ax_xy, 'equal');
xlabel(ax_xy, 'Along-track y [m]', 'FontSize', 7, 'Color', [0.69 0.80 0.91]);
ylabel(ax_xy, 'Radial x [m]',      'FontSize', 7, 'Color', [0.69 0.80 0.91]);
title(ax_xy, 'HILL FRAME — RADIAL vs ALONG-TRACK', ...
    'Color', [0 0.78 1], 'FontName', 'Courier New', 'FontSize', 8);

%% PLOT ARTISTS — 2D Hill frame
% Full trajectory trail (static, drawn on load)
h_xy_trail = plot(ax_xy, nan, nan, 'Color', [0 0.78 1]*0.5, ...
                  'LineWidth', 0.8, 'LineStyle', '-');
% Live position dot
h_xy_dot   = plot(ax_xy, nan, nan, 'o', 'Color', [1 0.55 0], ...
                  'MarkerFaceColor', [1 0.55 0], 'MarkerSize', 7);
% Origin cross (target position)
plot(ax_xy, 0, 0, '+', 'Color', [0 0.78 1], 'MarkerSize', 10, 'LineWidth', 1.5);

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
    'Position', [0.74 0.010 0.055 0.045], 'String', 'LOAD JSON', ...
    'BackgroundColor', [0.04 0.08 0.14], 'ForegroundColor', [0 0.78 1], ...
    'FontName', 'Courier New', 'FontSize', 8);

% btn_play  = start / resume playback
btn_play = uicontrol(fig, 'Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', [0.80 0.010 0.055 0.045], 'String', 'PLAY', ...
    'BackgroundColor', [0.04 0.08 0.14], 'ForegroundColor', [0.22 1 0.08], ...
    'FontName', 'Courier New', 'FontSize', 8);

% btn_pause = pause playback
btn_pause = uicontrol(fig, 'Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', [0.86 0.010 0.055 0.045], 'String', 'PAUSE', ...
    'BackgroundColor', [0.04 0.08 0.14], 'ForegroundColor', [1 0.55 0], ...
    'FontName', 'Courier New', 'FontSize', 8);

% btn_reset = rewind to first frame
btn_reset = uicontrol(fig, 'Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', [0.92 0.010 0.055 0.045], 'String', 'RESET', ...
    'BackgroundColor', [0.04 0.08 0.14], 'ForegroundColor', [0.69 0.80 0.91], ...
    'FontName', 'Courier New', 'FontSize', 8);

% btn_zoom = toggle between full-orbit ECI view and proximity view centred on target
btn_zoom = uicontrol(fig, 'Style', 'togglebutton', 'Units', 'normalized', ...
    'Position', [0.03 0.010 0.09 0.045], 'String', 'PROX ZOOM', ...
    'BackgroundColor', [0.04 0.08 0.14], 'ForegroundColor', [0.29 0.42 0.53], ...
    'FontName', 'Courier New', 'FontSize', 8);

% sl_spd  = playback speed slider (1x to 200x)
uicontrol(fig, 'Style', 'text', 'Units', 'normalized', ...
    'Position', [0.74 0.135 0.05 0.025], 'String', 'Speed', ...
    'HorizontalAlignment', 'left', 'BackgroundColor', sc_bg, ...
    'ForegroundColor', [0.69 0.80 0.91], 'FontName', 'Courier New', 'FontSize', 8);
sl_spd = uicontrol(fig, 'Style', 'slider', 'Units', 'normalized', ...
    'Position', [0.74 0.110 0.20 0.022], ...
    'Min', 1, 'Max', 200, 'Value', 1, ...
    'BackgroundColor', [0.06 0.13 0.22]);
vh_spd = uicontrol(fig, 'Style', 'text', 'Units', 'normalized', ...
    'Position', [0.93 0.110 0.05 0.022], 'String', '1x', ...
    'HorizontalAlignment', 'left', 'BackgroundColor', sc_bg, ...
    'ForegroundColor', [0.22 1 0.08], 'FontName', 'Courier New', 'FontSize', 8);

set(sl_spd, 'Callback', @(~,~) update_speed());
    function update_speed()
        vh_spd.String = sprintf('%dx', round(sl_spd.Value));
    end

% sl_exag = ECI-panel offset exaggeration (1x = true scale, default)
% Only inflates the DRAWN chaser offset in the ECI panel (trail, live dot,
% radius line, LOS line). Telemetry, the Hill-frame panel and the JSON data
% are always the true, unscaled separation.
uicontrol(fig, 'Style', 'text', 'Units', 'normalized', ...
    'Position', [0.74 0.085 0.10 0.025], 'String', 'ECI Exaggeration', ...
    'HorizontalAlignment', 'left', 'BackgroundColor', sc_bg, ...
    'ForegroundColor', [0.69 0.80 0.91], 'FontName', 'Courier New', 'FontSize', 8);
sl_exag = uicontrol(fig, 'Style', 'slider', 'Units', 'normalized', ...
    'Position', [0.74 0.060 0.20 0.022], ...
    'Min', 1, 'Max', 10000, 'Value', 1, ...
    'BackgroundColor', [0.06 0.13 0.22]);
vh_exag = uicontrol(fig, 'Style', 'text', 'Units', 'normalized', ...
    'Position', [0.93 0.060 0.05 0.022], 'String', '1x', ...
    'HorizontalAlignment', 'left', 'BackgroundColor', sc_bg, ...
    'ForegroundColor', [1 0.55 0], 'FontName', 'Courier New', 'FontSize', 8);

set(sl_exag, 'Callback', @(~,~) redraw_eci_display());

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

%% ECI DISPLAY REDRAW
% redraw_eci_display: recomputes the ECI chaser TRAIL and axis fit at the
% current exaggeration factor. Called on JSON load, on EXAG slider move,
% and on PROX ZOOM toggle. Live per-frame dot position is handled inside
% pb_update() so it also reacts to the slider during playback.

    function redraw_eci_display()
        if isempty(pb_pos_t_eci), return; end
        exag = sl_exag.Value;
        vh_exag.String = sprintf('%.0fx', exag);

        pos_c_disp = pb_pos_t_eci + exag * pb_offset_eci;   % N x 3 [km]
        set(h_trail_c_eci, 'XData', pos_c_disp(:,1)', ...
                           'YData', pos_c_disp(:,2)', ...
                           'ZData', pos_c_disp(:,3)');

        true_max_sep_km = max(vecnorm(pb_offset_eci, 2, 2));
        if true_max_sep_km < 1
            sep_str = sprintf('%.0f m', true_max_sep_km * 1000);
        else
            sep_str = sprintf('%.2f km', true_max_sep_km);
        end
        title(ax_eci, sprintf('ECI FRAME — ABSOLUTE ORBITS  (offset \\times%.0f, true sep \\leq %s)', ...
            exag, sep_str), 'Color', [0 0.78 1], 'FontName', 'Courier New', 'FontSize', 9);

        if ~prox_zoom
            all_pos_eci = [pb_pos_t_eci; pos_c_disp];
            r_max_eci   = max(vecnorm(all_pos_eci, 2, 2)) * 1.15;
            axis(ax_eci, [-r_max_eci r_max_eci -r_max_eci r_max_eci -r_max_eci r_max_eci]);
        end
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

        % convert chaser to ECI (TRUE, unscaled — used for telemetry)
        pos_c_m  = lvlh_to_eci(pos_t_m, vel_t_m, rho_m);
        pos_c_km = pos_c_m / 1000;      % chaser ECI position [km]

        % ECI display position (exaggerated by the EXAG slider, default 1x = true)
        exag       = sl_exag.Value;
        offset_km  = pos_c_km - pos_t_km;
        pos_c_disp = pos_t_km + exag * offset_km;

        % Update 3D — ECI frame: absolute target & (possibly exaggerated) chaser
        set(h_sc_t_eci,  'XData', pos_t_km(1), 'YData', pos_t_km(2), 'ZData', pos_t_km(3));
        set(h_rad_t_eci, 'XData', [0 pos_t_km(1)], ...
                         'YData', [0 pos_t_km(2)], ...
                         'ZData', [0 pos_t_km(3)]);
        set(h_sc_c_eci,  'XData', pos_c_disp(1), 'YData', pos_c_disp(2), 'ZData', pos_c_disp(3));
        set(h_rad_c_eci, 'XData', [0 pos_c_disp(1)], ...
                         'YData', [0 pos_c_disp(2)], ...
                         'ZData', [0 pos_c_disp(3)]);
        set(h_los_eci,   'XData', [pos_t_km(1) pos_c_disp(1)], ...
                         'YData', [pos_t_km(2) pos_c_disp(2)], ...
                         'ZData', [pos_t_km(3) pos_c_disp(3)]);

        % Update Plot: live dot on x-y phase plot (true Hill-frame scale)
        set(h_xy_dot, 'XData', rho_m(2), 'YData', rho_m(1));

        % Telemetry (always TRUE, unscaled values)
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

        % Proximity zoom: center the ECI view on the target, auto-scale to
        % 3x the current DISPLAYED (possibly exaggerated) separation
        if prox_zoom
            d_km = max(norm(offset_km) * exag * 3, 0.05);   % [km], floor 50 m
            axis(ax_eci, [pos_t_km(1)-d_km pos_t_km(1)+d_km, ...
                          pos_t_km(2)-d_km pos_t_km(2)+d_km, ...
                          pos_t_km(3)-d_km pos_t_km(3)+d_km]);
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

        % Precompute full-mission ECI positions (target + TRUE chaser offset)
        pos_t_all = pb_Xt(:, 1:3) / 1000;      % N × 3 [km]
        n_rows    = size(pb_Xt, 1);
        pos_c_all = zeros(n_rows, 3);
        for k = 1:n_rows
            pt = pb_Xt(k, 1:3)';
            vt = pb_Xt(k, 4:6)';
            rc = pb_Rho(k, 1:3)';
            pos_c_all(k, :) = lvlh_to_eci(pt, vt, rc)' / 1000;   % [km]
        end
        pb_pos_t_eci  = pos_t_all;
        pb_offset_eci = pos_c_all - pos_t_all;   % true ECI-aligned offset [km]

        set(h_trail_t_eci, 'XData', pos_t_all(:,1)', ...
                           'YData', pos_t_all(:,2)', ...
                           'ZData', pos_t_all(:,3)');
        redraw_eci_display();

        % Draw full Hill-frame trail on the radial-vs-along-track plot [m]
        rho_y = pb_Rho(:, 2);   % along-track
        rho_x = pb_Rho(:, 1);   % radial
        set(h_xy_trail, 'XData', rho_y', 'YData', rho_x');
        axis(ax_xy, 'auto');

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
        end
        redraw_eci_display();   % restores full-orbit fit when turning off
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
