% sim_adrios_combinations.m
% Non-interactive enumeration of ALL controller combinations for the ADRIOS
% full mission (Closing × Flyaround × Docking).
%
% Phase configs mirror sim_adrios_full_mission.m exactly; this driver just
% removes the operator gates and loops over every combination, tabulating
% total Δv and total mission duration.
%
% Combinations:  Closing {A1, B, C} × Flyaround {Run2, Run5} × Docking {Run2, Run3}
%                = 3 × 2 × 2 = 12.
%
% Efficiency: closing depends only on its own choice (3 runs); flyaround on
% (closing, S2) → 6 runs; docking on (closing, S2, S1) → 12 runs. The
% expensive flyaround QPs (N = 92) are therefore solved 6×, not 12×.

clear; clc;

%% Path setup + constants
script_dir = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(script_dir, '..', 'tools', 'matlab2tikz-master', 'src')));
constants;   % n, T, a, i, alt, ...

%% Shared parameters (identical to sim_adrios_full_mission.m)
dt_S4 = 120;   N_S4 = 20;
dt_S2 = T/92;  P_S2 = 92;  N_S2 = 92;
dt_S1 = 10;    N_S1 = 24;   % min. feasible horizon for Run 3 term.eq @~150 m FKP
u_max    = 1e-2;
r_KOS_S4 = 50;
b_target = 75;
los_deg  = 30;
rho_c    = 1e-11;
eps_r5   = 1e-4;
P_avail  = 1e-5;

Ac_cwh = [0      0     0    1     0    0  ;
          0      0     0    0     1    0  ;
          0      0     0    0     0    1  ;
          3*n^2  0     0    0     2*n  0  ;
          0      0     0   -2*n   0    0  ;
          0      0    -n^2  0     0    0 ];
Bc_cwh = [zeros(3,3); eye(3)];

Q_dare = diag([1e-2, 1e-2, 1e-2, 1e0, 1e0, 1e0]);
R_dare = diag([1e4,  1e4,  1e4]);

[Ad_S4, Bd_S4] = discretize(Ac_cwh, Bc_cwh, dt_S4);
[~, P_clf_S4]  = lqr_controller(Ad_S4, Bd_S4, Q_dare, R_dare);
[Ad_S1, Bd_S1] = discretize(Ac_cwh, Bc_cwh, dt_S1);
[~, P_clf_S1]  = lqr_controller(Ad_S1, Bd_S1, Q_dare, R_dare);

%% S2 prediction matrices (choice-independent) — build once
[A_d_S2, B_d_S2] = discretize(Ac_cwh, Bc_cwh, dt_S2);

x_star_S2 = zeros(6, P_S2);
x_star_S2(:,1) = [b_target; 0; 0; 0; -2*n*b_target; 0];
for k_s = 1:P_S2-1
    x_star_S2(:, k_s+1) = A_d_S2 * x_star_S2(:, k_s);
end

Sx_S2 = zeros(6*N_S2, 6);  Su_S2 = zeros(6*N_S2, 3*N_S2);  Ad_pow = eye(6);
for ii = 1:N_S2
    Ad_pow = Ad_pow * A_d_S2;
    Sx_S2((ii-1)*6+1:ii*6, :) = Ad_pow;
    for jj = 1:ii
        Su_S2((ii-1)*6+1:ii*6, (jj-1)*3+1:jj*3) = A_d_S2^(ii-jj) * B_d_S2;
    end
end
pos_idx = zeros(1, 3*N_S2);
for ki = 1:N_S2, pos_idx((ki-1)*3+1:ki*3) = (ki-1)*6+1:(ki-1)*6+3; end

mats0.A_d = A_d_S2;  mats0.B_d = B_d_S2;
mats0.Sx  = Sx_S2;   mats0.Su  = Su_S2;
mats0.Sx_pos_all = Sx_S2(pos_idx, :);  mats0.Su_pos_all = Su_S2(pos_idx, :);
mats0.A_u = [eye(3*N_S2); -eye(3*N_S2)];  mats0.b_u = u_max * ones(6*N_S2, 1);

%% Controller labels
lab_S4 = {'A1 energy+ball', 'B dissipative', 'C energy+CLF'};
lab_S2 = {'R2 periodic-eq', 'R5 augmented'};
lab_S1 = {'R2 EMPC-CLF',    'R3 EMPC-termeq'};

%% Enumerate all combinations
n_comb = 3 * 2 * 2;
Res = struct('c4',{},'c2',{},'c1',{}, 'dv_S4',{},'dv_S2',{},'dv_dep',{},'dv_S1',{}, ...
             'dv_tot',{}, 'dur_S4',{},'dur_S2',{},'dur_S1',{},'dur_tot',{}, 'conv',{});
row = 0;
fprintf('\nEnumerating %d combinations...\n\n', n_comb);

for i4 = 1:3
    % ─── Closing ─────────────────────────────────────────────────────────
    con_S4 = struct();
    con_S4.u_max = u_max; con_S4.y_min_active = false; con_S4.los_cone.active = false;
    con_S4.use_avg_power = false; con_S4.avg_power.P_avail = P_avail;
    con_S4.use_passive_safety = false;
    params_S4 = struct();
    params_S4.x0 = [0;10000;0;0;0;0]; params_S4.x_switch = [0;300;0;0;0;0];
    params_S4.dt = dt_S4; params_S4.u_max = u_max; params_S4.t_final = 5*T;
    params_S4.r_KOS = r_KOS_S4; params_S4.verbose = false;
    switch i4
        case 1
            cost_S4 = struct('Q',zeros(6),'R',eye(3),'P',zeros(6));
            con_S4.use_terminal_ball = true; con_S4.terminal_ball.r_switch = 200;
            params_S4.N = N_S4;
        case 2
            cost_S4 = struct('Q',rho_c*eye(6),'R',eye(3),'P',zeros(6));
            con_S4.use_terminal_ball = false; params_S4.N = 40;
        case 3
            cost_S4 = struct('Q',zeros(6),'R',eye(3),'P',P_clf_S4);
            con_S4.use_terminal_ball = false; params_S4.N = N_S4;
    end
    res_S4 = run_closing(cost_S4, con_S4, params_S4);
    conv4 = ~isnan(res_S4.conv_step);
    if conv4, x_end_S4 = res_S4.X(:,res_S4.conv_step); steps_S4 = res_S4.conv_step;
    else,     x_end_S4 = res_S4.X(:,end);              steps_S4 = res_S4.n_ctrl;  end
    dv_S4 = res_S4.dv; dur_S4 = steps_S4 * dt_S4;

    for i2 = 1:2
        % ─── Flyaround ───────────────────────────────────────────────────
        mats = mats0; rcfg = struct();
        switch i2
            case 1
                mats.H = 2*eye(3*N_S2);
                rcfg.use_nmc_manifold = true;  rcfg.use_periodic_terminal = true;
            case 2
                Qbar = kron(eye(N_S2), eps_r5*eye(6));
                Hr5  = 2*(Su_S2'*Qbar*Su_S2 + eye(3*N_S2));
                mats.H = (Hr5+Hr5')/2;
                mats.F_mat = 2*Su_S2'*Qbar*Sx_S2;  mats.f_ref_mat = -2*Su_S2'*Qbar;
                rcfg.use_nmc_manifold = false; rcfg.use_periodic_terminal = false;
        end
        rcfg.x_star = x_star_S2; rcfg.n_sim_override = []; rcfg.use_band = false;
        params_S2 = struct();
        params_S2.x0 = x_end_S4; params_S2.rho_min = 50; params_S2.rho_max = 200;
        params_S2.u_max = u_max; params_S2.N = N_S2; params_S2.P = P_S2;
        params_S2.n = n; params_S2.dt = dt_S2; params_S2.n_sim = 6*P_S2;
        res_S2 = run_flyaround(rcfg, mats, params_S2);

        dv_S2 = res_S2.dv_total; n_sim_S2 = res_S2.n_sim; b_est = res_S2.b_est;
        % V-bar crossing + departure burn (mirror mission)
        last_half = max(1, n_sim_S2 - round(P_S2/2));
        xs = abs(res_S2.x_log(1, last_half:end));  ys = res_S2.x_log(2, last_half:end);
        xs(ys < b_est*0.5) = inf;  [~, krel] = min(xs);  k_vbar = last_half + krel - 1;
        x_dep = res_S2.x_log(:, k_vbar);
        dv_dep = abs(x_dep(4));  x_FKP = [0; x_dep(2); 0; 0; 0; 0];
        dur_S2 = (k_vbar-1) * dt_S2;

        for i1 = 1:2
            % ─── Docking ─────────────────────────────────────────────────
            con_S1 = struct();
            con_S1.u_max = u_max; con_S1.y_min_active = true; con_S1.los_cone.active = true;
            con_S1.los_cone.half_angle = deg2rad(los_deg);
            con_S1.use_passive_safety = false; con_S1.terminal_eq = false;
            params_S1 = struct();
            params_S1.x0 = x_FKP; params_S1.x_dock = zeros(6,1);
            params_S1.dt = dt_S1; params_S1.N = N_S1; params_S1.u_max = u_max;
            params_S1.t_final = 3*T; params_S1.r_KOS = 5; params_S1.conv_tol = 1.0;
            params_S1.verbose = false;
            switch i1
                case 1, cost_S1 = struct('Q',zeros(6),'R',eye(3),'P',P_clf_S1);
                case 2, cost_S1 = struct('Q',zeros(6),'R',eye(3),'P',zeros(6));
                        con_S1.terminal_eq = true;
            end
            res_S1 = run_docking(cost_S1, con_S1, params_S1);
            conv1 = ~isnan(res_S1.conv_step);
            if conv1, steps_S1 = res_S1.conv_step; else, steps_S1 = res_S1.n_ctrl; end
            dv_S1 = res_S1.dv; dur_S1 = steps_S1 * dt_S1;

            row = row + 1;
            Res(row).c4 = lab_S4{i4}; Res(row).c2 = lab_S2{i2}; Res(row).c1 = lab_S1{i1};
            Res(row).dv_S4 = dv_S4; Res(row).dv_S2 = dv_S2;
            Res(row).dv_dep = dv_dep; Res(row).dv_S1 = dv_S1;
            Res(row).dv_tot = dv_S4 + dv_S2 + dv_dep + dv_S1;
            Res(row).dur_S4 = dur_S4; Res(row).dur_S2 = dur_S2; Res(row).dur_S1 = dur_S1;
            Res(row).dur_tot = dur_S4 + dur_S2 + dur_S1;
            Res(row).conv = conv4 && conv1;
            fprintf('[%2d/%d] S4=%-14s S2=%-14s S1=%-14s | dv=%6.3f m/s | dur=%5.2f h%s\n', ...
                row, n_comb, lab_S4{i4}, lab_S2{i2}, lab_S1{i1}, ...
                Res(row).dv_tot, Res(row).dur_tot/3600, ternary(Res(row).conv,'',' [!conv]'));
        end
    end
end

%% Results table
T_all = table( ...
    {Res.c4}', {Res.c2}', {Res.c1}', ...
    round([Res.dv_S4]',4), round([Res.dv_S2]',4), round([Res.dv_dep]',4), round([Res.dv_S1]',4), ...
    round([Res.dv_tot]',4), round([Res.dur_tot]'/3600,3), ...
    'VariableNames', {'Closing','Flyaround','Docking', ...
    'dvS4','dvS2','dvDep','dvS1','dv_total_ms','dur_total_h'});
fprintf('\n══════════════════ ADRIOS combination sweep ══════════════════\n');
disp(T_all);

[~, i_best_dv]  = min([Res.dv_tot]);
[~, i_best_dur] = min([Res.dur_tot]);
fprintf('Min Δv : #%d  (%.4f m/s)  — %s | %s | %s\n', i_best_dv, Res(i_best_dv).dv_tot, ...
    Res(i_best_dv).c4, Res(i_best_dv).c2, Res(i_best_dv).c1);
fprintf('Min dur: #%d  (%.2f h)     — %s | %s | %s\n', i_best_dur, Res(i_best_dur).dur_tot/3600, ...
    Res(i_best_dur).c4, Res(i_best_dur).c2, Res(i_best_dur).c1);

%% Exports
out_dir = fullfile(script_dir, '..', 'results', 'figures', 's_full_mission');
[~,~] = mkdir(out_dir);

% LaTeX (booktabs)
tex_path = fullfile(out_dir, 'combinations_table.tex');
ftex = fopen(tex_path, 'w');
fprintf(ftex, '%% Auto-generated by sim_adrios_combinations.m\n');
fprintf(ftex, '\\begin{tabular}{clllrr}\n\\toprule\n');
fprintf(ftex, '\\# & Closing & Flyaround & Docking & $\\Delta v_{\\mathrm{tot}}$ [m/s] & Duration [h] \\\\\n\\midrule\n');
for r = 1:n_comb
    fprintf(ftex, '%d & %s & %s & %s & %.4f & %.2f \\\\\n', r, ...
        Res(r).c4, Res(r).c2, Res(r).c1, Res(r).dv_tot, Res(r).dur_tot/3600);
end
fprintf(ftex, '\\bottomrule\n\\end{tabular}\n');
fclose(ftex);
fprintf('\nLaTeX table -> %s\n', tex_path);

% JSON
json_path = fullfile(script_dir, '..', 'exports', 'scenarios', 'sim_adrios_combinations.json');
jd.metadata = struct('date', char(datetime('now','Format','yyyy-MM-dd HH:mm')), ...
    'n_combinations', n_comb, 'dt_S4', dt_S4, 'dt_S2', dt_S2, 'dt_S1', dt_S1, ...
    'note', 'Total duration is dominated by the ~6-orbit flyaround maintenance phase.');
jd.combinations = Res;
fid = fopen(json_path, 'w'); fprintf(fid, '%s', jsonencode(jd, 'PrettyPrint', true)); fclose(fid);
fprintf('JSON       -> %s\n\n', json_path);

%% Local helper
function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end