function cfg = ATLAS_int()
%ATLAS_INT  Interactive configuration wizard for ATLAS v1.0.7
%
%  Usage:
%      cfg = ATLAS_int();
%      results = ATLAS(cfg);
%
%  Scoring in ATLAS: ARETE × MOIRA
%    ARETE — Above-threshold Resonance Extent and Target Efficacy
%    MOIRA — Measure of Intercepted Resonance and Absorption

fprintf('\n');
fprintf('  +================================================+\n');
fprintf('  |   ATLAS v1.0.7 — Interactive Setup             |\n');
fprintf('  |   Scoring: ARETE × MOIRA                       |\n');
fprintf('  +================================================+\n\n');

cfg = default_cfg();

%% [0] Propagation Mode Selection
fprintf('  [0] PROPAGATION MODE\n');
fprintf('      Available datasets:\n');
fprintf('        (1) O1 — Fundamental O-mode\n');
fprintf('        (2) O2 — Second harmonic O-mode\n');
fprintf('        (3) X1 — Fundamental X-mode\n');
fprintf('        (4) X2 — Second harmonic X-mode\n');
fprintf('        (A) All — Cross-mode comparison\n');
fprintf('        (C) Custom selection\n');
ms = upper(strtrim(input('      [1/2/3/4/A/C] (default 1): ', 's')));

all_csv  = {'data_kinetic_O1.csv','data_kinetic_O2.csv','data_kinetic_X1.csv','data_kinetic_X2.csv'};
all_lbls = {'O1','O2','X1','X2'};

switch ms
    case '2'
        cfg.csvPath = all_csv{2}; fprintf('      -> O2\n');
    case '3'
        cfg.csvPath = all_csv{3}; fprintf('      -> X1\n');
    case '4'
        cfg.csvPath = all_csv{4}; fprintf('      -> X2\n');
    case 'A'
        cfg.csvPaths = all_csv;
        cfg.mode_labels = all_lbls;
        fprintf('      -> Cross-mode: O1, O2, X1, X2\n');
    case 'C'
        fprintf('      Select modes (e.g. "1 3" for O1 + X1): ');
        sel_str = strtrim(input('', 's'));
        sel_idx = str2num(sel_str); %#ok<ST2NM>
        sel_idx = sel_idx(sel_idx >= 1 & sel_idx <= 4);
        if isempty(sel_idx), sel_idx = 1; end
        if numel(sel_idx) == 1
            cfg.csvPath = all_csv{sel_idx};
            fprintf('      -> %s\n', all_lbls{sel_idx});
        else
            cfg.csvPaths = all_csv(sel_idx);
            cfg.mode_labels = all_lbls(sel_idx);
            fprintf('      -> Cross-mode: %s\n', strjoin(all_lbls(sel_idx), ', '));
        end
    otherwise
        cfg.csvPath = all_csv{1}; fprintf('      -> O1\n');
end

%% [1] Operation Mode
fprintf('  [1] OPERATION MODE\n');
fprintf('      (A) ECCD — NTM suppression (Upper Launcher)\n');
fprintf('      (B) ECRH — Plasma heating (Equatorial Launcher)\n');
fprintf('      (C) Resonance Landscape Mapper (raw visualization)\n');
op = upper(strtrim(input('      [A/B/C] (default A): ', 's')));
if strcmp(op, 'B')
    cfg.operation_mode = 'ecrh';
    fprintf('      -> ECRH\n');
elseif strcmp(op, 'C')
    cfg.operation_mode = 'landscape';
    fprintf('      -> Landscape Mapper\n');
else
    cfg.operation_mode = 'eccd';
    fprintf('      -> ECCD\n');
end
is_eccd = strcmp(cfg.operation_mode, 'eccd');
is_landscape = strcmp(cfg.operation_mode, 'landscape');

%% Landscape: actuator triplets, then skip to summary
if is_landscape
    % Landscape requires single dataset
    if isfield(cfg, 'csvPaths')
        fprintf('      (Landscape uses single dataset — taking first: %s)\n', cfg.mode_labels{1});
        cfg.csvPath = cfg.csvPaths{1};
        cfg = rmfield(cfg, 'csvPaths');
        cfg = rmfield(cfg, 'mode_labels');
    end
    fprintf('\n  [2] ACTUATOR TRIPLETS (up to 4)\n');
    fprintf('      Enter {B0 [T], phi [deg], f [GHz]} per line.\n');
    fprintf('      Empty line to finish.\n');
    acts = [];
    for k = 1:6
        s = strtrim(input(sprintf('      [%d/6]: ', k), 's'));
        if isempty(s), break; end
        vals = str2num(s); %#ok<ST2NM>
        if numel(vals) == 3
            acts(end+1, :) = vals; %#ok<AGROW>
            fprintf('            -> B0=%.2f T, phi=%.1f°, f=%.0f GHz\n', vals(1), vals(2), vals(3));
        else
            fprintf('            (skipped — need exactly 3 values)\n');
        end
    end
    if isempty(acts)
        error('ATLAS:NoActuators', 'At least one actuator triplet required.');
    end
    cfg.landscape_actuators = acts;
    cfg.targeting_mode = 'landscape';
    cfg.ExportFigures = true;

    % Color axis toggle
    fprintf('\n  [3] COLOR AXIS\n');
    fprintf('      (1) n_{res} — resonant electron density [m^{-3}]\n');
    fprintf('      (2) f_{res} — resonant fraction n_{res}/n_r\n');
    cq_ls = strtrim(input('      [1/2] (default 1): ', 's'));
    if strcmp(cq_ls, '2')
        cfg.color_quantity = 'fres';
        fprintf('      -> f_res (resonant fraction)\n');
    else
        cfg.color_quantity = 'nres';
        fprintf('      -> n_res (resonant density)\n');
    end

    [~, ds_name] = fileparts(cfg.csvPath);
    fprintf('\n  +------------------------------------------------+\n');
    fprintf('  |  Dataset    : %-32s|\n', ds_name);
    fprintf('  |  Operation  : %-32s|\n', 'LANDSCAPE MAPPER');
    fprintf('  |  Actuators  : %-32s|\n', sprintf('%d triplet(s)', size(acts,1)));
    fprintf('  |  Color axis  : %-32s|\n', upper(cfg.color_quantity));
    fprintf('  +------------------------------------------------+\n\n');
    return;
end

%% [2] Targeting Mode (ECCD only)
if is_eccd
    fprintf('\n  [2] TARGETING MODE\n');
    fprintf('      (F) Forward — declare target, find best actuators\n');
    fprintf('      (R) Reverse — declare {B0,f}, find best flux surfaces\n');
    tm = upper(strtrim(input('      [F/R] (default F): ', 's')));
    if strcmp(tm, 'R')
        cfg.targeting_mode = 'reverse';
        fprintf('      -> Reverse\n');
    else
        cfg.targeting_mode = 'forward';
        fprintf('      -> Forward\n');
    end
else
    cfg.targeting_mode = 'ecrh';
end

is_forward = strcmp(cfg.targeting_mode, 'forward');
is_reverse = strcmp(cfg.targeting_mode, 'reverse');

%% [3] Target / Actuator inputs
if is_forward
    fprintf('\n  [3] TARGET\n');
    fprintf('      (P) Point target (r, theta)\n');
    fprintf('      (F) Flux surface\n');
    tt = upper(strtrim(input('      [P/F] (default F): ', 's')));
    if strcmp(tt, 'P')
        cfg.target_type = 'point';
        r_in = input(sprintf('      r [m] (default %.2f): ', cfg.r_user));
        if ~isempty(r_in) && isfinite(r_in) && r_in > 0, cfg.r_user = r_in; end
        th_in = input(sprintf('      theta [deg] (default %.1f): ', cfg.theta_user));
        if ~isempty(th_in) && isfinite(th_in), cfg.theta_user = th_in; end
        fprintf('      -> Point: r=%.2f m, theta=%.0f deg\n', cfg.r_user, cfg.theta_user);
    else
        cfg.target_type = 'flux';
        cfg.r_user = prompt_q_surface(cfg, true);
    end

    % Actuator constraints (with gate question)
    cfg = prompt_actuator_constraints(cfg);

elseif is_reverse
    fprintf('\n  [3] REVERSE: ACTUATOR SET\n');
    cfg.B0_reverse = parse_actuator_input('      B0 [T] (default 5.0)', 5.0);
    cfg.f_reverse = parse_actuator_input('      f [GHz] (default 170)', 170);

    fprintf('\n  [3b] Q-SURFACE FILTER (optional)\n');
    q_val = prompt_q_surface(cfg, false);
    if ~isnan(q_val), cfg.q_target = q_val; end

else  % ECRH
    fprintf('\n  [3] HEATING ZONE\n');
    ecrh_r_default = 0.40;
    r_in = input(sprintf('      r_max [m] (default %.2f): ', ecrh_r_default));
    if ~isempty(r_in) && isfinite(r_in) && r_in > 0
        cfg.r_user = r_in;
    else
        cfg.r_user = ecrh_r_default;
    end
    fprintf('      -> r <= %.2f m\n', cfg.r_user);

    % Actuator constraints (with gate question)
    cfg = prompt_actuator_constraints(cfg);
end

%% [4] Color Axis Quantity
fprintf('\n  [4] COLOR AXIS\n');
fprintf('      (1) n_{res} — resonant electron density [m^{-3}]\n');
fprintf('      (2) f_{res} — resonant fraction n_{res}/n_r\n');
cq = strtrim(input('      [1/2] (default 1): ', 's'));
if strcmp(cq, '2')
    cfg.color_quantity = 'fres';
    fprintf('      -> f_res (resonant fraction)\n');
else
    cfg.color_quantity = 'nres';
    fprintf('      -> n_res (resonant density)\n');
end

%% Summary
fprintf('\n  +------------------------------------------------+\n');
if isfield(cfg,'csvPaths')
    fprintf('  |  Modes      : %-32s|\n', strjoin(cfg.mode_labels, ', '));
else
    [~, fn_only] = fileparts(cfg.csvPath);
    fprintf('  |  Dataset    : %-32s|\n', fn_only);
end
fprintf('  |  Operation  : %-32s|\n', upper(cfg.operation_mode));
if is_eccd
    fprintf('  |  Targeting  : %-32s|\n', upper(cfg.targeting_mode));
end
if is_forward
    fprintf('  |  Target     : %-32s|\n', cfg.target_type);
end
fprintf('  |  eta         : %-32s|\n', sprintf('%.2f', cfg.eta));
fprintf('  |  Color axis  : %-32s|\n', upper(cfg.color_quantity));
if isfield(cfg,'q_target') && isfinite(cfg.q_target)
    fprintf('  |  q_target    : %-32s|\n', sprintf('%.3f', cfg.q_target));
end
fprintf('  +------------------------------------------------+\n\n');

end


%% ========== LOCAL HELPERS ==========

function val = parse_actuator_input(label, default_val)
    raw = strtrim(input(sprintf('%s: ', label), 's'));
    if isempty(raw), val = default_val; return; end
    cp = strsplit(raw, ':');
    if numel(cp)==2
        lo=str2double(cp{1}); hi=str2double(cp{2});
        if isfinite(lo)&&isfinite(hi)&&lo<=hi
            val=struct('range',[lo,hi]); return;
        end
    end
    try val=eval(raw);
        if isnumeric(val)&&~isempty(val)&&all(isfinite(val)), return; end
    catch
    end
    val = default_val;
end


function r_or_q = prompt_q_surface(cfg, is_required)
%PROMPT_Q_SURFACE  Shared q-surface selection menu.
%  is_required=true  -> must pick one (forward flux: sets r_user)
%  is_required=false -> can skip (reverse: sets q_target or NaN)

    fprintf('      Surface presets (by increasing r):\n');
    fprintf('        (1) q = 4/3  (r = %.2f m)\n', cfg.a*sqrt(max(0,(4/3-cfg.Q_c)/(cfg.Q_95-cfg.Q_c))));
    fprintf('        (2) q = 3/2  (r = %.2f m)\n', cfg.a*sqrt(max(0,(1.5-cfg.Q_c)/(cfg.Q_95-cfg.Q_c))));
    fprintf('        (3) q = 2    (r = %.2f m)\n', cfg.a*sqrt(max(0,(2-cfg.Q_c)/(cfg.Q_95-cfg.Q_c))));
    fprintf('        (4) Custom r [m]\n');
    if is_required
        qs = strtrim(input('      [1/2/3/4] (default 2): ', 's'));
    else
        fprintf('        (Enter) Skip — no filter\n');
        qs = strtrim(input('      [1/2/3/4/Enter]: ', 's'));
    end

    q_val = NaN; r_custom = NaN;
    switch qs
        case '1', q_val = 4/3;
        case '2', q_val = 1.5;
        case '3', q_val = 2.0;
        case '4'
            r_in = input('      r [m]: ');
            if ~isempty(r_in) && isfinite(r_in) && r_in > 0 && r_in < cfg.a
                r_custom = r_in;
                q_val = cfg.Q_c + (cfg.Q_95-cfg.Q_c)*(r_in/cfg.a)^2;
            end
        otherwise
            if is_required, q_val = 1.5; end  % safe default
    end

    if isnan(q_val) && isnan(r_custom)
        r_or_q = NaN;
        fprintf('      -> No surface filter\n');
    elseif ~isnan(r_custom)
        % Custom r: user gave r directly
        if is_required
            r_or_q = r_custom;
            fprintf('      -> r = %.4f m (q = %.3f)\n', r_custom, q_val);
        else
            r_or_q = q_val;
            fprintf('      -> q-filter: q = %.3f (r = %.4f m)\n', q_val, r_custom);
        end
    else
        frac = (q_val - cfg.Q_c) / (cfg.Q_95 - cfg.Q_c);
        if frac >= 0 && frac <= 1
            r_val = cfg.a * sqrt(frac);
            if is_required
                r_or_q = r_val;
                fprintf('      -> q = %.3f, r = %.4f m\n', q_val, r_val);
            else
                r_or_q = q_val;
                fprintf('      -> q-filter: q = %.3f (r = %.4f m)\n', q_val, r_val);
            end
        else
            fprintf('      -> q = %.3f outside profile range. Using q = 3/2.\n', q_val);
            frac = (1.5 - cfg.Q_c) / (cfg.Q_95 - cfg.Q_c);
            r_val = cfg.a * sqrt(max(0, frac));
            if is_required, r_or_q = r_val; else, r_or_q = 1.5; end
        end
    end
end


function cfg = prompt_actuator_constraints(cfg)
%PROMPT_ACTUATOR_CONSTRAINTS  Gate question + per-actuator input.
    gate = upper(strtrim(input('\n      Constrain actuators? [y/n] (default n): ', 's')));
    if strcmp(gate, 'Y')
        fprintf('      (Enter to skip any = scan all values)\n');
        cfg.B0_select = parse_actuator_input('        B0 [T]', NaN);
        cfg.phi_select = parse_actuator_input('        phi [deg]', NaN);
        cfg.f_select = parse_actuator_input('        f [GHz]', NaN);
    else
        cfg.B0_select = NaN; cfg.phi_select = NaN; cfg.f_select = NaN;
        fprintf('      -> No constraints (full actuator scan)\n');
    end
end


function cfg = default_cfg()
    cfg = struct();
    cfg.csvPath          = 'data_kinetic_O1.csv';
    cfg.operation_mode   = '';
    cfg.targeting_mode   = '';
    cfg.target_type      = 'flux';
    cfg.r_user           = 1.20;
    cfg.theta_user       = 45.0;
    cfg.q_target         = NaN;
    cfg.eta              = 0.80;
    cfg.threshold_percentile = 100;
    cfg.color_percentile = 100;

    % Actuator constraints
    cfg.B0_select = NaN; cfg.phi_select = NaN; cfg.f_select = NaN;
    cfg.B0_reverse = 5.0; cfg.f_reverse = 170;

    % Geometry
    cfg.a = 1.95; cfg.R_0 = 6.2; cfg.k_e = 1.7;
    cfg.Q_c = 1.0; cfg.Q_95 = 4.0;
    cfg.beta_p = 0.65; cfg.l_i = 0.85;

    % Kernel
    cfg.DeltaR = 0.06; cfg.sigma_theta = 2.0;
    cfg.sigma_theta_adaptive = true; cfg.kappa = 0.04;
    cfg.delta_theta_req = 5.0; cfg.theta_gap_floor = 5;
    cfg.tolRT_r = 0.01;

    % Launchers
    cfg.EL_R = 8.5; cfg.EL_Z = 0;
    cfg.UL_R = 7.1; cfg.UL_Z = 4.1;

    % LEAP
    cfg.beam_waist = 0.025; cfg.beam_freq_GHz = 170;
    cfg.pass_dr = 0.05; cfg.pass_search_r = 0.15;
    cfg.pass_exclude_deg = 5.0;
    cfg.phi_co_sign = -1;

    % Display
    cfg.RankIndices = [1,2,3];
    cfg.ShowBeam = true; cfg.ShowQSurfaces = true;
    cfg.ShowSeparatrix = true; cfg.ShowLegend = true;
    cfg.ShowAllPoints = true; cfg.ShowInfoBox = false;
    cfg.ShowGrid = false; cfg.BackgroundParula = false;
    cfg.FillPlasma = false; cfg.InvertColormap = true;
    cfg.MarkerSize = 40; cfg.MarkerSizeBG = 12;
    cfg.r_max_display = 1.90;
    cfg.min_points_for_Seff = 5;


    % Export
    cfg.ExportFigures = true; cfg.ExportDPI = 300;
    cfg.ExportFormat = 'png';
    cfg.color_quantity = 'nres';

    % Operational (informational)
end