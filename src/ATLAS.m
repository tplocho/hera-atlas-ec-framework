function results = ATLAS(cfg)
% ATLAS - Actuator-Target Landscape Analysis and Scoring  v1.0.7
%
%  Post-processing tool for EC resonance datasets produced by HERA.m.
%  Ranks actuator configurations {B0, phi, f} against user-specified targets
%  on the poloidal cross-section and produces publication-quality figures.
%
%  Usage:  results = ATLAS(cfg)              % scripted
%          cfg = ATLAS_int();                % interactive setup
%          results = ATLAS(cfg)
%
%  Three code paths:
%    ECCD Forward  — declare target, find best actuator configs
%    ECCD Reverse  — declare {B0,f}, find best flux surfaces + optimal phi
%    ECRH          — enclosed surface heating
%
%  Scoring: Score = ARETE × MOIRA (universal, all modes)
%    ARETE — Above-threshold Resonance Extent and Target Efficacy
%    MOIRA — Measure of Intercepted Resonance and Absorption
%
%  Ref: Bornatici et al., Nucl. Fusion 23, 1153 (1983);
%       Henderson et al., Nucl. Fusion 48, 054013 (2008).

VERSION = '1.0.7';
fprintf('\n  ATLAS v%s — Actuator-Target Landscape Analysis and Scoring\n\n', VERSION);

%% ================================================================
%%  [1] PARSE CFG, RESOLVE MODE
%% ================================================================

% --- Operation mode ---
if isfield(cfg,'operation_mode') && ~isempty(cfg.operation_mode)
    operation_mode = lower(cfg.operation_mode);
else
    error('ATLAS:NoMode', 'cfg.operation_mode must be ''eccd'', ''ecrh'', or ''landscape''.');
end
is_eccd = strcmp(operation_mode, 'eccd');
is_ecrh = strcmp(operation_mode, 'ecrh');
is_landscape = strcmp(operation_mode, 'landscape');
fprintf('  Operation: %s\n', upper(operation_mode));

% --- Landscape mode: early branch (no scan, no ranking) ---
if is_landscape
    results = run_landscape(cfg, VERSION);
    return;
end

% --- Targeting mode (ECCD only) ---
if is_eccd
    if isfield(cfg,'targeting_mode') && ~isempty(cfg.targeting_mode)
        targeting_mode = lower(cfg.targeting_mode);
    else
        error('ATLAS:NoTargeting', 'cfg.targeting_mode must be ''forward'' or ''reverse'' for ECCD.');
    end
    is_forward = strcmp(targeting_mode, 'forward');
    is_reverse = strcmp(targeting_mode, 'reverse');
    fprintf('  Targeting: %s\n', targeting_mode);
else
    targeting_mode = 'ecrh';
    is_forward = false;  is_reverse = false;
end

% --- Target type (forward only) ---
if is_forward
    if isfield(cfg,'target_type'), target_type = lower(cfg.target_type);
    else, target_type = 'flux'; end
    is_point = strcmp(target_type, 'point');
else
    target_type = 'flux';
    is_point = false;
end

% --- Compact operation label for filenames ---
op_label = compact_op_label(operation_mode, targeting_mode);

% --- Geometry ---
a = cfg.a;  R_0 = cfg.R_0;  k_e = cfg.k_e;
Q_c = cfg.Q_c;  Q_95 = cfg.Q_95;
if isfield(cfg,'beta_p'), beta_p=cfg.beta_p; else, beta_p=0.65; end
if isfield(cfg,'l_i'), l_i=cfg.l_i; else, l_i=0.85; end

delta_sh = @(r) (a^2./(2*R_0)).*(beta_p+l_i/2).*(1-(r./a).^2);
kappa_fn = @(r) 1+(k_e-1)*(r./a).^2;
Rfun = @(r,th) R_0 + r.*cos(th) + delta_sh(r);
Zfun = @(r,th) kappa_fn(r).*r.*sin(th);

% --- Target parameters ---
r_user = cfg.r_user;
if isfield(cfg,'theta_user'), theta_user = cfg.theta_user; else, theta_user = NaN; end
eta = cfg.eta;
q_target = Q_c + (Q_95-Q_c)*(r_user/a)^2;

% --- Q-surface filter ---
if isfield(cfg,'q_target') && isfinite(cfg.q_target)
    q_target_surface = cfg.q_target;
    frac_q = (q_target_surface-Q_c)/(Q_95-Q_c);
    if frac_q >= 0 && frac_q <= 1
        r_q_target = a*sqrt(frac_q);
        use_q_filter = true;
        fprintf('  Q-filter: q=%.3f -> r=%.3f m\n', q_target_surface, r_q_target);
    else
        use_q_filter = false;
        r_q_target = NaN;
    end
else
    q_target_surface = NaN; use_q_filter = false; r_q_target = NaN;
end

% --- Kernel (point mode) ---
DeltaR = cfg.DeltaR;
if isfield(cfg,'sigma_theta'), sigma_theta=cfg.sigma_theta; else, sigma_theta=2.0; end
if isfield(cfg,'sigma_theta_adaptive'), adaptive=cfg.sigma_theta_adaptive; else, adaptive=true; end
if isfield(cfg,'kappa'), kappa=cfg.kappa; else, kappa=0.04; end
if isfield(cfg,'delta_theta_req'), delta_theta_req=cfg.delta_theta_req; else, delta_theta_req=5.0; end
if isfield(cfg,'theta_gap_floor'), theta_gap_floor=cfg.theta_gap_floor; else, theta_gap_floor=5; end
tol_r = 0.01; if isfield(cfg,'tolRT_r'), tol_r=cfg.tolRT_r; end

% --- Launchers ---
if isfield(cfg,'EL_R'), R_EL=cfg.EL_R; else, R_EL=R_0+a+0.2; end
if isfield(cfg,'EL_Z'), Z_EL=cfg.EL_Z; else, Z_EL=0; end
if isfield(cfg,'UL_R'), R_UL=cfg.UL_R; else, R_UL=7.1; end
if isfield(cfg,'UL_Z'), Z_UL=cfg.UL_Z; else, Z_UL=4.1; end
if isfield(cfg,'EL_steer_range'), EL_steer=cfg.EL_steer_range; else, EL_steer=[-35,35]; end
if isfield(cfg,'UL_steer_range'), UL_steer=cfg.UL_steer_range; else, UL_steer=[-20,20]; end
EL_boresight = atan2d(0-Z_EL, R_EL-R_0);
UL_boresight = atan2d(0-Z_UL, R_UL-R_0);

if is_eccd
    R_active=R_UL; Z_active=Z_UL; active_bore=UL_boresight; active_steer=UL_steer;
    active_label='UL'; theta_L_sign=-1;
    if isfield(cfg,'phi_co_sign'), phi_co_sign=cfg.phi_co_sign; else, phi_co_sign=-1; end
else
    R_active=R_EL; Z_active=Z_EL; active_bore=EL_boresight; active_steer=EL_steer;
    active_label='EL'; theta_L_sign=1;
    phi_co_sign = 0;  % no co-current filter for ECRH
end
fprintf('  Launcher: %s at (%.1f, %.1f) m\n', active_label, R_active, Z_active);

% --- MOIRA parameters ---
if isfield(cfg,'beam_waist'), beam_w0=cfg.beam_waist; else, beam_w0=0.025; end
if isfield(cfg,'beam_freq_GHz'), beam_freq=cfg.beam_freq_GHz; else, beam_freq=170; end
beam_lambda = 299792458/(beam_freq*1e9);
beam_s_R = pi*beam_w0^2/beam_lambda;
if isfield(cfg,'pass_dr'), pass_dr=cfg.pass_dr; else, pass_dr=0.05; end
if isfield(cfg,'pass_search_r'), pass_search_r=cfg.pass_search_r; else, pass_search_r=0.15; end
if isfield(cfg,'pass_exclude_deg'), pass_exclude_deg=cfg.pass_exclude_deg; else, pass_exclude_deg=5.0; end

% --- Display ---
if isfield(cfg,'color_percentile'), color_pct=cfg.color_percentile; else, color_pct=95; end
if isfield(cfg,'threshold_percentile'), thr_pct=cfg.threshold_percentile; else, thr_pct=99; end
if isfield(cfg,'RankIndices'), RankIndices=cfg.RankIndices; else, RankIndices=[1,2,3]; end
if isfield(cfg,'r_max_display'), r_max_disp=cfg.r_max_display; else, r_max_disp=1.90; end
if isfield(cfg,'min_points_for_Seff'), min_pts=cfg.min_points_for_Seff; else, min_pts=5; end
if isfield(cfg,'color_quantity'), color_qty=lower(cfg.color_quantity); else, color_qty='nres'; end
use_fres = strcmp(color_qty, 'fres');

% --- Operational (informational) ---
if isfield(cfg,'B0_nom'), B0_nom=cfg.B0_nom; else, B0_nom=5.0; end
if isfield(cfg,'f_ref'), f_ref=cfg.f_ref; else, f_ref=170; end

%% ================================================================
%%  [2-5] LOAD, THRESHOLD, GRID, SCAN — per mode
%% ================================================================

% Resolve multi-mode vs single-mode
if isfield(cfg,'csvPaths') && iscell(cfg.csvPaths)
    csv_list = cfg.csvPaths;
    if isfield(cfg,'mode_labels'), mode_labels = cfg.mode_labels;
    else, mode_labels = arrayfun(@(i) sprintf('Mode%d',i), 1:numel(csv_list), 'Uni', false); end
    multi_mode = true;
else
    csv_list = {cfg.csvPath};
    mode_labels = {''};
    multi_mode = false;
end
n_modes = numel(csv_list);
if multi_mode
    fprintf('  Cross-mode: %d datasets [%s]\n', n_modes, strjoin(mode_labels, ', '));
end

% Storage for per-mode data (needed for plotting)
mode_data = struct('r',{},'theta',{},'R_geom',{},'Z_geom',{},'nres',{},'fres',{},'ImN2',{}, ...
    'B0',{},'phi',{},'f',{},'nres_ref',{},'nres_thr',{},'label',{});

S_all = {};  % cell array of per-mode scan results

for im = 1:n_modes
    csvPath = csv_list{im};
    mlabel = mode_labels{im};
    if multi_mode, fprintf('\n  --- Mode: %s ---\n', mlabel); end

    if ~isfile(csvPath), error('ATLAS:FileNotFound','File not found: %s', csvPath); end
    fprintf('  Loading: %s\n', csvPath);
    tic_load = tic;
    if isfield(cfg,'csv_robust') && cfg.csv_robust
        [T, n_skip] = robust_readtable(csvPath);
        fprintf('  %d rows in %.1f s', height(T), toc(tic_load));
        if n_skip>0, fprintf(' (%d malformed rows skipped)', n_skip); end
        fprintf('\n');
    else
        T = readtable(csvPath);
        fprintf('  %d rows in %.1f s\n', height(T), toc(tic_load));
    end

    % Validate columns
    req = ["r_m_","theta_deg_","R_m_","B0_T_","phi_deg_","f_GHz_","n_res"];
    missing = req(~ismember(req, string(T.Properties.VariableNames)));
    if ~isempty(missing), error('ATLAS:MissingCols','[%s] Missing: %s', mlabel, strjoin(missing,', ')); end

    r_m     = T.r_m_;
    theta_m = T.theta_deg_;
    R_m     = T.R_m_;
    B0_m    = T.B0_T_;
    phi_m   = T.phi_deg_;
    f_m     = T.f_GHz_;
    nres_m  = T.n_res;
    if ismember('ImN2', string(T.Properties.VariableNames))
        ImN2_m = T.ImN2;
    else
        ImN2_m = zeros(height(T),1);
    end
    if ismember('f_res', string(T.Properties.VariableNames))
        fres_m = T.f_res;
    else
        fres_m = zeros(height(T),1);
    end
    clear T;

    % Filter invalid rows
    gv = isfinite(R_m) & isfinite(nres_m) & isfinite(r_m) & isfinite(theta_m);
    if sum(~gv) > 0
        fprintf('  Removed %d invalid rows.\n', sum(~gv));
        r_m=r_m(gv); theta_m=theta_m(gv); R_m=R_m(gv);
        B0_m=B0_m(gv); phi_m=phi_m(gv); f_m=f_m(gv);
        nres_m=nres_m(gv); ImN2_m=ImN2_m(gv); fres_m=fres_m(gv);
    end
    r_cut = r_m <= r_max_disp;
    if sum(~r_cut) > 0
        fprintf('  Excluded %d points with r > %.2f m.\n', sum(~r_cut), r_max_disp);
        r_m=r_m(r_cut); theta_m=theta_m(r_cut); R_m=R_m(r_cut);
        B0_m=B0_m(r_cut); phi_m=phi_m(r_cut); f_m=f_m(r_cut);
        nres_m=nres_m(r_cut); ImN2_m=ImN2_m(r_cut); fres_m=fres_m(r_cut);
    end
    N_m = numel(r_m);
    R_geom_m = Rfun(r_m, deg2rad(theta_m));
    Z_geom_m = Zfun(r_m, deg2rad(theta_m));

    % Threshold (per-mode)
    nres_pos_m = nres_m(nres_m > 0);
    if isempty(nres_pos_m)
        fprintf('  [%s] No positive n_res — skipping.\n', mlabel); continue;
    end
    nres_ref_m = prctile(nres_pos_m, thr_pct);
    if nres_ref_m <= 0, nres_ref_m = max(nres_pos_m); end
    nres_thr_m = eta * nres_ref_m;
    fprintf('  nres_ref=%.2e (p%d), nres_thr=%.2e\n', nres_ref_m, thr_pct, nres_thr_m);

    % Store per-mode data for plotting
    md = struct('r',r_m,'theta',theta_m,'R_geom',R_geom_m,'Z_geom',Z_geom_m, ...
        'nres',nres_m,'fres',fres_m,'ImN2',ImN2_m,'B0',B0_m,'phi',phi_m,'f',f_m, ...
        'nres_ref',nres_ref_m,'nres_thr',nres_thr_m,'label',mlabel);
    mode_data(im) = md;

    % Actuator grid (per-mode — different datasets may have different grids)
    B0u_m = unique(B0_m);  phiu_m = unique(phi_m);  fu_m = unique(f_m);

    if is_forward || is_ecrh
        if isfield(cfg,'B0_select')
            [B0s,~,~]=resolve_param(cfg.B0_select,B0u_m,'B0','T');
        else
            B0s=B0u_m;
        end
        if isfield(cfg,'phi_select') && ~(isscalar(cfg.phi_select) && isnan(cfg.phi_select))
            [phis,~,~]=resolve_param(cfg.phi_select,phiu_m,'phi','deg');
            phi_user_set = true;
        else
            phis=phiu_m; phi_user_set=false;
        end
        if isfield(cfg,'f_select')
            [fs,~,~]=resolve_param(cfg.f_select,fu_m,'f','GHz');
        else
            fs=fu_m;
        end
    elseif is_reverse
        if isfield(cfg,'B0_reverse')
            [B0s,~,~]=resolve_param(cfg.B0_reverse,B0u_m,'B0','T');
        else
            error('ATLAS:NoB0','cfg.B0_reverse required.');
        end
        if isfield(cfg,'f_reverse')
            [fs,~,~]=resolve_param(cfg.f_reverse,fu_m,'f','GHz');
        else
            error('ATLAS:NoF','cfg.f_reverse required.');
        end
        phis=phiu_m; phi_user_set=false;
    end

    nBm=numel(B0s); nPm=numel(phis); nFm=numel(fs);
    B0_masks_m=false(N_m,nBm); for i=1:nBm, B0_masks_m(:,i)=abs(B0_m-B0s(i))<1e-6; end
    phi_masks_m=false(N_m,nPm); for i=1:nPm, phi_masks_m(:,i)=abs(phi_m-phis(i))<1e-6; end
    f_masks_m=false(N_m,nFm);   for i=1:nFm, f_masks_m(:,i)=abs(f_m-fs(i))<1e-6; end

    fprintf('  Actuator grid: %d B0 x %d phi x %d f = %d combos\n', nBm, nPm, nFm, nBm*nPm*nFm);

    % Scan
    tic_scan = tic;
    if is_forward
        if is_point
            fprintf('  Forward scan: point (r=%.2f, theta=%.0f)\n', r_user, theta_user);
        else
            fprintf('  Forward scan: flux (r=%.2f, q=%.2f)\n', r_user, q_target);
        end
        Sm = scan_forward(r_m, theta_m, R_geom_m, Z_geom_m, nres_m, ImN2_m, ...
            B0s, phis, fs, B0_masks_m, phi_masks_m, f_masks_m, ...
            r_user, theta_user, is_point, q_target, ...
            nres_thr_m, nres_ref_m, eta, DeltaR, sigma_theta, adaptive, kappa, ...
            delta_theta_req, theta_gap_floor, tol_r, min_pts, ...
            R_active, Z_active, active_bore, active_steer, Rfun, Zfun, ...
            a, R_0, phi_co_sign, phi_user_set, ...
            pass_dr, pass_search_r, pass_exclude_deg, beam_w0, beam_s_R, N_m);
    elseif is_reverse
        fprintf('  Reverse scan: {B0,f} -> phi\n');
        Sm = scan_reverse(r_m, theta_m, R_geom_m, Z_geom_m, nres_m, ImN2_m, ...
            B0s, phis, fs, B0_masks_m, phi_masks_m, f_masks_m, ...
            nres_thr_m, nres_ref_m, theta_gap_floor, tol_r, ...
            R_active, Z_active, active_bore, active_steer, Rfun, Zfun, ...
            a, R_0, Q_c, Q_95, phi_co_sign, sigma_theta, adaptive, kappa, ...
            pass_dr, pass_search_r, pass_exclude_deg, beam_w0, beam_s_R, N_m);
    else
        fprintf('  ECRH scan: r<=%.2f\n', r_user);
        Sm = scan_ecrh(r_m, theta_m, R_geom_m, Z_geom_m, nres_m, ImN2_m, ...
            B0s, phis, fs, B0_masks_m, phi_masks_m, f_masks_m, ...
            r_user, nres_thr_m, nres_ref_m, theta_gap_floor, tol_r, min_pts, ...
            R_active, Z_active, Rfun, Zfun, a, R_0, ...
            pass_dr, pass_search_r, pass_exclude_deg, beam_w0, beam_s_R, N_m);
    end
    fprintf('  Scan: %.2f s, %d results\n', toc(tic_scan), Sm.nv);

    % Tag results with mode index
    if Sm.nv > 0
        Sm.mode_idx = repmat(im, Sm.nv, 1);
        Sm.mode_label = repmat({mlabel}, Sm.nv, 1);
        S_all{end+1} = Sm; %#ok<AGROW>
    end
end

% --- Merge all mode results into one S struct ---
if isempty(S_all)
    S = alloc_results(0, is_reverse); S.nv = 0;
    S.mode_idx = []; S.mode_label = {};
else
    fn_merge = fieldnames(S_all{1});
    S = struct();
    for fi = 1:numel(fn_merge)
        fname = fn_merge{fi};
        if strcmp(fname, 'nv'), continue; end
        parts = cellfun(@(s) s.(fname), S_all, 'Uni', false);
        S.(fname) = vertcat(parts{:});
    end
    S.nv = sum(cellfun(@(s) s.nv, S_all));
end

% Cross-mode diagnostics
if multi_mode && S.nv > 0
    fprintf('\n  --- Cross-mode merge diagnostics ---\n');
    for im = 1:n_modes
        mlbl = mode_labels{im};
        midx = find(cellfun(@(x) strcmp(x, mlbl), S.mode_label));
        if isempty(midx)
            fprintf('    [%s] 0 configs\n', mlbl);
        else
            arete_m = S.ARETE(midx);
            fprintf('    [%s] %d configs, max ARETE_raw=%.2e, non-zero=%d\n', ...
                mlbl, numel(midx), max(arete_m), sum(arete_m > 0));
        end
    end
    fprintf('  Total merged: %d configs\n', S.nv);
end

% Set workspace aliases from first mode (used by fallback + single-mode plotting)
if ~isempty(mode_data)
    r_all = mode_data(1).r; theta_all = mode_data(1).theta;
    R_geom_all = mode_data(1).R_geom; Z_geom_all = mode_data(1).Z_geom;
    nres_all = mode_data(1).nres; fres_all = mode_data(1).fres;
    ImN2_all = mode_data(1).ImN2;
    B0_all = mode_data(1).B0; phi_all = mode_data(1).phi; f_all = mode_data(1).f;
    nres_ref = mode_data(1).nres_ref; nres_thr = mode_data(1).nres_thr;
    N_data = numel(r_all);
    B0u = unique(B0_all); phiu = unique(phi_all); fu = unique(f_all);
end

%% ================================================================
%%  [6] SCORE = ARETE x MOIRA, SORT, DEDUP
%% ================================================================
if S.nv == 0
    % --- FALLBACK: no results ---
    fprintf('\n  *** No viable configurations found. ***\n');
    results = run_fallback(cfg, operation_mode, targeting_mode, is_point, ...
        r_user, theta_user, q_target, use_q_filter, q_target_surface, r_q_target, ...
        r_all, theta_all, R_geom_all, Z_geom_all, nres_all, ImN2_all, ...
        B0_all, phi_all, f_all, B0u, phiu, fu, ...
        nres_thr, nres_ref, eta, DeltaR, sigma_theta, adaptive, kappa, ...
        delta_theta_req, theta_gap_floor, tol_r, min_pts, ...
        R_active, Z_active, active_bore, active_steer, Rfun, Zfun, ...
        a, R_0, Q_c, Q_95, phi_co_sign, ...
        pass_dr, pass_search_r, pass_exclude_deg, beam_w0, beam_s_R, ...
        B0_nom, f_ref, N_data, VERSION);
    return;
end

% Trim pre-allocated arrays
fn = fieldnames(S);
nv_trim = S.nv;
for i = 1:numel(fn)
    v = S.(fn{i});
    if (isnumeric(v) || iscell(v)) && numel(v) > nv_trim
        S.(fn{i}) = v(1:nv_trim);
    end
end

keep = S.ARETE > 0;
n_zero = S.nv - sum(keep);
if n_zero > 0
    fprintf('  Filtered %d configs with ARETE = 0.\n', n_zero);
    fn_all = fieldnames(S);
    for i = 1:numel(fn_all)
        v = S.(fn_all{i});
        if (isnumeric(v) || iscell(v)) && numel(v) == S.nv
            S.(fn_all{i}) = v(keep);
        end
    end
    S.nv = sum(keep);
end

% Filter MOIRA >= 1.0 (unphysical: zero beam-path interference = aiming at empty space)
keep_moira = S.MOIRA < (1.0 - 1e-6);
n_moira1 = S.nv - sum(keep_moira);
if n_moira1 > 0
    fprintf('  Filtered %d configs with MOIRA >= 1.0 (unphysical).\n', n_moira1);
    fn_all = fieldnames(S);
    for i = 1:numel(fn_all)
        v = S.(fn_all{i});
        if (isnumeric(v) || iscell(v)) && numel(v) == S.nv
            S.(fn_all{i}) = v(keep_moira);
        end
    end
    S.nv = sum(keep_moira);
end

if S.nv == 0
    fprintf('\n  *** No viable configs (ARETE=0 or MOIRA>=1.0 for all). Running fallback. ***\n');
    results = run_fallback(cfg, operation_mode, targeting_mode, is_point, ...
        r_user, theta_user, q_target, use_q_filter, q_target_surface, r_q_target, ...
        r_all, theta_all, R_geom_all, Z_geom_all, nres_all, ImN2_all, ...
        B0_all, phi_all, f_all, B0u, phiu, fu, ...
        nres_thr, nres_ref, eta, DeltaR, sigma_theta, adaptive, kappa, ...
        delta_theta_req, theta_gap_floor, tol_r, min_pts, ...
        R_active, Z_active, active_bore, active_steer, Rfun, Zfun, ...
        a, R_0, Q_c, Q_95, phi_co_sign, ...
        pass_dr, pass_search_r, pass_exclude_deg, beam_w0, beam_s_R, ...
        B0_nom, f_ref, N_data, VERSION);
    return;
end

% Apply q-filter (reverse mode) — BEFORE normalization so that ARETE is
% normalized within the target-surface population, not across all surfaces.
if is_reverse && use_q_filter
    q_keep = abs(S.r_opt - r_q_target) <= tol_r;
    n_pass = sum(q_keep);
    fprintf('  Q-filter: %d/%d pass (|r - r_q| <= %.3f m)\n', n_pass, S.nv, tol_r);
    if n_pass == 0
        fprintf('  *** No surfaces match q=%.2f (r=%.3f m). Running fallback. ***\n', q_target_surface, r_q_target);
        results = run_fallback(cfg, operation_mode, targeting_mode, is_point, ...
            r_user, theta_user, q_target, use_q_filter, q_target_surface, r_q_target, ...
            r_all, theta_all, R_geom_all, Z_geom_all, nres_all, ImN2_all, ...
            B0_all, phi_all, f_all, B0u, phiu, fu, ...
            nres_thr, nres_ref, eta, DeltaR, sigma_theta, adaptive, kappa, ...
            delta_theta_req, theta_gap_floor, tol_r, min_pts, ...
            R_active, Z_active, active_bore, active_steer, Rfun, Zfun, ...
            a, R_0, Q_c, Q_95, phi_co_sign, ...
            pass_dr, pass_search_r, pass_exclude_deg, beam_w0, beam_s_R, ...
            B0_nom, f_ref, N_data, VERSION);
        return;
    end
    fn_all = fieldnames(S);
    for i = 1:numel(fn_all)
        v = S.(fn_all{i});
        if (isnumeric(v) || iscell(v)) && numel(v) == S.nv
            S.(fn_all{i}) = v(q_keep);
        end
    end
    S.nv = n_pass;
end

% Normalize ARETE to (0, 1]: best config = 1.0, others scale proportionally
ARETE_max = max(S.ARETE);
if ARETE_max > 0
    S.ARETE = S.ARETE / ARETE_max;
    fprintf('  ARETE normalized: max_raw=%.2e -> [%.3f, 1.000]\n', ARETE_max, min(S.ARETE));
end

% Universal scoring
score = S.ARETE .* S.MOIRA;

% Sort descending
[~, si] = sort(score, 'descend');
fn_sort = {'B0','phi','f','r_opt','theta_opt','R_opt','Z_opt', ...
           'nres_peak','N_above','N_total','ARETE','MOIRA','theta_ranges', ...
           'coverage_deg','sigma_theta_eff','mode_idx','mode_label'};
if is_reverse
    fn_sort = [fn_sort, {'q_at_target'}];
end
for i = 1:numel(fn_sort)
    if isfield(S, fn_sort{i})
        S.(fn_sort{i}) = S.(fn_sort{i})(si);
    end
end
S.score = score(si);

% Determine dedup strategy for the ranking table.
%   Reverse + q-filter : phi is sole free variable -> expand all phi
%   Reverse, no q      : best phi per surface -> dedup by {B0, f, r_opt}
%   Forward/ECRH 1 set : phi is sole free variable -> expand all phi
%   Forward/ECRH N sets: best phi per actuator set -> dedup by {B0, f}
Bf_tuples = unique(round([S.B0, S.f], 4), 'rows');
n_Bf = size(Bf_tuples, 1);

n_full = S.nv;
true_rank = (1:n_full)';

% Store full ranking for CSV before dedup
S_full = S;
S_full.true_rank = true_rank;

% Build dedup key — S is already sorted by score descending, so
% unique(...,'stable') retains the best-scoring row per key.
if is_reverse && ~use_q_filter
    % Reverse without q: best phi per flux surface
    if multi_mode
        dedup_key = [round(S.B0,4), round(S.f,2), round(S.r_opt,4), S.mode_idx];
    else
        dedup_key = [round(S.B0,4), round(S.f,2), round(S.r_opt,4)];
    end
    fprintf('  Table mode: best phi per surface (%d unique surfaces)\n', ...
        size(unique(round(S.r_opt,4)), 1));
elseif is_reverse || (n_Bf == 1)
    % Reverse with q / Forward-ECRH single {B0,f}: expand all phi
    if multi_mode
        dedup_key = [round(S.B0,4), round(S.phi,2), round(S.f,2), S.mode_idx];
    else
        dedup_key = [round(S.B0,4), round(S.phi,2), round(S.f,2)];
    end
    fprintf('  Table mode: full phi expansion (phi is sole free variable)\n');
else
    % Forward/ECRH with multiple {B0,f}: collapse phi
    if multi_mode
        dedup_key = [round(S.B0,4), round(S.f,2), S.mode_idx];
    else
        dedup_key = [round(S.B0,4), round(S.f,2)];
    end
    fprintf('  Table mode: best phi per {B0, f} set (%d unique tuples)\n', n_Bf);
end
[~, dedup_idx] = unique(dedup_key, 'rows', 'stable');

% Apply dedup for display
for i = 1:numel(fn_sort)
    if isfield(S, fn_sort{i})
        S.(fn_sort{i}) = S.(fn_sort{i})(dedup_idx);
    end
end
S.score = S.score(dedup_idx);
S.true_rank = true_rank(dedup_idx);
S.nv = numel(dedup_idx);

fprintf('  Ranking: %d total, %d after dedup\n', n_full, S.nv);

%% ================================================================
%%  [7] CSV EXPORT
%% ================================================================
ts = datestr(now, 'yyyymmdd_HHMMSS');

% Build mode tag for filenames
if multi_mode
    mode_tag = strjoin(mode_labels, '-');
else
    [~, fn_csv] = fileparts(cfg.csvPath);
    % Extract mode label from filename (e.g. 'data_kinetic_O1' -> 'O1')
    tok = regexp(fn_csv, '(O[12]|X[12])$', 'match');
    if ~isempty(tok), mode_tag = tok{1};
    else, mode_tag = 'single'; end
end

csv_fn = sprintf('ATLAS_%s_%s_%s.csv', op_label, mode_tag, ts);
fid = fopen(csv_fn, 'w');
if fid > 0
    % Build header
    hdr = 'Rank';
    if multi_mode, hdr = [hdr, ',Mode']; end
    if is_reverse
        hdr = [hdr, ',r_m,q,B0_T,phi_deg,f_GHz,theta_opt_deg,theta_L_deg,ARETE,MOIRA,Score,nres_peak,N_above,N_total'];
    else
        hdr = [hdr, ',B0_T,phi_deg,f_GHz,theta_opt_deg,theta_L_deg,ARETE,MOIRA,Score,nres_peak,N_above,N_total'];
    end
    fprintf(fid, '%s\n', hdr);

    for ic = 1:n_full
        thL = theta_L_sign*(atan2d(S_full.Z_opt(ic)-Z_active, R_active-S_full.R_opt(ic))-active_bore);

        % Rank + optional mode
        fprintf(fid, '%d', ic);
        if multi_mode, fprintf(fid, ',%s', S_full.mode_label{ic}); end

        if is_reverse
            q_ic = Q_c + (Q_95-Q_c)*(S_full.r_opt(ic)/a)^2;
            fprintf(fid, ',%.4f,%.3f', S_full.r_opt(ic), q_ic);
        end

        fprintf(fid, ',%.2f,%.1f,%.1f,%.1f,%.1f,%.4f,%.4f,%.4f,%.4e,%d,%d\n', ...
            S_full.B0(ic), S_full.phi(ic), S_full.f(ic), ...
            S_full.theta_opt(ic), thL, S_full.ARETE(ic), S_full.MOIRA(ic), S_full.score(ic), ...
            S_full.nres_peak(ic), S_full.N_above(ic), S_full.N_total(ic));
    end
    fclose(fid);
    fprintf('  CSV: %s (%d rows)\n', csv_fn, n_full);
end

%% ================================================================
%%  [8] PUBLICATION-QUALITY PLOTTING
%% ================================================================

% --- Geometry overlays ---
thg = linspace(0, 2*pi, 361);
R_sep = Rfun(r_max_disp + 0.05, thg);  Z_sep = Zfun(r_max_disp + 0.05, thg);
R_tgt_c = Rfun(r_user, thg);  Z_tgt_c = Zfun(r_user, thg);

qv = [4/3, 3/2, 2]; ql = {'q = 4/3','q = 3/2','q = 2'};
qc = {[.4 .4 .6],[.7 .15 .15],[.3 .6 .3]};
qr = NaN(size(qv)); qRc = cell(size(qv)); qZc = cell(size(qv));
for iq = 1:numel(qv)
    fr = (qv(iq)-Q_c)/(Q_95-Q_c);
    if fr >= 0 && fr <= 1
        qr(iq) = a*sqrt(fr);
        qRc{iq} = Rfun(qr(iq),thg); qZc{iq} = Zfun(qr(iq),thg);
    end
end

% --- Pre-compute actuator keys per mode for fast plotting extraction ---
if multi_mode
    act_keys_per_mode = cell(n_modes, 1);
    for im = 1:n_modes
        md = mode_data(im);
        act_keys_per_mode{im} = round(md.B0,6)*1e12 + round(md.phi,6)*1e6 + round(md.f,4);
    end
else
    act_keys_all = round(B0_all,6)*1e12 + round(phi_all,6)*1e6 + round(f_all,4);
end

% --- Select configs to plot ---
if multi_mode
    % Top 1 from each mode that survived scoring.
    % Modes with negligible ARETE (below floor) are skipped with a console warning.
    arete_floor = 0.01;  % minimum normalized ARETE to qualify for plot
    ri_list = [];
    modes_seen = {};
    for k = 1:S.nv
        mlbl = S.mode_label{k};
        if ~ismember(mlbl, modes_seen)
            if S.ARETE(k) >= arete_floor
                ri_list(end+1) = k; %#ok<AGROW>
                modes_seen{end+1} = mlbl; %#ok<AGROW>
            else
                % Mode's best config is below floor — skip with warning
                modes_seen{end+1} = mlbl; %#ok<AGROW>
                fprintf('  [%s] Skipped — ARETE below threshold (best normalized ARETE = %.2e).\n', mlbl, S.ARETE(k));
                fprintf('         The current actuator selection is not well-supported by this propagation mode.\n');
            end
        end
        if numel(modes_seen) >= n_modes, break; end
    end
    % Modes eliminated entirely by zero-ARETE filter (no surviving configs at all)
    for im = 1:n_modes
        if ~ismember(mode_labels{im}, modes_seen)
            fprintf('  [%s] Skipped — no viable configurations (ARETE = 0 for all actuators).\n', mode_labels{im});
            fprintf('         The current actuator selection is not supported by this propagation mode.\n');
        end
    end
    fprintf('  Multi-mode: %d mode(s) with viable results plotted.\n', numel(ri_list));
else
    ri_list = RankIndices(RankIndices <= S.nv);
end
if isempty(ri_list), ri_list = 1; end
n_plots = numel(ri_list);
figs = gobjects(n_plots, 1);

for ip = 1:n_plots
    ri = ri_list(ip);
    act_B0 = S.B0(ri);  act_phi = S.phi(ri);  act_f = S.f(ri);
    act_ARETE = S.ARETE(ri); act_MOIRA = S.MOIRA(ri);
    act_score = S.score(ri); act_rank = S.true_rank(ri);
    R_tgt = S.R_opt(ri);  Z_tgt = S.Z_opt(ri);
    th_opt = S.theta_opt(ri);  r_opt = S.r_opt(ri);
    if multi_mode
        act_mode_idx = S.mode_idx(ri);
        act_mode_lbl = S.mode_label{ri};
    else
        act_mode_idx = 1;
        act_mode_lbl = '';
    end

    if multi_mode
        fprintf('  [Plot %d/%d] Rank #%d [%s]: B0=%.2f, phi=%.1f, f=%.1f\n', ...
            ip, n_plots, act_rank, act_mode_lbl, act_B0, act_phi, act_f);
    else
        fprintf('  [Plot %d/%d] Rank #%d: B0=%.2f, phi=%.1f, f=%.1f\n', ...
            ip, n_plots, act_rank, act_B0, act_phi, act_f);
    end

    % Extract actuator data from the correct mode dataset
    if multi_mode
        md = mode_data(act_mode_idx);
        act_key = round(act_B0,6)*1e12 + round(act_phi,6)*1e6 + round(act_f,4);
        mask_act = (act_keys_per_mode{act_mode_idx} == act_key);
        nres_act = md.nres(mask_act);  fres_act = md.fres(mask_act);
        Rg = md.R_geom(mask_act);
        Zd = md.Z_geom(mask_act); rd_act = md.r(mask_act);
    else
        act_key = round(act_B0,6)*1e12 + round(act_phi,6)*1e6 + round(act_f,4);
        mask_act = (act_keys_all == act_key);
        nres_act = nres_all(mask_act);  fres_act = fres_all(mask_act);
        Rg = R_geom_all(mask_act);
        Zd = Z_geom_all(mask_act); rd_act = r_all(mask_act);
    end

    % Resolve color quantity
    if use_fres
        cdata_act = fres_act;
        cbar_label = 'f_{res} = n_{res}/n_r';
    else
        cdata_act = nres_act;
        cbar_label = 'n_{res} [m^{-3}]';
    end

    % Color scale: robust per-actuator percentile
    nres_act_pos = nres_act(isfinite(nres_act) & nres_act > 0);
    cdata_pos = cdata_act(isfinite(cdata_act) & cdata_act > 0);
    if use_fres
        cax_top = min(1.0, prctile(cdata_pos, color_pct));
        if isempty(cdata_pos) || cax_top <= 0, cax_top = 1.0; end
    else
        if isempty(cdata_pos), cax_top = 1;
        else, cax_top = prctile(cdata_pos, color_pct); if cax_top <= 0, cax_top = max(cdata_pos); end; end
    end

    % --- Figure ---
    figs(ip) = figure('Color','w', 'Position',[80+30*ip, 60+30*ip, 1100, 850]);
    hold on; axis equal;

    % Colormap selection
    %cmap = parula(256);
    cmap = parula_white(256, 0.15);

    if isfield(cfg,'InvertColormap') && cfg.InvertColormap
        cmap = flipud(cmap);
    end

    % Plasma interior fill (drawn first = behind everything)
    if isfield(cfg,'FillPlasma') && cfg.FillPlasma
        patch(R_sep, Z_sep, cmap(1,:), 'EdgeColor', 'none', 'HandleVisibility', 'off');
    end

    % Separatrix
    if isfield(cfg,'ShowSeparatrix') && cfg.ShowSeparatrix
        plot(R_sep, Z_sep, 'k-', 'LineWidth', 1.8, 'HandleVisibility','off');
    end

    % Scatter: mode-aware rendering (drawn BEFORE overlays so surfaces appear on top)
    % rd_act already set in data extraction above
    ms_scatter = cfg.MarkerSize * 0.4;
    if is_ecrh && ~isempty(nres_act_pos)
        % ECRH: grayscale background, parula only inside enclosed surface
        in_zone = (rd_act <= r_user);
        out_zone = ~in_zone;

        % Outside: grayscale (white for low, dark gray for high — matches parula_white)
        if any(out_zone)
            gvals = max(0.30, min(1.0 - 0.70*(cdata_act(out_zone)/cax_top), 1.0));
            Rg_out = Rg(out_zone); Zd_out = Zd(out_zone);
            [~, sib] = sortrows([gvals(:), Rg_out(:), Zd_out(:)], [-1 2 3]);
            scatter(Rg_out(sib), Zd_out(sib), ms_scatter, repmat(gvals(sib), 1, 3), 'filled', ...
                'MarkerFaceAlpha', 0.35, 'MarkerEdgeColor','none', 'HandleVisibility','off');
        end

        % Inside: parula color (in-footprint alpha)
        if any(in_zone)
            cdata_in = min(max(cdata_act(in_zone),0), cax_top);
            Rg_in = Rg(in_zone); Zd_in = Zd(in_zone);
            [~, siw] = sortrows([cdata_in(:), Rg_in(:), Zd_in(:)], [1 2 3]);
            scatter(Rg_in(siw), Zd_in(siw), ms_scatter, cdata_in(siw), 'filled', ...
                'MarkerFaceAlpha', 0.95, 'MarkerEdgeColor','none', 'HandleVisibility','off');
        end
    elseif ~isempty(nres_act_pos)
        % ECCD: full parula scatter
        cdata_clipped = min(max(cdata_act,0), cax_top);
        [~, siw] = sortrows([cdata_clipped(:), Rg(:), Zd(:)], [1 2 3]);
        scatter(Rg(siw), Zd(siw), ms_scatter, cdata_clipped(siw), 'filled', ...
            'MarkerFaceAlpha', 0.85, 'MarkerEdgeColor', 'none', 'HandleVisibility','off');
    end

    % Q-surfaces (drawn on top of scatter for visibility)
    if isfield(cfg,'ShowQSurfaces') && cfg.ShowQSurfaces
        for iq = 1:numel(qv)
            if isnan(qr(iq)), continue; end
            plot(qRc{iq}, qZc{iq}, '--', 'Color', qc{iq}, 'LineWidth', 1.5, ...
                'DisplayName', sprintf('%s (r = %.2f m)', ql{iq}, qr(iq)));
        end
    end

    % Target overlay (drawn on top of scatter for visibility)
    if is_reverse
        R_opt_c = Rfun(r_opt, thg); Z_opt_c = Zfun(r_opt, thg);
        q_opt = Q_c + (Q_95-Q_c)*(r_opt/a)^2;
        plot(R_opt_c, Z_opt_c, '-', 'Color', [1 .65 0], 'LineWidth', 2.5, ...
            'DisplayName', sprintf('q = %s (r = %.2f m)', q_rational(q_opt), r_opt));
    elseif is_ecrh
        patch(R_tgt_c, Z_tgt_c, [1 .65 0], 'FaceAlpha', 0.10, ...
            'EdgeColor', 'none', 'HandleVisibility', 'off');
        plot(R_tgt_c, Z_tgt_c, '-', 'Color', [1 .65 0], 'LineWidth', 3.0, ...
            'DisplayName', sprintf('r \\leq %.2f m', r_user));
    else
        if is_point
            % No surface ring for point mode
        else
            plot(R_tgt_c, Z_tgt_c, '-', 'Color', [1 .65 0], 'LineWidth', 2.5, ...
                'DisplayName', sprintf('q = %s surface', q_rational(q_target)));
        end
    end

    % Colorbar
    colormap(cmap);
    cb = colorbar; cb.Label.String = cbar_label;
    cb.Label.FontSize = 12;
    caxis([0, cax_top]);

    % Target marker
    thL_opt = theta_L_sign*(atan2d(Z_tgt-Z_active, R_active-R_tgt)-active_bore);
    % Visible marker on the actual plot (small, HandleVisibility off)
    plot(R_tgt, Z_tgt, 'p', 'MarkerSize', 11, 'LineWidth', 2, ...
        'Color', [1 1 1], 'MarkerFaceColor', [.85 .10 .10], ...
        'HandleVisibility', 'off');
    % Legend-only entry: NaN marker carries the DisplayName
    plot(NaN, NaN, 'p', 'MarkerSize', 14, 'LineWidth', 2, ...
        'Color', [1 1 1], 'MarkerFaceColor', [.85 .10 .10], ...
        'DisplayName', sprintf('Target (\\theta_L = %.1f°)', thL_opt));

    % Magnetic axis
    plot(R_0, 0, '+', 'MarkerSize', 12, 'LineWidth', 1.5, 'Color', [.4 .4 .4], 'HandleVisibility','off');

    % === BEAM + LAUNCHERS ===
    if isfield(cfg,'ShowBeam') && cfg.ShowBeam
        bcol = [.85 .10 .10];

        blen = hypot(R_active-R_tgt, Z_active-Z_tgt);
        if blen > 0.1
            bd = [R_tgt-R_active, Z_tgt-Z_active]/blen;
            bp = [-bd(2), bd(1)];  % perpendicular direction

            % --- Beam path sampling ---
            n_env = 60;
            s_env = linspace(0, blen, n_env);

            R_env = R_active + bd(1)*s_env;
            Z_env = Z_active + bd(2)*s_env;

            % --- Beam footprint: cone + deposition at same opacity ---
            if ~is_ecrh
                fp_r = r_opt;  fp_th = th_opt;
                sig_eff_plot = S.sigma_theta_eff(ri);
                if sig_eff_plot > 0
                    th_arc = linspace(fp_th - 2*sig_eff_plot, fp_th + 2*sig_eff_plot, 51);
                    th_arc_rad = deg2rad(th_arc);
                    r_in = max(fp_r - DeltaR, 0.001);  r_out = fp_r + DeltaR;
                    R_dep_in  = Rfun(r_in,  th_arc_rad);  Z_dep_in  = Zfun(r_in,  th_arc_rad);
                    R_dep_out = Rfun(r_out, th_arc_rad);  Z_dep_out = Zfun(r_out, th_arc_rad);

                    % Find outermost arc corners via beam-perpendicular projection
                    th_lo_rad = deg2rad(fp_th - 2*sig_eff_plot);
                    th_hi_rad = deg2rad(fp_th + 2*sig_eff_plot);
                    Rc = [Rfun(r_out,th_hi_rad), Rfun(r_out,th_lo_rad), ...
                          Rfun(r_in, th_hi_rad), Rfun(r_in, th_lo_rad)];
                    Zc = [Zfun(r_out,th_hi_rad), Zfun(r_out,th_lo_rad), ...
                          Zfun(r_in, th_hi_rad), Zfun(r_in, th_lo_rad)];
                    proj = (Rc - R_active) .* bp(1) + (Zc - Z_active) .* bp(2);
                    [~, idx_L] = max(proj);  [~, idx_R] = min(proj);

                    fp_alpha = 0.25;

                    % Cone: launcher -> outermost corners
                    patch([R_active, Rc(idx_L), Rc(idx_R)], ...
                          [Z_active, Zc(idx_L), Zc(idx_R)], ...
                          bcol, 'FaceAlpha', fp_alpha, 'EdgeColor','none', 'HandleVisibility','off');

                    % Deposition arc (doubled alpha for visibility)
                    patch([R_dep_out, fliplr(R_dep_in)], [Z_dep_out, fliplr(Z_dep_in)], ...
                        bcol, 'FaceAlpha', 2*fp_alpha, 'EdgeColor','none', 'HandleVisibility','off');

                    % Deposition boundaries (dashed, more opaque)
                    plot(R_dep_in,  Z_dep_in,  '--', 'Color', [bcol, 0.55], 'LineWidth', 1.0, 'HandleVisibility','off');
                    plot(R_dep_out, Z_dep_out, '--', 'Color', [bcol, 0.55], 'LineWidth', 1.0, 'HandleVisibility','off');
                end
            end

            % --- 1D beam line (solid, drawn last = on top) ---
            plot(R_env, Z_env, '-', 'Color', bcol, 'LineWidth', 2.8, 'HandleVisibility','off');
        else
            plot([R_active, R_tgt], [Z_active, Z_tgt], '-', 'Color', bcol, 'LineWidth', 2.0, ...
                'HandleVisibility','off');
        end

        % --- Active launcher mirror graphic ---
        ma = atan2(Z_tgt-Z_active, R_tgt-R_active);  mw = 0.28;
        plot([R_active-mw*sin(ma), R_active+mw*sin(ma)], ...
             [Z_active+mw*cos(ma), Z_active-mw*cos(ma)], ...
             '-', 'Color', [.30 .30 .35], 'LineWidth', 5, 'HandleVisibility','off');
        % Label: EL to the right at mirror level, UL above mirror horizontally centered
        if is_eccd  % active = UL
            text(R_active+0.10, Z_active+0.25, active_label, 'FontSize', 9, ...
                'FontWeight','bold', 'Color', bcol, 'HorizontalAlignment','center');
        else        % active = EL
            text(R_active+0.20, Z_active, active_label, 'FontSize', 9, ...
                'FontWeight','bold', 'Color', bcol, 'VerticalAlignment','middle');
        end

        % --- Inactive launcher mirror graphic (gray, at boresight angle) ---
        if is_eccd
            R_inact = R_EL; Z_inact = Z_EL; inact_lbl = 'EL';
        else
            R_inact = R_UL; Z_inact = Z_UL; inact_lbl = 'UL';
        end
        gray_col = [0.65 0.65 0.68];
        ma_inact = atan2(0 - Z_inact, R_0 - R_inact);
        plot([R_inact-mw*sin(ma_inact), R_inact+mw*sin(ma_inact)], ...
             [Z_inact+mw*cos(ma_inact), Z_inact-mw*cos(ma_inact)], ...
             '-', 'Color', gray_col, 'LineWidth', 4, 'HandleVisibility','off');
        if is_eccd  % inactive = EL
            text(R_inact+0.20, Z_inact, sprintf('%s:OFF', inact_lbl), ...
                'FontSize', 9, 'FontWeight','bold', 'Color', gray_col, ...
                'VerticalAlignment','middle');
        else        % inactive = UL
            text(R_inact+0.10, Z_inact+0.25, sprintf('%s:OFF', inact_lbl), ...
                'FontSize', 9, 'FontWeight','bold', 'Color', gray_col, ...
                'HorizontalAlignment','center');
        end
    end

    % --- Title + subtitle (actuator set) ---
    if is_reverse
        ttl_str = 'ECCD reverse — NTM suppression';
    elseif is_ecrh
        ttl_str = 'ECRH — Plasma heating';
    else
        ttl_str = 'ECCD forward — NTM suppression';
    end
    if multi_mode && ~isempty(act_mode_lbl)
        sub_tag = act_mode_lbl;
    else
        sub_tag = mode_tag;
    end
    act_sub = sprintf('[%s]   B_0 = %.2f T,   \\phi = %.1f°,   f = %.0f GHz', ...
        sub_tag, act_B0, act_phi, act_f);
    title(ttl_str, 'FontSize', 12, 'FontWeight', 'normal', 'Interpreter', 'tex');
    subtitle(act_sub, 'FontSize', 9, 'Color', [.40 .40 .45], 'Interpreter', 'tex');

    xlabel('R [m]', 'FontSize', 14);
    ylabel('Z [m]', 'FontSize', 14);
    if isfield(cfg,'ShowGrid') && cfg.ShowGrid
        grid on; set(gca, 'GridAlpha', 0.08);
    end
    set(gca, 'FontSize', 11);
    xlim([4.0, 9.5]); ylim([-3.5, 4.5]);

    if isfield(cfg,'ShowLegend') && cfg.ShowLegend
        [lg, icons] = legend('Location','southeast','FontSize', 8);
        lg.ItemTokenSize = [30, 30];
        % Enlarge pentagram in legend using documented icons output
        pent_icons = findobj(icons, 'Marker', 'pentagram');
        for lm = 1:numel(pent_icons)
            pent_icons(lm).MarkerSize = 14;
        end
    end
    hold off; drawnow;

    % Export
    if isfield(cfg,'ExportFigures') && cfg.ExportFigures
        if multi_mode
            export_fig(figs(ip), sprintf('ATLAS_%s_%s_rank%d_%s', op_label, act_mode_lbl, act_rank, actuator_tag(act_B0, act_phi, act_f)), cfg);
        else
            export_fig(figs(ip), sprintf('ATLAS_%s_%s_rank%d_%s', op_label, mode_tag, act_rank, actuator_tag(act_B0, act_phi, act_f)), cfg);
        end
    end
end

results.figures = figs;

% --- Summary matrix ---
if S.nv >= 2
    fig_matrix = create_matrix(S, n_plots, ri_list, is_reverse, is_point, is_ecrh, multi_mode, ...
        q_target, r_user, theta_user, operation_mode, active_label, ...
        R_active, Z_active, active_bore, theta_L_sign, Q_c, Q_95, a, use_q_filter, cfg);
    results.figures_matrix = fig_matrix;
    if isfield(cfg,'ExportFigures') && cfg.ExportFigures
        export_fig(fig_matrix, sprintf('ATLAS_%s_%s_matrix', op_label, mode_tag), cfg);
    end
end

%% ================================================================
%%  [9] PACK RESULTS
%% ================================================================
results.version = VERSION;
results.operation_mode = operation_mode;
results.targeting_mode = targeting_mode;
results.target_type = target_type;
results.nres_ref = nres_ref;
results.nres_thr = nres_thr;
results.eta = eta;
results.scan = S;
results.scan_full = S_full;
results.cfg = cfg;
results.csv_file = csv_fn;

% .mat export
mat_fn = sprintf('ATLAS_%s_%s_%s_%s.mat', operation_mode, targeting_mode, mode_tag, ts);
save(mat_fn, 'results', '-v7.3');
fprintf('  MAT: %s\n', mat_fn);
fprintf('\n  ATLAS complete.\n\n');

end  % ===== END MAIN FUNCTION =====


%% ================================================================
%%  LANDSCAPE MODE
%% ================================================================

function results = run_landscape(cfg, VERSION)
%RUN_LANDSCAPE  Resonance Landscape Mapper — raw n_res visualization.
%  Produces poloidal cross-section plots for up to 6 user-specified
%  actuator triplets {B0, phi, f}, across one or multiple propagation modes.
%  No ranking, no beam, no CSV.  Both launchers displayed as OFF.

    fprintf('\n  === Resonance Landscape Mapper ===\n');

    % Geometry
    a = cfg.a;  R_0 = cfg.R_0;  k_e = cfg.k_e;
    if isfield(cfg,'beta_p'), beta_p=cfg.beta_p; else, beta_p=0.65; end
    if isfield(cfg,'l_i'), l_i=cfg.l_i; else, l_i=0.85; end
    delta_sh = @(r) (a^2./(2*R_0)).*(beta_p+l_i/2).*(1-(r./a).^2);
    kappa_fn = @(r) 1+(k_e-1)*(r./a).^2;
    Rfun = @(r,th) R_0 + r.*cos(th) + delta_sh(r);
    Zfun = @(r,th) kappa_fn(r).*r.*sin(th);

    % Launchers (displayed as OFF)
    if isfield(cfg,'EL_R'), R_EL=cfg.EL_R; else, R_EL=R_0+a+0.2; end
    if isfield(cfg,'EL_Z'), Z_EL=cfg.EL_Z; else, Z_EL=0; end
    if isfield(cfg,'UL_R'), R_UL=cfg.UL_R; else, R_UL=7.1; end
    if isfield(cfg,'UL_Z'), Z_UL=cfg.UL_Z; else, Z_UL=4.1; end

    % Display
    if isfield(cfg,'color_percentile'), color_pct=cfg.color_percentile; else, color_pct=95; end
    if isfield(cfg,'r_max_display'), r_max_disp=cfg.r_max_display; else, r_max_disp=1.90; end
    if isfield(cfg,'MarkerSize'), ms_scatter=cfg.MarkerSize*0.4; else, ms_scatter=16; end
    if isfield(cfg,'color_quantity'), color_qty=lower(cfg.color_quantity); else, color_qty='nres'; end
    use_fres = strcmp(color_qty, 'fres');

    % Resolve mode list
    if isfield(cfg,'csvPaths') && iscell(cfg.csvPaths)
        csv_list = cfg.csvPaths;
        if isfield(cfg,'mode_labels'), mode_labels = cfg.mode_labels;
        else, mode_labels = arrayfun(@(i) sprintf('Mode%d',i), 1:numel(csv_list), 'Uni', false); end
    else
        csv_list = {cfg.csvPath};
        [~, fn_csv] = fileparts(cfg.csvPath);
        tok = regexp(fn_csv, '(O[12]|X[12])$', 'match');
        if ~isempty(tok), mode_labels = tok(1); else, mode_labels = {'landscape'}; end
    end
    n_modes = numel(csv_list);

    % Get actuator triplets from cfg
    if ~isfield(cfg, 'landscape_actuators') || isempty(cfg.landscape_actuators)
        error('ATLAS:NoActuators', 'cfg.landscape_actuators required for landscape mode.');
    end
    acts = cfg.landscape_actuators;
    n_acts = min(size(acts, 1), 6);

    fprintf('  Modes: %d,  Actuators: %d  ->  %d plots\n', n_modes, n_acts, n_modes*n_acts);

    % Geometry overlays (shared across all plots)
    Q_c = cfg.Q_c;  Q_95 = cfg.Q_95;
    thg = linspace(0, 2*pi, 361);
    R_sep = Rfun(r_max_disp + 0.05, thg);  Z_sep = Zfun(r_max_disp + 0.05, thg);
    qv = [4/3, 3/2, 2]; ql = {'4/3','3/2','2'};
    qc = {[.4 .4 .6],[.7 .15 .15],[.3 .6 .3]};
    qr = NaN(size(qv)); qRc = cell(size(qv)); qZc = cell(size(qv));
    for iq = 1:numel(qv)
        fr = (qv(iq)-Q_c)/(Q_95-Q_c);
        if fr >= 0 && fr <= 1
            qr(iq) = a*sqrt(fr);
            qRc{iq} = Rfun(qr(iq),thg); qZc{iq} = Zfun(qr(iq),thg);
        end
    end

    UL_bore_angle = atan2(0 - Z_UL, R_0 - R_UL);
    EL_bore_angle = atan2(0 - Z_EL, R_0 - R_EL);

    ts = datestr(now, 'yyyymmdd_HHMMSS');
    figs = gobjects(n_modes * n_acts, 1);
    fig_count = 0;

    for im = 1:n_modes
        csvPath = csv_list{im};
        mtag = mode_labels{im};
        if ~isfile(csvPath)
            fprintf('  [%s] File not found: %s — skipping.\n', mtag, csvPath);
            continue;
        end
        fprintf('\n  --- %s ---\n', mtag);
        fprintf('  Loading: %s\n', csvPath);
        T = readtable(csvPath);
        fprintf('  %d rows\n', height(T));

        r_all = T.r_m_; theta_all = T.theta_deg_; R_all = T.R_m_;
        B0_all = T.B0_T_; phi_all = T.phi_deg_; f_all = T.f_GHz_;
        nres_all = T.n_res;
        if ismember('f_res', string(T.Properties.VariableNames))
            fres_all_ls = T.f_res;
        else
            fres_all_ls = zeros(height(T),1);
        end
        clear T;

        gv = isfinite(R_all) & isfinite(nres_all) & isfinite(r_all) & (r_all <= r_max_disp);
        r_all=r_all(gv); theta_all=theta_all(gv); R_all=R_all(gv);
        B0_all=B0_all(gv); phi_all=phi_all(gv); f_all=f_all(gv);
        nres_all=nres_all(gv); fres_all_ls=fres_all_ls(gv);
        R_geom = Rfun(r_all, deg2rad(theta_all));
        Z_geom = Zfun(r_all, deg2rad(theta_all));

        act_keys = round(B0_all,6)*1e12 + round(phi_all,6)*1e6 + round(f_all,4);

        for ia = 1:n_acts
            aB0 = acts(ia, 1); aphi = acts(ia, 2); af = acts(ia, 3);
            fprintf('  [%s %d/%d] B0=%.2f T, phi=%.1f°, f=%.0f GHz\n', mtag, ia, n_acts, aB0, aphi, af);

            akey = round(aB0,6)*1e12 + round(aphi,6)*1e6 + round(af,4);
            mask = (act_keys == akey);
            if sum(mask) == 0
                fprintf('    No data for this actuator set in %s — skipping.\n', mtag);
                continue;
            end
            nres_a = nres_all(mask); fres_a = fres_all_ls(mask);
            Rg = R_geom(mask); Zd = Z_geom(mask);

            % Resolve color quantity
            if use_fres
                cdata_a = fres_a;
                cbar_label_ls = 'f_{res} = n_{res}/n_r';
            else
                cdata_a = nres_a;
                cbar_label_ls = 'n_{res} [m^{-3}]';
            end

            cdata_pos = cdata_a(cdata_a > 0);
            if use_fres
                cax_top_ls = min(1.0, prctile(cdata_pos, color_pct));
                if isempty(cdata_pos) || cax_top_ls <= 0, cax_top_ls = 1.0; end
            else
                if isempty(cdata_pos), cax_top_ls = 1;
                else, cax_top_ls = prctile(cdata_pos, color_pct); if cax_top_ls<=0, cax_top_ls=max(cdata_pos); end; end
            end

            fig_count = fig_count + 1;
            figs(fig_count) = figure('Color','w', 'Position', [80+40*fig_count, 60+40*fig_count, 1100, 850]);
            hold on; axis equal;

            % Colormap selection
            %cmap = parula(256);
            cmap = parula_white(256, 0.15);

            if isfield(cfg,'InvertColormap') && cfg.InvertColormap
                cmap = flipud(cmap);
            end

            % Plasma interior fill
            if isfield(cfg,'FillPlasma') && cfg.FillPlasma
                patch(R_sep, Z_sep, cmap(1,:), 'EdgeColor', 'none', 'HandleVisibility', 'off');
            end

            plot(R_sep, Z_sep, 'k-', 'LineWidth', 1.8, 'HandleVisibility','off');

            cdata_clipped = min(max(cdata_a, 0), cax_top_ls);
            [~, siw] = sortrows([cdata_clipped(:), Rg(:), Zd(:)], [1 2 3]);
            scatter(Rg(siw), Zd(siw), ms_scatter, cdata_clipped(siw), 'filled', ...
                'MarkerFaceAlpha', 0.85, 'MarkerEdgeColor','none', 'HandleVisibility','off');

            % Q-surfaces (drawn on top of scatter)
            for iq = 1:numel(qv)
                if isfinite(qr(iq))
                    plot(qRc{iq}, qZc{iq}, '--', 'Color', qc{iq}, 'LineWidth', 2.2, ...
                        'DisplayName', sprintf('q = %s (r = %.2f m)', ql{iq}, qr(iq)));
                end
            end

            colormap(cmap);
            cb = colorbar; cb.Label.String = cbar_label_ls; cb.Label.FontSize = 12;
            caxis([0, cax_top_ls]);

            plot(R_0, 0, '+', 'MarkerSize', 12, 'LineWidth', 1.5, 'Color', [.4 .4 .4], 'HandleVisibility','off');

            gray_col = [0.65 0.65 0.68]; mw = 0.28;
            ma_ul = UL_bore_angle;
            plot([R_UL-mw*sin(ma_ul), R_UL+mw*sin(ma_ul)], ...
                 [Z_UL+mw*cos(ma_ul), Z_UL-mw*cos(ma_ul)], ...
                 '-', 'Color', gray_col, 'LineWidth', 4, 'HandleVisibility','off');
            text(R_UL+0.10, Z_UL+0.25, 'UL:OFF', 'FontSize', 9, 'FontWeight','bold', 'Color', gray_col, ...
                'HorizontalAlignment','center');
            ma_el = EL_bore_angle;
            plot([R_EL-mw*sin(ma_el), R_EL+mw*sin(ma_el)], ...
                 [Z_EL+mw*cos(ma_el), Z_EL-mw*cos(ma_el)], ...
                 '-', 'Color', gray_col, 'LineWidth', 4, 'HandleVisibility','off');
            text(R_EL+0.20, Z_EL, 'EL:OFF', 'FontSize', 9, 'FontWeight','bold', 'Color', gray_col, ...
                'VerticalAlignment','middle');

            title('Resonance Landscape', 'FontSize', 12, 'FontWeight', 'normal');
            act_sub = sprintf('[%s]   B_0 = %.2f T,   \\phi = %.1f°,   f = %.0f GHz', mtag, aB0, aphi, af);
            subtitle(act_sub, 'FontSize', 9, 'Color', [.40 .40 .45], 'Interpreter', 'tex');

            xlabel('R [m]', 'FontSize', 14);
            ylabel('Z [m]', 'FontSize', 14);
            set(gca, 'FontSize', 11);
            xlim([4.0, 9.5]); ylim([-3.5, 4.5]);
            lg = legend('Location','southeast','FontSize', 8);
            lg.ItemTokenSize = [30, 30];
            hold off; drawnow;

            if isfield(cfg,'ExportFigures') && cfg.ExportFigures
                export_fig(figs(fig_count), sprintf('ATLAS_landscape_%s_%s', mtag, actuator_tag(aB0, aphi, af)), cfg);
            end
        end
    end

    figs = figs(1:fig_count);
    results.version = VERSION;
    results.operation_mode = 'landscape';
    results.figures = figs;
    results.cfg = cfg;
    results.actuators = acts(1:n_acts, :);
    mat_fn = sprintf('ATLAS_landscape_%s.mat', ts);
    save(mat_fn, 'results', '-v7.3');
    fprintf('  MAT: %s\n', mat_fn);
    fprintf('\n  Landscape complete. %d figure(s) produced.\n\n', fig_count);
end


%% ================================================================
%%  SCAN FUNCTIONS
%% ================================================================

function S = scan_forward(r_all, theta_all, R_geom, Z_geom, nres_all, ImN2_all, ...
    B0_scan, phi_scan, f_scan, B0_masks, phi_masks, f_masks, ...
    r_user, theta_user, is_point, q_target, ...
    nres_thr, nres_ref, eta, DeltaR, sigma_theta, adaptive, kappa, ...
    delta_theta_req, theta_gap_floor, tol_r, min_pts, ...
    R_active, Z_active, active_bore, active_steer, Rfun, Zfun, ...
    a, R_0, phi_co_sign, phi_user_set, ...
    pass_dr, pass_search_r, pass_exclude_deg, beam_w0, beam_s_R, ...
    N_data)

    nB=numel(B0_scan); nP=numel(phi_scan); nF=numel(f_scan);
    n_comb = nB*nP*nF;

    % Pre-allocate
    S = alloc_results(n_comb, false);

    % Window mask (target-dependent, actuator-independent)
    if is_point
        mask_window = (abs(r_all-r_user)<=3*DeltaR) & ...
            (min(abs(theta_all-theta_user), 360-abs(theta_all-theta_user)) <= 3*sigma_theta+30);
    else
        mask_window = (abs(r_all-r_user) <= tol_r);
    end

    nv = 0;
    for iB = 1:nB
        mB = B0_masks(:,iB);
        for iP = 1:nP
            % Co-current filter: skip counter-ECCD phi, UNLESS user explicitly set phi
            if ~phi_user_set && phi_co_sign ~= 0 && (phi_co_sign*phi_scan(iP)) < 0, continue; end
            mBP = mB & phi_masks(:,iP);
            if ~any(mBP), continue; end
            for iF = 1:nF
                mask_full = mBP & f_masks(:,iF);
                mask_act = mask_full & mask_window;
                Nt = sum(mask_act);
                if Nt < min_pts, continue; end

                nres_s = nres_all(mask_act);
                ths = theta_all(mask_act);
                rs = r_all(mask_act);
                Rsg = R_geom(mask_act);
                Zsg = Z_geom(mask_act);

                [nres_pk, imax] = max(nres_s);

                % Adaptive sigma
                if adaptive
                    dist_m = sqrt((Rsg(imax)-R_active)^2+(Zsg(imax)-Z_active)^2);
                    sig_eff = rad2deg(sqrt(deg2rad(sigma_theta)^2+(dist_m*kappa/r_user)^2));
                else, sig_eff = sigma_theta; end

                % Kernel weights
                if is_point
                    [w_w, inwin] = compute_kernel_weights(rs, ths, r_user, theta_user, DeltaR, sig_eff);
                    if sum(inwin) < min_pts, continue; end
                    nres_w=nres_s(inwin); ths_w=ths(inwin); w_w=w_w(inwin);
                else
                    nres_w=nres_s; ths_w=ths; w_w=ones(size(nres_s))/Nt;
                end

                above = nres_w >= nres_thr;
                Na = sum(above);

                % ARETE computation
                if is_point
                    % ARETE-p = (nres_eff/nres_ref) x Psi_rob
                    nres_eff_sum = sum(w_w .* nres_w .* double(above));
                    nres_eff = nres_eff_sum / max(nres_ref, eps);

                    if nres_pk > 0
                        rob_thr = eta * nres_pk;
                        robust_mask = nres_w >= rob_thr;
                        if any(robust_mask)
                            rob_ranges = compute_theta_ranges(sort(ths_w(robust_mask)), theta_gap_floor);
                            dth_eff = compute_coverage_deg(rob_ranges);
                        else, dth_eff = 0; end
                    else, dth_eff = 0; end
                    Psi_rob = min(1, dth_eff/max(delta_theta_req, eps));
                    ARETE_val = nres_eff * Psi_rob;

                    % Coverage (for display)
                    if is_point
                        w_pk = max(w_w);
                        if w_pk > 0, above_sig = above & (w_w > 0.01*w_pk);
                        else, above_sig = above; end
                    else, above_sig = above; end
                    th_ranges = compute_theta_ranges(sort(ths_w(above_sig)), theta_gap_floor);
                else
                    % ARETE-f = (N_ab/N_tot) x (Cov_theta/360)
                    above_sig = above;
                    th_ranges = compute_theta_ranges(sort(ths_w(above_sig)), theta_gap_floor);
                    th_L_ranges = map_thetaP_to_thetaL(th_ranges, r_user, R_active, Z_active, Rfun, Zfun);
                    cov_theta = compute_coverage_deg(th_L_ranges);
                    ARETE_val = (Na/max(Nt,1)) * (cov_theta/360);
                end

                % MOIRA
                moira_val = compute_MOIRA(R_active, Z_active, Rsg(imax), Zsg(imax), ths(imax), ...
                    R_geom(mask_full), Z_geom(mask_full), nres_all(mask_full), ImN2_all(mask_full), ...
                    a, R_0, nres_ref, pass_dr, pass_search_r, pass_exclude_deg, beam_w0, beam_s_R);


                % Store
                nv = nv + 1;
                S.B0(nv)=B0_scan(iB); S.phi(nv)=phi_scan(iP); S.f(nv)=f_scan(iF);
                S.r_opt(nv)=rs(imax); S.theta_opt(nv)=ths(imax);
                S.R_opt(nv)=Rsg(imax); S.Z_opt(nv)=Zsg(imax);
                S.nres_peak(nv)=nres_pk; S.N_above(nv)=Na; S.N_total(nv)=Nt;
                S.ARETE(nv)=max(0,min(1,ARETE_val)); S.MOIRA(nv)=moira_val;
                S.theta_ranges{nv}=th_ranges;
                if is_point, S.coverage_deg(nv)=dth_eff;
                else, S.coverage_deg(nv)=cov_theta; end
                S.sigma_theta_eff(nv)=sig_eff;
            end
        end
    end
    S.nv = nv;
end


function S = scan_reverse(r_all, theta_all, R_geom, Z_geom, nres_all, ImN2_all, ...
    B0_scan, phi_scan, f_scan, B0_masks, phi_masks, f_masks, ...
    nres_thr, nres_ref, theta_gap_floor, tol_r, ...
    R_active, Z_active, active_bore, active_steer, Rfun, Zfun, ...
    a, R_0, Q_c, Q_95, phi_co_sign, sigma_theta, adaptive, kappa, ...
    pass_dr, pass_search_r, pass_exclude_deg, beam_w0, beam_s_R, ...
    N_data)

    nB=numel(B0_scan); nP=numel(phi_scan); nF=numel(f_scan);
    max_results = min(nB*nP*nF*500, 5e6);

    S = alloc_results(max_results, true);
    nv = 0;

    for iB = 1:nB
        mB = B0_masks(:,iB);
        for iP = 1:nP
            if phi_co_sign ~= 0 && (phi_co_sign*phi_scan(iP)) < 0, continue; end
            mBP = mB & phi_masks(:,iP);
            if ~any(mBP), continue; end
            for iF = 1:nF
                mask_act = mBP & f_masks(:,iF);
                if sum(mask_act) == 0, continue; end

                nres_s_act = nres_all(mask_act);
                rs_act = r_all(mask_act);
                ths_act = theta_all(mask_act);
                Rsg_act = R_geom(mask_act);
                Zsg_act = Z_geom(mask_act);
                ImN2_act = ImN2_all(mask_act);

                r_unique = unique(round(rs_act, 4));

                for ir = 1:numel(r_unique)
                    r_s = r_unique(ir);
                    on_surface = abs(rs_act - r_s) < tol_r;
                    if sum(on_surface) < 2, continue; end

                    nres_ss = nres_s_act(on_surface);
                    ths_s = ths_act(on_surface);
                    Rsg_s = Rsg_act(on_surface);
                    Zsg_s = Zsg_act(on_surface);

                    above_s = nres_ss >= nres_thr;
                    if ~any(above_s), continue; end

                    % Vectorized launcher coverage
                    viable_idx = find(above_s);
                    R_viable = Rsg_s(viable_idx);
                    Z_viable = Zsg_s(viable_idx);
                    thL_abs = atan2d(Z_viable-Z_active, R_active-R_viable);
                    dev = thL_abs - active_bore;
                    reachable = (dev >= active_steer(1)) & (dev <= active_steer(2));
                    if ~any(reachable), continue; end

                    viable_idx = viable_idx(reachable);
                    ths_viable = ths_s(viable_idx);
                    nres_viable = nres_ss(viable_idx);
                    thL_reached = thL_abs(reachable);


                    % Theta ranges & coverage (reachable — for display and MOIRA target selection)
                    th_rng = compute_theta_ranges(sort(ths_viable), theta_gap_floor);
                    th_L_rng = map_thetaP_to_thetaL(th_rng, r_s, R_active, Z_active, Rfun, Zfun);
                    cov = compute_coverage_deg(th_L_rng);

                    % ARETE-f (aligned with forward scan: ALL above-threshold,
                    %          coverage in launcher-angle space)
                    ths_above_all = ths_s(above_s);
                    th_rng_arete = compute_theta_ranges(sort(ths_above_all), theta_gap_floor);
                    th_L_rng_arete = map_thetaP_to_thetaL(th_rng_arete, r_s, R_active, Z_active, Rfun, Zfun);
                    cov_theta_L = compute_coverage_deg(th_L_rng_arete);
                    N_total_surf = sum(on_surface);
                    N_above_surf = sum(above_s);
                    ARETE_val = (N_above_surf/max(N_total_surf,1)) * (cov_theta_L/360);

                    % MOIRA for all viable, pick best by nres_norm * MOIRA
                    nres_norm = nres_viable / nres_ref;
                    moira_viable = ones(size(nres_viable));
                    for iv = 1:numel(nres_viable)
                        R_v = Rfun(r_s, deg2rad(ths_viable(iv)));
                        Z_v = Zfun(r_s, deg2rad(ths_viable(iv)));
                        moira_viable(iv) = compute_MOIRA(R_active, Z_active, R_v, Z_v, ths_viable(iv), ...
                            Rsg_act, Zsg_act, nres_s_act, ImN2_act, ...
                            a, R_0, nres_ref, pass_dr, pass_search_r, pass_exclude_deg, beam_w0, beam_s_R);
                    end

                    [~, best] = max(nres_norm .* moira_viable);
                    th_best = ths_viable(best);
                    R_best = Rfun(r_s, deg2rad(th_best));
                    Z_best = Zfun(r_s, deg2rad(th_best));

                    if adaptive
                        dist_m = sqrt((R_best-R_active)^2+(Z_best-Z_active)^2);
                        sig_eff = rad2deg(sqrt(deg2rad(sigma_theta)^2+(dist_m*kappa/r_s)^2));
                    else, sig_eff = sigma_theta; end

                    nv = nv + 1;
                    S.B0(nv)=B0_scan(iB); S.phi(nv)=phi_scan(iP); S.f(nv)=f_scan(iF);
                    S.r_opt(nv)=r_s; S.theta_opt(nv)=th_best;
                    S.R_opt(nv)=R_best; S.Z_opt(nv)=Z_best;
                    S.nres_peak(nv)=nres_viable(best); S.N_above(nv)=N_above_surf; S.N_total(nv)=N_total_surf;
                    S.ARETE(nv)=max(0,min(1,ARETE_val)); S.MOIRA(nv)=moira_viable(best);
                    S.theta_ranges{nv}=th_rng; S.coverage_deg(nv)=cov;
                    S.sigma_theta_eff(nv)=sig_eff;
                    S.q_at_target(nv) = Q_c + (Q_95-Q_c)*(r_s/a)^2;

                    if nv >= max_results, break; end
                end
                if nv >= max_results, break; end
            end
            if nv >= max_results, break; end
        end
        if nv >= max_results, break; end
    end
    S.nv = nv;
end


function S = scan_ecrh(r_all, theta_all, R_geom, Z_geom, nres_all, ImN2_all, ...
    B0_scan, phi_scan, f_scan, B0_masks, phi_masks, f_masks, ...
    r_user, nres_thr, nres_ref, theta_gap_floor, tol_r, min_pts, ...
    R_active, Z_active, Rfun, Zfun, a, R_0, ...
    pass_dr, pass_search_r, pass_exclude_deg, beam_w0, beam_s_R, ...
    N_data)

    nB=numel(B0_scan); nP=numel(phi_scan); nF=numel(f_scan);
    n_comb = nB*nP*nF;

    S = alloc_results(n_comb, false);
    mask_window = (r_all <= r_user + tol_r);

    % Total enclosed area for f_spatial normalization (elliptical cross-section)
    A_disc = pi * r_user^2;  % approximate disc area in (r,theta) space

    nv = 0;
    for iB = 1:nB
        mB = B0_masks(:,iB);
        for iP = 1:nP
            mBP = mB & phi_masks(:,iP);
            if ~any(mBP), continue; end
            for iF = 1:nF
                mask_full = mBP & f_masks(:,iF);
                mask_act = mask_full & mask_window;
                Nt = sum(mask_act);
                if Nt < min_pts, continue; end

                nres_s = nres_all(mask_act);
                ths = theta_all(mask_act);
                rs = r_all(mask_act);
                Rsg = R_geom(mask_act);
                Zsg = Z_geom(mask_act);

                [nres_pk, imax] = max(nres_s);
                above = nres_s >= nres_thr;
                Na = sum(above);

                % HeatEff
                if Na > 0 && nres_ref > 0
                    heat_eff = sum(nres_s(above)) / (Nt * nres_ref);
                else, heat_eff = 0; end

                % f_spatial: convex hull area of above-threshold points / disc area
                % Uses raw (R,Z) coordinates — no binning
                if Na >= 3
                    R_ab = Rsg(above); Z_ab = Zsg(above);
                    try
                        [~, A_hull] = convhull(R_ab, Z_ab);
                        % Normalize by total disc area in (R,Z) space
                        A_disc_RZ = pi * r_user * (r_user * 1.7);  % approx for elongation
                        f_spatial = min(1, A_hull / A_disc_RZ);
                    catch
                        % Collinear or degenerate — use radial spread
                        f_spatial = range(rs(above)) / r_user;
                    end
                elseif Na > 0
                    f_spatial = range(rs(above)) / max(r_user, eps);
                else
                    f_spatial = 0;
                end

                ARETE_val = heat_eff * f_spatial;

                % MOIRA
                moira_val = compute_MOIRA(R_active, Z_active, Rsg(imax), Zsg(imax), ths(imax), ...
                    R_geom(mask_full), Z_geom(mask_full), nres_all(mask_full), ImN2_all(mask_full), ...
                    a, R_0, nres_ref, pass_dr, pass_search_r, pass_exclude_deg, beam_w0, beam_s_R);

                thL_abs = atan2d(Zsg(imax)-Z_active, R_active-Rsg(imax));

                nv = nv + 1;
                S.B0(nv)=B0_scan(iB); S.phi(nv)=phi_scan(iP); S.f(nv)=f_scan(iF);
                S.r_opt(nv)=rs(imax); S.theta_opt(nv)=ths(imax);
                S.R_opt(nv)=Rsg(imax); S.Z_opt(nv)=Zsg(imax);
                S.nres_peak(nv)=nres_pk; S.N_above(nv)=Na; S.N_total(nv)=Nt;
                S.ARETE(nv)=max(0,min(1,ARETE_val)); S.MOIRA(nv)=moira_val;
                S.theta_ranges{nv}=zeros(0,2); S.coverage_deg(nv)=0;
            end
        end
    end
    S.nv = nv;
end


%% ================================================================
%%  ALLOCATION HELPER
%% ================================================================

function S = alloc_results(n, has_q)
    S.B0=NaN(n,1); S.phi=NaN(n,1); S.f=NaN(n,1);
    S.r_opt=NaN(n,1); S.theta_opt=NaN(n,1);
    S.R_opt=NaN(n,1); S.Z_opt=NaN(n,1);
    S.nres_peak=NaN(n,1); S.N_above=zeros(n,1); S.N_total=zeros(n,1);
    S.ARETE=NaN(n,1); S.MOIRA=ones(n,1);
    S.theta_ranges=cell(n,1);
    S.coverage_deg=NaN(n,1); S.sigma_theta_eff=NaN(n,1);
    if has_q, S.q_at_target=NaN(n,1); end
    S.nv = 0;
end


%% ================================================================
%%  FALLBACK MECHANISM
%% ================================================================

function results = run_fallback(cfg, operation_mode, targeting_mode, is_point, ...
    r_user, theta_user, q_target, use_q_filter, q_target_surface, r_q_target, ...
    r_all, theta_all, R_geom, Z_geom, nres_all, ImN2_all, ...
    B0_all, phi_all, f_all, B0u, phiu, fu, ...
    nres_thr, nres_ref, eta, DeltaR, sigma_theta, adaptive, kappa, ...
    delta_theta_req, theta_gap_floor, tol_r, min_pts, ...
    R_active, Z_active, active_bore, active_steer, Rfun, Zfun, ...
    a, R_0, Q_c, Q_95, phi_co_sign, ...
    pass_dr, pass_search_r, pass_exclude_deg, beam_w0, beam_s_R, ...
    B0_nom, f_ref, N_data, VERSION)
%RUN_FALLBACK  Propose alternative configurations when primary scan fails.
%
%  Forward: proposes top 5 actuator sets from unconstrained scan.
%  Reverse: (A) best surfaces with user's {B0,f}, (B) nearest {B0,f} for user's q.

    results.version = VERSION;
    results.fallback = true;
    results.operation_mode = operation_mode;
    results.targeting_mode = targeting_mode;

    is_reverse = strcmp(targeting_mode, 'reverse');

    if is_reverse && use_q_filter
        fprintf('\n  --- Fallback (A): Your {B0,f} -> best reachable surfaces (any q) ---\n');
        % Re-run reverse scan WITHOUT q-filter
        % (already done in the main flow — the q_keep filter removed them)
        % For simplicity, report that the user should try without q-filter
        fprintf('  Suggestion: re-run without q_target to see all viable surfaces.\n');

        fprintf('\n  --- Fallback (B): q=%.2f surface -> nearest {B0,f} combos ---\n', q_target_surface);
        % Find all {B0, f} combos that produce any above-threshold data near r_q_target
        B0_range = max(B0u) - min(B0u); if B0_range==0, B0_range=1; end
        f_range = max(fu) - min(fu); if f_range==0, f_range=1; end

        proposal = [];
        for iB = 1:numel(B0u)
            for iF = 1:numel(fu)
                mask = abs(B0_all-B0u(iB))<1e-6 & abs(f_all-fu(iF))<1e-6 & ...
                       abs(r_all-r_q_target)<2*DeltaR;
                if sum(mask) < 2, continue; end
                nres_s = nres_all(mask);
                Na = sum(nres_s >= nres_thr);
                if Na == 0, continue; end
                dist = abs(B0u(iB)-B0_nom)/B0_range + abs(fu(iF)-f_ref)/f_range;
                proposal = [proposal; B0u(iB), fu(iF), Na, sum(mask), dist]; %#ok<AGROW>
            end
        end

        if ~isempty(proposal)
            [~, si] = sort(proposal(:,5), 'ascend');
            n_show = min(5, size(proposal,1));
            fprintf('  %8s  %10s  %8s  %8s  %8s\n', 'B0 [T]', 'f [GHz]', 'N_above', 'N_total', 'dist');
            for k = 1:n_show
                j = si(k);
                fprintf('  %8.2f  %10.1f  %8d  %8d  %8.3f\n', ...
                    proposal(j,1), proposal(j,2), proposal(j,3), proposal(j,4), proposal(j,5));
            end
            results.fallback_B = proposal(si(1:n_show), :);
        else
            fprintf('  No {B0,f} combos produce above-threshold data near q=%.2f.\n', q_target_surface);
            results.fallback_B = [];
        end

    else
        % Forward fallback: unconstrained top 5
        fprintf('\n  --- Fallback: Top 5 actuator sets (unconstrained) ---\n');
        % Quick scan: for each unique {B0, phi, f}, count N_above near target
        keys = round([B0_all, phi_all, f_all], 6);
        [uk, ~, ic] = unique(keys, 'rows');
        n_combos = size(uk, 1);
        combo_Na = zeros(n_combos, 1);
        combo_Mpk = zeros(n_combos, 1);

        for j = 1:n_combos
            mask = (ic == j);
            if is_point
                mask = mask & (abs(r_all-r_user)<=3*DeltaR) & ...
                    (min(abs(theta_all-theta_user),360-abs(theta_all-theta_user))<=30);
            else
                mask = mask & (abs(r_all-r_user)<=tol_r);
            end
            if sum(mask) < 2, continue; end
            nres_s = nres_all(mask);
            combo_Na(j) = sum(nres_s >= nres_thr);
            combo_Mpk(j) = max(nres_s);
        end

        has_data = combo_Na > 0;
        if any(has_data)
            idx = find(has_data);
            [~, si] = sort(combo_Na(idx), 'descend');
            n_show = min(5, numel(si));
            fprintf('  %8s  %10s  %10s  %8s  %10s\n', 'B0 [T]', 'phi [deg]', 'f [GHz]', 'N_above', 'nres_peak');
            for k = 1:n_show
                j = idx(si(k));
                fprintf('  %8.2f  %10.1f  %10.1f  %8d  %10.2e\n', ...
                    uk(j,1), uk(j,2), uk(j,3), combo_Na(j), combo_Mpk(j));
            end
            results.fallback_forward = [uk(idx(si(1:n_show)),:), combo_Na(idx(si(1:n_show))), combo_Mpk(idx(si(1:n_show)))];
        else
            fprintf('  No actuator set produces above-threshold n_res at this target.\n');
            fprintf('  Consider lowering eta (currently %.2f) or checking the dataset.\n', eta);
            results.fallback_forward = [];
        end
    end

    results.cfg = cfg;
end


%% ================================================================
%%  HELPER FUNCTIONS
%% ================================================================

function [sv, is_fixed, fv] = resolve_param(uv, au, name, unit)
    if isstruct(uv) && isfield(uv, 'range')
        rng = uv.range;
        sel = au(au >= rng(1) & au <= rng(2));
        if isempty(sel)
            [~,i1]=min(abs(au-rng(1))); [~,i2]=min(abs(au-rng(2)));
            sel = unique(au([i1,i2]));
        end
        sv = sel(:)'; is_fixed = (numel(sv)==1); fv = sv(1);
        fprintf('    %s range [%.2f,%.2f] %s: %d values\n', name, rng(1), rng(2), unit, numel(sv));
    elseif isnumeric(uv) && numel(uv) > 1
        snapped = NaN(size(uv));
        for k=1:numel(uv), [~,idx]=min(abs(au-uv(k))); snapped(k)=au(idx); end
        sv = unique(snapped(:)','stable'); is_fixed = (numel(sv)==1); fv = sv(1);
        fprintf('    %s: %s %s (%d values)\n', name, mat2str(sv,3), unit, numel(sv));
    elseif isnumeric(uv) && isscalar(uv) && isfinite(uv)
        [~,idx]=min(abs(au-uv)); fv=au(idx); sv=fv; is_fixed=true;
        fprintf('    %s fixed: %.3f %s\n', name, fv, unit);
    else
        fv=NaN; sv=au; is_fixed=false;
        fprintf('    %s: all (%d values)\n', name, numel(au));
    end
end


function rL = map_thetaP_to_thetaL(ranges, r_user, Rm, Zm, Rfun, Zfun)
    if isempty(ranges), rL=zeros(0,2); return; end
    n=size(ranges,1); rL=zeros(n,2);
    for i=1:n
        th_lo=ranges(i,1); th_hi=ranges(i,2);
        if th_hi<th_lo, th_hi=th_hi+360; end
        th_samp=linspace(th_lo,th_hi,21);
        R_s=Rfun(r_user,deg2rad(th_samp)); Z_s=Zfun(r_user,deg2rad(th_samp));
        tL=atan2d(Z_s-Zm, Rm-R_s);
        rL(i,:)=[min(tL),max(tL)];
    end
end


function ranges = compute_theta_ranges(th_sorted, gap_floor)
    if isempty(th_sorted), ranges=zeros(0,2); return; end
    if numel(th_sorted)==1, ranges=[th_sorted(1),th_sorted(1)]; return; end
    th_sorted=th_sorted(:); d=diff(th_sorted); d_pos=d(d>0);
    if numel(d_pos)>=3, gt=max(3*median(d_pos),gap_floor); else, gt=gap_floor; end
    gt=min(gt,45);
    ranges=zeros(numel(th_sorted),2); nr=0; bs=th_sorted(1);
    for i=2:numel(th_sorted)
        if d(i-1)>gt, nr=nr+1; ranges(nr,:)=[bs,th_sorted(i-1)]; bs=th_sorted(i); end
    end
    nr=nr+1; ranges(nr,:)=[bs,th_sorted(end)]; ranges=ranges(1:nr,:);
    if nr>=2
        wg=(360-th_sorted(end))+th_sorted(1);
        if wg<=gt, ranges(1,1)=ranges(end,1); ranges(end,:)=[]; end
    end
end


function [w, inwin] = compute_kernel_weights(r, theta, r_aim, theta_aim, sigma_R, sigma_theta)
    dr=r-r_aim; w_R=exp(-(dr.^2)/(2*sigma_R^2));
    dtheta=abs(theta-theta_aim); dtheta=min(dtheta,360-dtheta);
    w_theta=exp(-(dtheta.^2)/(2*sigma_theta^2));
    w=w_R.*w_theta;
    inwin=(abs(dr)<=3*sigma_R)&(dtheta<=3*sigma_theta);
    w_sum=sum(w); if w_sum>0, w=w/w_sum; end
end


function cov_deg = compute_coverage_deg(theta_ranges)
    if isempty(theta_ranges), cov_deg=0; return; end
    widths=theta_ranges(:,2)-theta_ranges(:,1);
    widths(widths<0)=widths(widths<0)+360;
    cov_deg=sum(widths);
end


function [moira_val, L_path, path_len] = compute_MOIRA( ...
        R_launcher, Z_launcher, R_target, Z_target, ~, ...
        R_data, Z_data, nres_data, ImN2_data, ...
        a, R_0, ~, pass_dr, search_r, exclude_deg, beam_w0, beam_s_R)

    dR=R_target-R_launcher; dZ=Z_target-Z_launcher;
    beam_len=sqrt(dR^2+dZ^2);
    if beam_len<1e-6, moira_val=1; L_path=0; path_len=0; return; end
    uR=dR/beam_len; uZ=dZ/beam_len;

    n_steps=max(2,ceil(beam_len/pass_dr));
    ds=beam_len/n_steps;
    t_mid=(0.5:n_steps)*ds;

    R_samp=R_launcher+uR*t_mid; Z_samp=Z_launcher+uZ*t_mid;
    r_samp=sqrt((R_samp-R_0).^2+Z_samp.^2);
    in_plasma=(r_samp<=a*1.05);

    r_tgt=sqrt((R_target-R_0)^2+Z_target^2);
    exclude_m=max(r_tgt*deg2rad(exclude_deg),0.05);
    dist_tgt=sqrt((R_samp-R_target).^2+(Z_samp-Z_target).^2);
    valid=in_plasma&(dist_tgt>exclude_m);
    path_len=sum(in_plasma)*ds;

    if ~any(valid), moira_val=1; L_path=0; return; end

    has_imn2=any(abs(ImN2_data)>0);
    if has_imn2, alpha_data=abs(nres_data).*abs(ImN2_data);
    else, alpha_data=abs(nres_data); end

    % Spatial grid
    cell_sz=max(search_r,0.05);
    Rmin=min(R_data)-cell_sz; Zmin=min(Z_data)-cell_sz;
    gi=floor((R_data-Rmin)/cell_sz)+1; gj=floor((Z_data-Zmin)/cell_sz)+1;
    ngi=max(gi); ngj=max(gj);
    grid_key=(gi-1)*ngj+gj;
    [sorted_keys,sort_idx]=sort(grid_key);
    [unique_keys,first_pos]=unique(sorted_keys,'first');
    [~,last_pos]=unique(sorted_keys,'last');
    cell_start=zeros(ngi*ngj,1); cell_end=zeros(ngi*ngj,1);
    cell_start(unique_keys)=first_pos; cell_end(unique_keys)=last_pos;

    w_s=beam_w0*sqrt(1+(t_mid/max(beam_s_R,eps)).^2);
    alpha_path=zeros(size(t_mid));
    valid_idx=find(valid);

    for ki=1:numel(valid_idx)
        k=valid_idx(ki);
        r_lookup=max(w_s(k),search_r); r_lookup2=r_lookup^2;
        ci=floor((R_samp(k)-Rmin)/cell_sz)+1;
        cj=floor((Z_samp(k)-Zmin)/cell_sz)+1;
        nc=ceil(r_lookup/cell_sz);
        ci_lo=max(1,ci-nc); ci_hi=min(ngi,ci+nc);
        cj_lo=max(1,cj-nc); cj_hi=min(ngj,cj+nc);
        cand=[];
        for cci=ci_lo:ci_hi
            for ccj=cj_lo:cj_hi
                ck=(cci-1)*ngj+ccj;
                if cell_start(ck)>0
                    cand=[cand;sort_idx(cell_start(ck):cell_end(ck))]; %#ok<AGROW>
                end
            end
        end
        if isempty(cand), continue; end
        dR2=(R_data(cand)-R_samp(k)).^2+(Z_data(cand)-Z_samp(k)).^2;
        near=dR2<=r_lookup2;
        if ~any(near), continue; end
        w2k=2*w_s(k)^2;
        weights_k=exp(-dR2(near)/w2k);
        w_sum=sum(weights_k);
        if w_sum>0, alpha_path(k)=sum(alpha_data(cand(near)).*weights_k)/w_sum; end
    end

    L_path=sum(alpha_path(valid))*ds;
    alpha_ref=prctile(alpha_data(alpha_data>0), 99);
    if isempty(alpha_ref) || alpha_ref<=0, alpha_ref=max(abs(alpha_data)); end
    if alpha_ref<=0, alpha_ref=1; end
    L_ref=alpha_ref*max(path_len,a);
    if L_ref>0, moira_val=exp(-L_path/L_ref); else, moira_val=1; end
    moira_val=max(0,min(1,moira_val));
end


function [T, n_skipped] = robust_readtable(csvPath)
    try, T=readtable(csvPath); n_skipped=0; return; catch, end
    fid=fopen(csvPath,'r');
    header=fgetl(fid);
    n_expected=numel(strfind(header,','))+1;
    tmp=[tempname '.csv']; fid_tmp=fopen(tmp,'w');
    fprintf(fid_tmp,'%s\n',header);
    n_skipped=0;
    while ~feof(fid)
        ln=fgetl(fid);
        if ~ischar(ln)||isempty(ln), continue; end
        if numel(strfind(ln,','))+1==n_expected
            fprintf(fid_tmp,'%s\n',ln);
        else, n_skipped=n_skipped+1; end
    end
    fclose(fid); fclose(fid_tmp);
    T=readtable(tmp); delete(tmp);
end


function str = fmt_union(ranges, inline, show_unit)
    if nargin<2, inline=false; end
    if nargin<3, show_unit=true; end
    if isempty(ranges), str='(none)'; return; end
    if show_unit, suffix=' deg'; else, suffix=''; end
    lo=arrayfun(@(x)mod(x+180,360)-180, ranges(:,1));
    hi=arrayfun(@(x)mod(x+180,360)-180, ranges(:,2));
    tol=5e-3; is_single=abs(lo-hi)<=tol;
    parts=cell(size(ranges,1),1);
    for i=1:size(ranges,1)
        if is_single(i), parts{i}=sprintf('{%.1f}',lo(i));
        else, parts{i}=sprintf('[%.1f, %.1f]',lo(i),hi(i)); end
    end
    if inline, str=[strjoin(parts,' U ') suffix];
    else
        str=parts{1};
        for i=2:numel(parts), str=sprintf('%s U %s',str,parts{i}); end
        str=[str suffix];
    end
end


function fig = create_matrix(S, n_plots, ri_list, is_reverse, is_point, is_ecrh, multi_mode, ...
        q_target, r_user, theta_user, operation_mode, active_label, ...
        R_active, Z_active, active_bore, theta_L_sign, Q_c, Q_95, a, use_q_filter, cfg)
%CREATE_MATRIX  Summary ranking table as a figure.

    % Build display row indices: top overall + best-per-mode if not already present
    base_n = min(10, S.nv);
    show_idx = (1:base_n)';
    if multi_mode
        for k = 1:numel(ri_list)
            if ~ismember(ri_list(k), show_idx)
                show_idx(end+1) = ri_list(k); %#ok<AGROW>
            end
        end
    end
    n_show = numel(show_idx);

    % Determine ARETE suffix
    if is_ecrh, arete_lbl = 'ARETE_{h}';
    elseif is_point, arete_lbl = 'ARETE_{p}';
    else, arete_lbl = 'ARETE_{f}'; end

    % Build column structure
    reverse_show_rq = is_reverse && ~use_q_filter;
    col_names = {'#'};
    cw_base = [.04];
    if multi_mode
        col_names{end+1} = 'Mode'; cw_base(end+1) = .07;
    end
    if reverse_show_rq
        col_names = [col_names, {'r [m]','q','B_0 [T]','\phi [°]','f [GHz]','\theta_L [°]'}];
        cw_base = [cw_base, .07, .06, .08, .08, .08, .09];
    elseif is_ecrh
        col_names = [col_names, {'B_0 [T]','\phi [°]','f [GHz]'}];
        cw_base = [cw_base, .13, .13, .13];
    else
        col_names = [col_names, {'B_0 [T]','\phi [°]','f [GHz]','\theta_L [°]'}];
        cw_base = [cw_base, .11, .11, .11, .12];
    end
    col_names = [col_names, {arete_lbl, 'MOIRA', 'Score'}];
    cw_base = [cw_base, .12, .12, .11];
    n_cols = numel(col_names);
    cw_base = cw_base / sum(cw_base);  % normalize to 1

    fig = figure('Color','w','Name','ATLAS — Ranking','Position',[120,120,1200,min(280+n_show*42,850)]);

    % Title: operation label + target info subtitle
    if is_reverse
        ttl_main = sprintf('ECCD reverse — Top %d', n_show);
        % Reference configuration
        if isfield(cfg,'B0_reverse'), B0_ref=cfg.B0_reverse; else, B0_ref=5.0; end
        if isfield(cfg,'f_reverse'), f_ref=cfg.f_reverse; else, f_ref=170; end
        if isstruct(B0_ref), ref_B0_str = sprintf('%.1f:%.1f T', B0_ref.range(1), B0_ref.range(2));
        elseif numel(B0_ref)>1, ref_B0_str = sprintf('[%.1f–%.1f] T', min(B0_ref), max(B0_ref));
        else, ref_B0_str = sprintf('%.2f T', B0_ref); end
        if isstruct(f_ref), ref_f_str = sprintf('%.0f:%.0f GHz', f_ref.range(1), f_ref.range(2));
        elseif numel(f_ref)>1, ref_f_str = sprintf('[%.0f–%.0f] GHz', min(f_ref), max(f_ref));
        else, ref_f_str = sprintf('%.0f GHz', f_ref); end
        ref_str = sprintf('B_0 = %s,  f = %s', ref_B0_str, ref_f_str);
        if isfield(cfg,'q_target') && isfinite(cfg.q_target)
            frac_q = (cfg.q_target - Q_c)/(Q_95 - Q_c);
            if frac_q >= 0 && frac_q <= 1
                r_q = a*sqrt(frac_q);
                ttl_sub = sprintf('%s  |  q = %s surface (r = %.2f m)   [%s]', ref_str, q_rational(cfg.q_target), r_q, active_label);
            else
                ttl_sub = sprintf('%s   [%s]', ref_str, active_label);
            end
        else
            ttl_sub = sprintf('%s   [%s]', ref_str, active_label);
        end
    elseif is_ecrh
        ttl_main = sprintf('ECRH — Top %d', n_show);
        ttl_sub = sprintf('r \\leq %.2f m   [%s]', r_user, active_label);
    elseif is_point
        ttl_main = sprintf('ECCD forward — Top %d', n_show);
        ttl_sub = sprintf('r = %.2f m,  \\theta = %.0f\\circ   [%s]', ...
            r_user, mod(theta_user+180,360)-180, active_label);
    else
        ttl_main = sprintf('ECCD forward — Top %d', n_show);
        ttl_sub = sprintf('q = %s surface  (r = %.2f m)   [%s]', q_rational(q_target), r_user, active_label);
    end
    sgtitle(ttl_main, 'FontSize', 12, 'FontWeight', 'bold');
    annotation('textbox', [0.0, 0.885, 1.0, 0.04], 'String', ttl_sub, ...
        'FontSize', 11, 'Color', [.40 .40 .45], 'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center', 'Interpreter', 'tex', 'FitBoxToText', 'off');

    hdr_col = [0.18 0.28 0.48];
    lm = 0.03; cw = cw_base * (1 - 2*lm);
    col_x = cumsum([lm, cw]);
    row_h = 0.055; hdr_h = 0.065;
    top_m = 0.15;

    ax = axes('Position',[0,0,1,1],'Visible','off'); hold on;

    % Header
    y = 1 - top_m;
    for c = 1:n_cols
        rectangle('Position',[col_x(c), y-hdr_h, cw(c), hdr_h], ...
            'FaceColor', hdr_col, 'EdgeColor', [.1 .1 .1], 'LineWidth', 1.2);
        text(col_x(c)+cw(c)/2, y-hdr_h/2, col_names{c}, ...
            'HorizontalAlignment','center','VerticalAlignment','middle', ...
            'FontSize', 10, 'FontWeight','bold', 'Color','w', 'Interpreter','tex');
    end

    y = y - hdr_h;
    for r = 1:n_show
        si = show_idx(r);  % actual index into S
        if mod(r,2)==0, rc=[.90 .93 .97]; else, rc=[.96 .97 1.0]; end
        is_sel = ismember(si, ri_list);
        ec = iif(is_sel, [.1 .6 .2], [.7 .7 .7]);
        elw = iif(is_sel, 2.0, 0.5);

        for c = 1:n_cols
            rectangle('Position',[col_x(c),y-row_h,cw(c),row_h], ...
                'FaceColor',rc,'EdgeColor',ec,'LineWidth',elw);
        end

        ci = 1;
        rk_str = sprintf('#%d', S.true_rank(si));
        if is_sel, rk_str = [rk_str ' *']; end
        txt(col_x(ci),cw(ci),y,row_h,rk_str,10,'bold',[0 0 0]);
        ci=ci+1;

        % Mode column (multi-mode only)
        if multi_mode
            mode_col = [0.15 0.35 0.60];
            txt(col_x(ci),cw(ci),y,row_h,S.mode_label{si},9,'bold',mode_col);
            ci=ci+1;
        end

        if reverse_show_rq
            txt(col_x(ci),cw(ci),y,row_h,sprintf('%.3f',S.r_opt(si)),9,'normal',[0 0 0]); ci=ci+1;
            q_r = Q_c+(Q_95-Q_c)*(S.r_opt(si)/a)^2;
            txt(col_x(ci),cw(ci),y,row_h,sprintf('%.2f',q_r),9,'normal',[0 0 0]); ci=ci+1;
        end

        txt(col_x(ci),cw(ci),y,row_h,sprintf('%.2f',S.B0(si)),9,'normal',[0 0 0]); ci=ci+1;
        txt(col_x(ci),cw(ci),y,row_h,sprintf('%.1f',S.phi(si)),9,'normal',[0 0 0]); ci=ci+1;
        txt(col_x(ci),cw(ci),y,row_h,sprintf('%.1f',S.f(si)),9,'normal',[0 0 0]); ci=ci+1;

        if ~is_ecrh
            thL = theta_L_sign*(atan2d(S.Z_opt(si)-Z_active,R_active-S.R_opt(si))-active_bore);
            txt(col_x(ci),cw(ci),y,row_h,sprintf('%.1f',thL),9,'normal',[0 0 0]); ci=ci+1;
        end

        % ARETE (color-coded)
        arete_v = S.ARETE(si);
        if arete_v >= 0.01, arete_str = sprintf('%.3f', arete_v);
        elseif arete_v > 0, arete_str = sprintf('%.1e', arete_v);
        else, arete_str = '0'; end
        if arete_v>=0.1, rc_arete=[0 .5 0]; elseif arete_v>=0.005, rc_arete=[.6 .5 0]; else, rc_arete=[.7 .1 .1]; end
        txt(col_x(ci),cw(ci),y,row_h,arete_str,10,'bold',rc_arete); ci=ci+1;

        % MOIRA (color-coded)
        lp = S.MOIRA(si);
        if lp>=0.9, rc_mo=[0 .5 0]; elseif lp>=0.7, rc_mo=[.6 .5 0]; else, rc_mo=[.7 .1 .1]; end
        txt(col_x(ci),cw(ci),y,row_h,sprintf('%.3f',lp),10,'normal',rc_mo); ci=ci+1;

        % Score
        sc_v = S.score(si);
        if sc_v >= 0.01, sc_str = sprintf('%.3f', sc_v);
        elseif sc_v > 0, sc_str = sprintf('%.1e', sc_v);
        else, sc_str = '0'; end
        txt(col_x(ci),cw(ci),y,row_h,sc_str,10,'bold',[0 0 0]);

        y = y - row_h;
    end

    % Footer
    text(0.5, 0.02, '* = plotted', 'HorizontalAlignment','center', ...
        'FontSize', 7, 'FontAngle','italic', 'Color', [.5 .5 .5]);
    xlim([0,1]); ylim([0,1]); hold off;
end


function txt(x, w, y, h, str, fs, fw, col)
    text(x+w/2, y-h/2, str, 'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'FontSize', fs, 'FontWeight', fw, 'Color', col);
end


function export_fig(fig_handle, base_name, cfg)
    dpi = cfg.ExportDPI;
    if isfield(cfg,'ExportFormat') && ismember(cfg.ExportFormat,{'png','both'})
        % Hide titles/subtitles for tight PNG export (kept in .mat figures)
        hidden = hide_titles(fig_handle);
        fn = sprintf('%s.png', base_name);
        try
            exportgraphics(fig_handle, fn, 'Resolution', dpi, 'BackgroundColor', 'white');
        catch
            % Fallback for MATLAB < R2020a: tighten manually
            set(fig_handle, 'PaperPositionMode', 'auto');
            drawnow;
            print(fig_handle, '-dpng', sprintf('-r%d', dpi), fn);
        end
        fprintf('    -> %s\n', fn);
        restore_titles(hidden);
    end
    if isfield(cfg,'ExportFormat') && ismember(cfg.ExportFormat,{'eps','both'})
        fn = sprintf('%s.eps', base_name);
        try
            exportgraphics(fig_handle, fn, 'BackgroundColor', 'white', 'ContentType', 'vector');
        catch
            set(fig_handle, 'PaperPositionMode', 'auto');
            print(fig_handle, '-depsc', sprintf('-r%d', dpi), fn);
        end
        fprintf('    -> %s\n', fn);
    end
end


function hidden = hide_titles(fig)
%HIDE_TITLES  Temporarily hide title/subtitle on poloidal and landscape plots.
%  Hides axes Title and Subtitle objects for tight PNG export.
%  Does NOT touch sgtitle or annotation textboxes — these belong to the
%  matrix figure and must remain visible in exported PNGs.
    hidden = gobjects(0);

    all_ax = findobj(fig, 'Type', 'axes');
    for i = 1:numel(all_ax)
        th = get(all_ax(i), 'Title');
        if ~isempty(th) && ~isempty(get(th, 'String'))
            set(th, 'Visible', 'off');
            hidden(end+1) = th; %#ok<AGROW>
        end
        try %#ok<TRYNC>  — Subtitle property exists R2020b+
            sh = get(all_ax(i), 'Subtitle');
            if ~isempty(sh) && ~isempty(get(sh, 'String'))
                set(sh, 'Visible', 'off');
                hidden(end+1) = sh; %#ok<AGROW>
            end
        end
    end

    drawnow;
end


function restore_titles(hidden)
%RESTORE_TITLES  Re-show objects hidden by hide_titles.
    for i = 1:numel(hidden)
        if isvalid(hidden(i))
            set(hidden(i), 'Visible', 'on');
        end
    end
    drawnow;
end


function result = iif(condition, true_val, false_val)
    if condition, result = true_val; else, result = false_val; end
end


function lbl = compact_op_label(operation_mode, targeting_mode)
%COMPACT_OP_LABEL  Map operation/targeting to compact filename token.
    if strcmp(operation_mode, 'eccd')
        if strcmp(targeting_mode, 'forward'), lbl = 'forward';
        elseif strcmp(targeting_mode, 'reverse'), lbl = 'reverse';
        else, lbl = 'eccd'; end
    elseif strcmp(operation_mode, 'ecrh')
        lbl = 'heating';
    else
        lbl = operation_mode;  % landscape, etc.
    end
end


function tag = actuator_tag(B0, phi, f)
%ACTUATOR_TAG  Compact string for filenames: e.g. '5T_15deg_170GHz'
    if B0 == round(B0)
        b_str = sprintf('%.0fT', B0);
    else
        b_str = sprintf('%.2fT', B0);
    end
    if phi == round(phi)
        p_str = sprintf('%.0fdeg', phi);
    else
        p_str = sprintf('%.1fdeg', phi);
    end
    f_str = sprintf('%.0fGHz', f);
    tag = sprintf('%s_%s_%s', b_str, p_str, f_str);
end



function str = q_rational(q_val)
%Q_RATIONAL  Return rational string for standard q values, decimal otherwise.
    tol = 1e-4;
    if abs(q_val - 4/3) < tol, str = '4/3';
    elseif abs(q_val - 3/2) < tol, str = '3/2';
    elseif abs(q_val - 2) < tol, str = '2';
    elseif abs(q_val - 5/3) < tol, str = '5/3';
    elseif abs(q_val - 5/2) < tol, str = '5/2';
    elseif abs(q_val - 3) < tol, str = '3';
    else, str = sprintf('%.2f', q_val);
    end
end



function cmap = parula_white(N, topFrac)
    if nargin < 1, N = 256; end
    if nargin < 2, topFrac = 0.15; end

    cmap = parula(N);

    k = max(1, round(N * topFrac));
    i0 = N - k + 1;

    c0 = cmap(i0,:);
    white = [1 1 1];

    t = linspace(0,1,k)';
    cmap(i0:end,:) = (1-t).*c0 + t.*white;
end