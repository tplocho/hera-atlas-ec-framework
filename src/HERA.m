% HERA - Hot-plasma Electron-cyclotron Resonance Analysis  v.2.0.0
%
%  EC resonance geometry code with 2D velocity-space integration
%  and mode-resolved polarization-coupled harmonic weighting.
%
%  Resonance curve (Bornatici et al., Nucl. Fusion 23, 1153, 1983, Sec. 2):
%    Dw(v_par,v_perp) = w - k_par*v_par - ell*wc/gamma(v_par,v_perp) = 0
%    => v_perp^2 = c^2[1 - ((w-k_par*v_par)/(ell*wc))^2] - v_par^2
%
%  J2D: mode-resolved coupling-weighted resonant measure [m^-3 s].
%    J2D = n_r * integral f_M * |Psi_ell^(mode)|^2 * 2*pi*c^2/(ell*wc*gamma) dvp
%
%    Coupling factor (Bornatici 1983, Sec. 3; Stix 1992, Ch. 10, Eq. 10-46):
%      Psi_ell = eR*J_{ell-1}(z) + eL*J_{ell+1}(z) + ez*(vp/vr)*J_ell(z)
%      z = k_perp*v_perp/omega_c
%      eR, eL, ez: branch-consistent polarization eigenvector.
%
%    This is a mode-resolved polarization-coupling factor built from the
%    harmonic Bessel structure contracted with the branch-consistent
%    polarization eigenvector.  The polarization is extracted from the
%    null eigenvector of the wave matrix Lambda built from the SAME
%    dielectric model that produced the branch N^2:
%      kinetic branch   -> kinetic Lambda (Harris/Stix tensor, complex N^2)
%      warm-fluid branch -> warm-fluid Lambda (pressure tensor)
%      cold branch       -> cold-plasma Lambda (Stix S,D,P)
%    The wave matrix uses the exact structure from the branch determinant
%    solver (L21=-K12, L31=L13, L32=-K23), not Hermitian conjugates.
%
%  n_res: local resonant electron density [m^-3].
%    n_res = sqrt(2*pi) * |k_par|*v_th * J2D (old formulation)
%    n_res: 2D line-shape integral (new formulation)
%    Ref: Bornatici (1983) Sec. 2.2.
%
%  Integration bounds from ECRH onset (Lopez et al., arXiv:2501.04619, 2025, Eq. 3).
%  Dielectric: nonrelativistic Harris/Stix.  Resonance: relativistic.
%    Ref: Prater, Phys. Plasmas 11, 2349 (2004); Stix (1992) Ch. 10.
%  Shafranov shift: Freidberg (2014) Sec. 6.6; Wesson (2011) Sec. 3.6.
%  Elongation: Miller et al., Phys. Plasmas 5, 973 (1998).
%  Booker quartic: Stix (1992) Eq. 1-56.  Bernstein guard: Laqua (2007).

%% ===== INPUTS ===== %%
a         = 1.95;
R_0       = 6.2;
ell       = 1;
Q_c       = 1;
Q_95      = 4;
k_e       = 1.7;
T_c_keV   = 10;
T_e_keV   = 1;
n_c       = 1e20;
n_e       = 1e19;

% Shafranov shift (Freidberg 2014, Sec. 6.6)
beta_p    = 0.65;
l_i       = 0.85;

snap_frequencies_to_grid = false;
print_live = true;

%% ===== TOGGLES ===== %%
fprintf('\nGeometry & model options (press Enter for defaults):\n');

% Temperature profile
tmp = upper(strtrim(input(['  Temperature profile:\n' ...
    '    (P) Pedestal -- H-mode with edge pedestal\n' ...
    '    (Q) Parabolic -- simple T(r) = T_c - (T_c-T_e)*(r/a)^2\n  [default P]: '],'s')));
if ismember(tmp, {'Q','PARABOLIC'})
    Te_profile = 'parabolic';
    fprintf('      -> Parabolic T_e(r)\n');
else
    Te_profile = 'pedestal';
    fprintf('      -> Pedestal T_e(r)\n');
end
clear tmp;

tmp = lower(strtrim(input('  Enable Shafranov shift? (y/n) [default y]: ','s')));
enable_shafranov = ~ismember(tmp, {'n','no'});

tmp = lower(strtrim(input('  Enable elongation kappa(r)? (y/n) [default y]: ','s')));
enable_elongation = ~ismember(tmp, {'n','no'});

tmp = lower(strtrim(input(['  Refractive index model:\n' ...
    '    (K) Kinetic  -- Harris/Stix + Z-function\n' ...
    '    (F) Fluid    -- warm-fluid pressure-tensor\n' ...
    '    (B) Booker   -- cold-plasma Stix\n  [default K]: '],'s')));
if ismember(tmp, {'b','booker'}),     refr_model = 'booker';
elseif ismember(tmp, {'f','fluid'}),  refr_model = 'warmfluid';
else,                                 refr_model = 'kinetic';
end
clear tmp;

% Wave mode
tmp = upper(strtrim(input('  Wave mode: Only (O) / Only (X) / All (A) [default O]: ','s')));
if ismember(tmp, {'X'}),                 wave_mode_sel = 'X';
elseif ismember(tmp, {'A','B','BOTH'}),  wave_mode_sel = 'A';
else,                                    wave_mode_sel = 'O';
end
clear tmp;

% EC harmonic index
tmp = strtrim(input('  EC harmonic index ell (1, 2, ...) [default 1]: ','s'));
if ~isempty(tmp)
    ell_input = str2double(tmp);
    if isfinite(ell_input) && ell_input >= 1 && ell_input == round(ell_input)
        ell = ell_input;
    end
end
clear tmp;
fprintf('  -> ell = %d\n', ell);

%% ===== COLLISION REG (Booker only) ===== %%
if strcmp(refr_model, 'booker')
    nu_input = input('  nu_frac (0=off, rec=1e-3) [default 1e-3]: ', 's');
    if isempty(strtrim(nu_input)), nu_frac = 1e-3;
    else, nu_frac = str2double(nu_input);
        if ~isfinite(nu_frac)||nu_frac<0, nu_frac=1e-3; end
    end
    tmp = lower(strtrim(input('  Adaptive freq refinement? (y/n) [default y]: ','s')));
    enable_adaptive_freq = ~ismember(tmp, {'n','no'});
    clear tmp;
else
    nu_frac = 0;  enable_adaptive_freq = false;
end

%% ===== PEDESTAL (density) ===== %%
ped_frac = 0.80;  r_ped = 0.94*a;  delta_c = 0.5;  delta_p = 10.0;

%% ===== PEDESTAL (temperature) ===== %%
% Same functional form as n_e pedestal.
% Ref: Doyle et al., Nucl. Fusion 47, S18 (2007);
%      Kinsey et al., Nucl. Fusion 51, 083001 (2011).
T_ped_keV  = 4.0;       % pedestal top [keV]
r_ped_T    = 0.92*a;     % pedestal location
delta_c_T  = 0.5;        % core shape
delta_p_T  = 12.0;       % pedestal steepness

%% ===== CONSTANTS (CODATA) ===== %%
c        = 299792458;
qe       = 1.602176634e-19;
me       = 9.1093837015e-31;
eps0     = 8.8541878128e-12;
kB       = 1.380649e-23;
keV_to_K = 1.160451812e7;

%% ===== SCANS ===== %%
% Radial resolution 0.025 m gives smooth LFS onset boundaries.
% Increase dr for faster scans; decrease for publication quality.
r_range       = 0:0.050:a;
theta_deg_rng = 0:1:360;
B_0_range     = 4.5:0.5:6.0;
phi_deg_rng   = -20:1:20;

BAND_LO_GHZ = 110;  BAND_HI_GHZ = 170;
base_f_GHz_grid = BAND_LO_GHZ:10:BAND_HI_GHZ;
base_w_grid     = base_f_GHz_grid*2*pi*1e9;
base_k_grid     = base_w_grid./c;

theta_rad_rng = deg2rad(theta_deg_rng);
phi_rad_rng   = deg2rad(phi_deg_rng);

N_harm = ell + 1;     % harmonic sum |n| <= N_harm
N_vpar = 64;          % resonance curve quadrature points (1D J2D)
K_mu   = 5;           % max |v|/v_th  (exp(-25) ~ 1e-11)
N_vp2d = 40;          % v_par  grid points (2D line-shape integral)
N_vr2d = 20;          % v_perp grid points (2D line-shape integral)

%% ===== FULL PARAMETRIC SCAN ===== %%
r_it = r_range(r_range > 0);
fprintf('Full parametric scan: %d radial points (r = %.2f to %.2f m)\n', ...
    numel(r_it), r_it(1), r_it(end));

%% ===== SETTINGS SUMMARY ===== %%
fprintf('T_e profile: %s\n', Te_profile);
fprintf('Shafranov: %s', mat2str(enable_shafranov));
if enable_shafranov
    fprintf('  (beta_p=%.2f, l_i=%.2f, Delta_0=%.3f m)\n', ...
        beta_p, l_i, (a^2/(2*R_0))*(beta_p+l_i/2));
else, fprintf('\n'); 
end
fprintf('Elongation: %s\n', mat2str(enable_elongation));
fprintf('Model: %s  |  ell=%d  |  Mode: %s\n', refr_model, ell, wave_mode_sel);
fprintf('Harmonic sum: |n| <= %d  |  1D quadrature: %d pts  |  2D grid: %dx%d\n', ...
    N_harm, N_vpar, N_vp2d, N_vr2d+1);

%% ===== OUTPUT CSV ===== %%
% Build default name from model/mode/ell selections
switch refr_model
    case 'kinetic',   model_tag = 'kinetic';
    case 'warmfluid', model_tag = 'fluid';
    case 'booker',    model_tag = 'booker';
    otherwise,        model_tag = refr_model;
end
if strcmp(wave_mode_sel, 'A'), mode_tag = 'OX';
else,                          mode_tag = wave_mode_sel;
end
default_csv = sprintf('data_%s_%s%d.csv', model_tag, mode_tag, ell);

filename = strtrim(input(sprintf('\nCSV filename [%s]: ', default_csv),'s'));
if isempty(filename), filename = default_csv; end
[~,nm,ex]=fileparts(filename); if isempty(ex), filename=[filename '.csv']; end
if isfile(filename)
    ow=lower(strtrim(input(sprintf('"%s" exists. Overwrite? (y/n) [default n]: ',filename),'s')));
    if ismember(ow,{'y','yes'}), delete(filename);
    else, filename=sprintf('%s_%s%s',nm,datestr(now,'yyyymmdd_HHMMSS'),ex); end
end
fprintf('Output: %s\n', filename);

fid = fopen(filename,'w');
fprintf(fid, 'r (m),theta (deg),R (m),B0 (T),phi (deg),f (GHz),');
fprintf(fid, 'v_par/v_th,v_perp/v_th,J2D,n_res,w/wc,N,N2,ImN2,k_par,mode,f_res\n');
fclose(fid);
fid = fopen(filename,'a');

if print_live
    W = [8 10 8 7 9 8 12 12 12 12 8 8 10 10 10 5 8];
    hdr = {'r','theta','R','B0','phi','f','vp/vth','vr/vth', ...
           'J2D','n_res','w/wc','N','N2','ImN2','kpar','mode','f_res'};
    printRow(hdr, W); fprintf('%s\n', repmat('-',1,sum(W)));
end

root_count = 0;  batch_ctr = 0;  kpar_min = 1e-12;

%% ===== MAIN LOOPS ===== %%
for r = r_it
    q_r  = Q_c + (Q_95-Q_c)*(r^2/a^2);
    n_ped = n_e + ped_frac*(n_c-n_e);
    n_r  = ne_pedestal(r, a, n_c, n_e, n_ped, r_ped, delta_c, delta_p);
    if strcmp(Te_profile, 'pedestal')
        T_r_keV = ne_pedestal(r, a, T_c_keV, T_e_keV, T_ped_keV, r_ped_T, delta_c_T, delta_p_T);
    else
        T_r_keV = T_c_keV - (T_c_keV - T_e_keV) * (r^2/a^2);
    end
    T_r  = T_r_keV * keV_to_K;
    vth  = sqrt(2*kB*T_r/me);

    for theta = theta_rad_rng

        if enable_shafranov
            Delta_sh = (a^2/(2*R_0))*(beta_p+l_i/2)*(1-r^2/a^2);
            delta_r  = -Delta_sh;
            delta_r_der = (beta_p+l_i/2)*r/R_0;
        else, delta_r=0; delta_r_der=0;
        end
        if enable_elongation
            k_r = 1+(k_e-1)*(r^2/a^2);
            r_k_r_der = k_r + 2*(k_e-1)*r^2/a^2;
        else, k_r=1; r_k_r_der=1;
        end

        a_ratio = r/(R_0-delta_r);
        j_r = k_r*cos(theta)*(cos(theta)-delta_r_der) + r_k_r_der*sin(theta)^2;
        R_pos = R_0 + r*cos(theta) - delta_r;

        for B_0 = B_0_range
            denom_Bt = 1+a_ratio*cos(theta);
            if abs(denom_Bt)<1e-14, continue; end
            B_t = B_0/denom_Bt;
            R0s = R_0-delta_r;
            denom_Bp = 1+a_ratio*cos(theta)*j_r;
            if abs(denom_Bp)<1e-14||abs(R0s)<1e-14, continue; end
            B_p = (r*B_t)/(q_r*R0s) ...
                  *(r_k_r_der*sqrt(1+(k_r^2-1)*cos(theta)^2))/denom_Bp;
            B   = hypot(B_t, B_p);

            omega_c = qe*B/me;
            omega_p = sqrt(qe^2*n_r/(me*eps0));
            nu_coll = nu_frac*omega_c;

            if enable_adaptive_freq
                f_c = omega_c/(2*pi*1e9);
                df = [-5,-3,-2,-1.5,-1,-.5,-.2,0,.2,.5,1,1.5,2,3,5];
                fe = f_c+df; fe=fe(fe>=BAND_LO_GHZ & fe<=BAND_HI_GHZ);
                fm = sort([base_f_GHz_grid,fe]); kp=true(size(fm));
                for jj=2:length(fm), if abs(fm(jj)-fm(jj-1))<0.05, kp(jj)=false; end; end
                w_grid = fm(kp)*2*pi*1e9; k_grid = w_grid./c;
            else
                w_grid = base_w_grid; k_grid = base_k_grid;
            end

            for phi = phi_rad_rng
                sphi = sin(phi);
                if abs(sphi)<1e-14, continue; end

                % --- Dispersion ---
                if strcmp(refr_model,'kinetic')
                    [N2_O,N2_X] = n_warm_plasma(w_grid,phi,omega_p,omega_c,vth,c,N_harm,wave_mode_sel);
                    N2_modes={N2_O,N2_X}; mode_labels={'O','X'};
                elseif strcmp(refr_model,'warmfluid')
                    [N2_O,N2_X] = n_warm_fluid(w_grid,phi,omega_p,omega_c,vth,c,wave_mode_sel);
                    N2_modes={N2_O,N2_X}; mode_labels={'O','X'};
                else  % booker
                    [N2_O,N2_X] = n_booker(w_grid,phi,omega_p,omega_c,nu_coll,wave_mode_sel);
                    N2_modes={N2_O,N2_X}; mode_labels={'O','X'};
                end

                switch wave_mode_sel
                    case 'O', mi_list=1; case 'X', mi_list=2; otherwise, mi_list=[1 2];
                end

                % --- Process each mode ---
                for mi = mi_list
                    N2g  = N2_modes{mi};
                    mch  = mode_labels{mi};
                    mask = isfinite(N2g) & real(N2g)>0;
                    if ~any(mask), continue; end

                    N_re   = real(sqrt(complex(N2g)));
                    ImN2g  = imag(N2g);
                    kpar_g = k_grid .* N_re .* sphi;
                    mask   = mask & (abs(kpar_g)>kpar_min);

                    idx_ok = find(mask);
                    if isempty(idx_ok), continue; end

                    for ii = idx_ok(:).'
                        w_ok    = w_grid(ii);
                        kp_ok   = kpar_g(ii);
                        N_ok    = N_re(ii);
                        N2_ok   = real(N2g(ii));      % for geometry and CSV
                        N2_full = N2g(ii);             % full complex for polarization
                        ImN2_ok = ImN2g(ii);

                        % --- Resonance curve bounds ---
                        % v_perp^2 >= 0 requires (Np^2+Y^2)*beta^2 - 2*Np*beta + (1-Y^2) = 0
                        % to have real roots.  Discriminant: Np^2 + Y^2 >= 1.
                        % Ref: Bornatici (1983) Sec. 2; Lopez (2025) Eq. 3.
                        Np = kp_ok*c/w_ok;
                        Y  = ell*omega_c/w_ok;
                        Dq = Np^2 + Y^2;
                        disc = Dq - 1;
                        if disc < 0, continue; end

                        sq = sqrt(disc);
                        b1 = (Np - Y*sq)/Dq;
                        b2 = (Np + Y*sq)/Dq;
                        if abs(b1)>=1 && abs(b2)>=1, continue; end

                        vp_lo = c*min(b1,b2);
                        vp_hi = c*max(b1,b2);
                        vp_lo = max(vp_lo, -K_mu*vth);
                        vp_hi = min(vp_hi,  K_mu*vth);
                        if vp_lo >= vp_hi, continue; end

                        % Skip negligible Maxwellian: if entire curve is
                        % beyond 5*v_th from origin, exp(-25) ~ 1e-11.
                        if vp_lo > 5*vth || vp_hi < -5*vth, continue; end

                        % --- J2D: mode-resolved coupling-weighted resonant measure ---
                        %
                        % Mode-resolved polarization-coupling factor based on
                        % harmonic Bessel structure (Bornatici 1983, Sec. 3;
                        % Stix 1992, Ch. 10, Eq. 10-46).
                        %
                        % Replaces generic W_ell = exp(-lam)*I_ell(lam) with the
                        % branch-specific coupling |Psi_ell^(mode)|^2,
                        % which depends on the wave polarization (O vs X) and
                        % the harmonic ell through ordinary Bessel functions
                        % evaluated at the actual v_perp on the resonance curve.
                        %
                        % Coupling function (Bornatici 1983, Sec. 3; Stix 1992, Ch. 10):
                        %   Psi_ell = eR*J_{ell-1}(z) + eL*J_{ell+1}(z) + ez*(vp/vr)*J_ell(z)
                        %   z = k_perp * v_perp / omega_c
                        %   eR, eL, ez: branch-consistent polarization eigenvector.
                        %
                        % Then: I(vp) = f_M * |Psi_ell|^2 * 2*pi*c^2/(ell*wc*gamma)
                        %
                        % At v_perp=0 endpoints:  J_ell(0)=0 for ell>=1, but
                        %   (vp/vr)*J_ell(z) has finite limit via J_ell(z) ~ z^ell/(2^ell*ell!).
                        %   For ell=1: limit = vp*k_perp/(2*wc).
                        %   For ell>=2: limit = 0.

                        vp_grid = linspace(vp_lo, vp_hi, N_vpar);

                        ell_wc = ell*omega_c;
                        ratio  = (w_ok - kp_ok.*vp_grid) ./ ell_wc;
                        vperp2 = c^2*(1 - ratio.^2) - vp_grid.^2;
                        vperp2 = max(vperp2, 0);
                        vperp  = sqrt(vperp2);

                        v2 = vp_grid.^2 + vperp2;
                        gamma_v = 1 ./ sqrt(max(1 - v2/c^2, 1e-30));
                        fM = (1/(pi^(1.5)*vth^3)) .* exp(-v2/vth^2);

                        % --- Branch-consistent polarization eigenvector ---
                        % Extracts (eR, eL, ez) from the wave matrix Lambda
                        % of the SAME model that produced the branch N^2.
                        %   kinetic  -> kinetic Lambda
                        %   warmfluid -> warm-fluid Lambda
                        %   booker/simple -> cold-plasma Lambda
                        % Ref: Stix (1992) Ch. 1, Sec. 1-4; Bornatici (1983) Sec. 3.
                        [eR, eL, ez_n] = branch_pol( ...
                            N2_full, w_ok, phi, omega_p, omega_c, vth, c, ...
                            refr_model, N_harm, nu_coll);

                        % --- Coupling function on resonance curve ---
                        % Psi_ell = eR*J_{ell-1}(z) + eL*J_{ell+1}(z) + ez*(vp/vr)*J_ell(z)
                        % z = k_perp*v_perp/omega_c
                        % Ref: Bornatici (1983) Sec. 3; Stix (1992) Ch. 10, Eq. 10-46.
                        Nperp = N_ok * cos(phi);
                        kperp_ok = (w_ok/c) * Nperp;
                        z = kperp_ok .* vperp / omega_c;

                        Jlm1 = besselj(ell-1, z);
                        Jlp1 = besselj(ell+1, z);
                        Jl   = besselj(ell, z);

                        % Stable vp/vr * J_ell(z) at v_perp=0
                        vp_Jl_over_vr = zeros(size(vp_grid));
                        nz = vperp > 1e-20;
                        vp_Jl_over_vr(nz) = vp_grid(nz) ./ vperp(nz) .* Jl(nz);
                        if ell == 1
                            ep = ~nz;
                            vp_Jl_over_vr(ep) = vp_grid(ep) * kperp_ok / (2*omega_c);
                        end

                        Psi = eR .* Jlm1 + eL .* Jlp1 + ez_n .* vp_Jl_over_vr;
                        W_mode = abs(Psi).^2;

                        % --- J2D: 1D resonance-curve integral ---
                        I_vp = fM .* W_mode .* (2*pi*c^2) ./ (ell_wc .* gamma_v);
                        J2D = n_r * trapz(vp_grid, I_vp);
                        if ~isfinite(J2D) || J2D < 0, J2D = 0; end

                        % --- n_res: 2D line-shape integral ---
                        %  Direct 2D velocity-space integral with Gaussian
                        %  resonance window.  Replaces the narrowband
                        %  approximation n_res = sqrt(2*pi)*|k_par|*v_th*J2D.
                        %  Saturation near n_r emerges naturally from the
                        %  finite velocity-space support of the Maxwellian;
                        %  no hard clipping to n_r is applied.
                        %  Ref: Bornatici (1983) Sec. 2.
                        sigma_w = abs(kp_ok) * vth;
                        vp_2d = linspace(-K_mu*vth, K_mu*vth, N_vp2d);
                        vr_2d = linspace(0, K_mu*vth, N_vr2d+1);
                        [VP2, VR2] = meshgrid(vp_2d, vr_2d);
                        v2_2d  = VP2.^2 + VR2.^2;
                        fM_2d  = (1/(pi^(1.5)*vth^3)) .* exp(-v2_2d/vth^2);
                        gam_2d = 1./sqrt(max(1-v2_2d/c^2, 1e-30));
                        Dw_2d  = w_ok - kp_ok.*VP2 - ell_wc./gam_2d;
                        G_2d   = exp(-Dw_2d.^2 ./ (2*sigma_w^2));
                        z_2d   = kperp_ok .* VR2 / omega_c;
                        Jlm1_2d = besselj(ell-1, z_2d);
                        Jlp1_2d = besselj(ell+1, z_2d);
                        Jl_2d   = besselj(ell, z_2d);
                        vpJl_vr_2d = zeros(size(VP2));
                        nz2 = VR2 > 1e-20;
                        vpJl_vr_2d(nz2) = VP2(nz2)./VR2(nz2).*Jl_2d(nz2);
                        if ell == 1
                            ez2 = ~nz2;
                            vpJl_vr_2d(ez2) = VP2(ez2)*kperp_ok/(2*omega_c);
                        end
                        Psi_2d = eR.*Jlm1_2d + eL.*Jlp1_2d + ez_n.*vpJl_vr_2d;
                        igr_2d = fM_2d .* abs(Psi_2d).^2 .* G_2d .* (2*pi) .* VR2;
                        n_res = n_r * trapz(vr_2d, trapz(vp_2d, igr_2d, 2));
                        if ~isfinite(n_res) || n_res < 0, n_res = 0; end
                        f_res = n_res / n_r;

                        % Dominant resonant electron (peak of 2D integrand)
                        [~, ipk2] = max(igr_2d(:));
                        [ir2, ip2] = ind2sub(size(igr_2d), ipk2);
                        v_par_dom  = vp_2d(ip2);
                        v_perp_dom = vr_2d(ir2);

                        % --- Output ---
                        root_count = root_count+1;
                        batch_ctr  = batch_ctr+1;

                        fprintf(fid, '%.3f,%.1f,%.3f,%.2f,%.1f,%.1f,', ...
                            r, rad2deg(theta), R_pos, B_0, rad2deg(phi), w_ok/(2*pi*1e9));
                        fprintf(fid, '%.4f,%.4f,%.6e,%.6e,%.4f,%.6f,%.6f,%.6e,%.6e,%s,%.6f\n', ...
                            v_par_dom/vth, v_perp_dom/vth, J2D, n_res, ...
                            w_ok/omega_c, N_ok, N2_ok, ImN2_ok, kp_ok, mch, f_res);

                        if mod(batch_ctr,2000)==0, fclose(fid); fid=fopen(filename,'a'); end

                        if print_live
                            row = {sprintf('%.3f',r), sprintf('%.0f',rad2deg(theta)), ...
                                   sprintf('%.3f',R_pos), sprintf('%.2f',B_0), ...
                                   sprintf('%.0f',rad2deg(phi)), sprintf('%.0f',w_ok/(2*pi*1e9)), ...
                                   sprintf('%.3f',v_par_dom/vth), sprintf('%.3f',v_perp_dom/vth), ...
                                   sprintf('%.2e',J2D), sprintf('%.2e',n_res), ...
                                   sprintf('%.3f',w_ok/omega_c), sprintf('%.4f',N_ok), ...
                                   sprintf('%.4f',N2_ok), sprintf('%.2e',ImN2_ok), ...
                                   sprintf('%.2e',kp_ok), mch, sprintf('%.4f',f_res)};
                            printRow(row, W);
                        end
                    end
                end  % mode
                clear mi_list;
            end  % phi
        end  % B0
    end  % theta
end  % r

fclose(fid);
fprintf('\nRoots: %d  |  File: %s\n', root_count, filename);


%% ==================== DISPERSION FUNCTIONS ==================== %%

function [N2_O, N2_X] = n_booker(w, theta, wp, wc, nu, wms)
    w_size=size(w); w=w(:); tiny=1e-14;
    ws=w; ws(abs(ws)<tiny)=tiny;
    we=ws+1i*nu;
    L=1-wp^2./(ws.*(we-wc)); R=1-wp^2./(ws.*(we+wc));
    S=.5*(R+L); P=1-wp^2./(ws.*we);
    st2=sin(theta)^2; ct2=cos(theta)^2;
    Aq=S.*ct2+P.*st2; Bq=(R.*L).*ct2+(P.*S).*(1+st2); Cq=P.*R.*L;
    disc=Bq.^2-4.*Aq.*Cq;
    if nu>0, F=sqrt(disc);
    else, disc(disc<0&disc>-1e-12.*max(1,abs(Bq).^2))=0; F=sqrt(max(real(disc),0)); end
    bad=abs(Aq)<1e-12;
    np=(Bq+F)./(2.*Aq); nm=(Bq-F)./(2.*Aq); np(bad)=NaN; nm(bad)=NaN;
    N2Op=real(P); S_r=real(S); RLS=NaN(size(w));
    g=abs(S_r)>1e-14; RLS(g)=real(R(g).*L(g))./S_r(g);
    N2_O=NaN(size(w)); N2_X=NaN(size(w));
    for j=1:numel(w)
        for val=[np(j),nm(j)]
            if ~isfinite(val), continue; end
            rv=real(val); if rv<=0, continue; end
            if nu==0 && abs(imag(val))>1e-10*max(1,abs(rv)), continue; end
            if abs(rv-N2Op(j))<=abs(rv-RLS(j)), N2_O(j)=rv; else, N2_X(j)=rv; end
        end
    end
    if strcmp(wms,'O'), N2_X(:)=NaN; elseif strcmp(wms,'X'), N2_O(:)=NaN; end
    N2_O=reshape(N2_O,w_size); N2_X=reshape(N2_X,w_size);
end


function [N2_O, N2_X] = n_booker_cold(w, theta, wp, wc)
    w_size=size(w); w=w(:); ws=w; ws(abs(ws)<1e-14)=1e-14;
    L=1-wp^2./(ws.*(ws-wc)); R=1-wp^2./(ws.*(ws+wc));
    S=.5*(R+L); P=1-wp^2./ws.^2;
    st2=sin(theta)^2; ct2=cos(theta)^2;
    Aq=S.*ct2+P.*st2; Bq=(R.*L).*ct2+(P.*S).*(1+st2); Cq=P.*R.*L;
    disc=Bq.^2-4.*Aq.*Cq; disc(disc<0&disc>-1e-12.*max(1,abs(Bq).^2))=0;
    F=sqrt(max(real(disc),0));
    np=(Bq+F)./(2.*Aq); nm=(Bq-F)./(2.*Aq);
    bad=abs(Aq)<1e-12; np(bad)=NaN; nm(bad)=NaN;
    isRP=@(x)isfinite(x)&real(x)>0&abs(imag(x))<1e-10*max(1,abs(real(x)));
    np(~isRP(np))=NaN; nm(~isRP(nm))=NaN; np=real(np); nm=real(nm);
    N2Op=real(P); S_r=real(S); RLS=NaN(size(w));
    g=abs(S_r)>1e-14; RLS(g)=real(R(g).*L(g))./S_r(g);
    N2_O=NaN(size(w)); N2_X=NaN(size(w));
    for j=1:numel(w)
        for val=[np(j),nm(j)]
            if ~isfinite(val)||val<=0, continue; end
            dO=abs(val-N2Op(j)); dX=abs(val-RLS(j)); if ~isfinite(dX),dX=inf;end
            if dO<=dX, N2_O(j)=val; else, N2_X(j)=val; end
        end
    end
    N2_O=reshape(N2_O,w_size); N2_X=reshape(N2_X,w_size);
end


function [N2_O, N2_X] = n_warm_fluid(w, theta, wp, wc, vth, cl, wms)
    w_size=size(w); w=w(:); Nw=numel(w);
    vt2=vth^2; wp2=wp^2; wc2=wc^2;
    st=sin(theta); ct=cos(theta); st2=st^2; ct2=ct^2;
    [O0,X0]=n_booker_cold(w,theta,wp,wc);
    N2_O=NaN(Nw,1); N2_X=NaN(Nw,1);
    switch wms, case 'O',sb=1; case 'X',sb=2; otherwise,sb=[1 2]; end
    for j=1:Nw
        wj=w(j); if wj<=0,continue;end; wj2=wj^2;
        seeds=[O0(j),X0(j)];
        for ib=sb
            g=seeds(ib); if ~isfinite(g)||g<=0,continue;end
            N2i=g;
            for it=1:15
                if ~isfinite(N2i)||N2i<=0,N2i=NaN;break;end
                Dv=wf_det(N2i,wj,wj2,wp2,wc,wc2,vt2,cl,st,ct,st2,ct2);
                dN=max(1e-6,abs(N2i)*1e-5);
                Dp=wf_det(N2i+dN,wj,wj2,wp2,wc,wc2,vt2,cl,st,ct,st2,ct2);
                dd=(Dp-Dv)/dN; if abs(dd)<1e-30,break;end
                del=-real(Dv)/real(dd);
                al=1; if abs(del)>.5*abs(N2i), al=.5*abs(N2i)/abs(del); end
                N2n=N2i+al*del; if real(N2n)<=0, N2n=N2i*.5; end; N2n=real(N2n);
                if abs(N2n-N2i)/max(1,abs(N2i))<1e-8, N2i=N2n; break; end
                N2i=N2n;
            end
            if isfinite(N2i)&&N2i>0
                if real(N2i)>max(3*g,5), continue; end  % Bernstein guard
                if ib==1, N2_O(j)=N2i; else, N2_X(j)=N2i; end
            end
        end
    end
    N2_O=reshape(N2_O,w_size); N2_X=reshape(N2_X,w_size);
end


function D = wf_det(N2,wj,wj2,wp2,wc,wc2,vt2,cl,st,ct,st2,ct2)
    Nv=sqrt(N2); km=Nv*wj/cl; k2v=(km^2)*vt2;
    Den=wj2*(wj2-k2v)-wc2*(wj2-k2v*st2);
    if abs(Den)<1e-30*wj2^2, D=NaN; return; end
    A=wj2-k2v*st2; B=k2v*st*ct; C=wj2-wc2-k2v*ct2;
    exx=1-wp2*A/Den; eyy=1-wp2*(wj2-k2v)/Den;
    exy=1i*(wc/wj)*wp2*A/Den; exz=-1i*wp2*B/Den;
    eyz=-1i*(wc/wj)*wp2*B/Den; ezz=1-wp2*C/Den;
    Np=Nv*st; Nr=Nv*ct;
    L11=exx-Np^2; L12=exy; L13=exz+Np*Nr;
    L22=eyy-N2; L23=eyz; L33=ezz-Nr^2;
    D=L11*(L22*L33-L23*(-eyz))-L12*((-exy)*L33-L23*L13)+L13*((-exy)*(-eyz)-L22*L13);
end


function [N2_O, N2_X] = n_warm_plasma(w, theta, wp, wc, vth, cl, Nh, wms)
% Kinetic dielectric (Harris/Stix) with O/X tracking.
% Ref: Stix (1992) Ch. 10; Harris (1961); Weideman (1994).
    w_size=size(w); w=w(:); Nw=numel(w);
    [O0,X0]=n_booker_cold(w,theta,wp,wc);
    N2_O=NaN(Nw,1); N2_X=NaN(Nw,1);
    st=sin(theta); ct=cos(theta);
    switch wms, case 'O',sb=1; case 'X',sb=2; otherwise,sb=[1 2]; end
    for j=1:Nw
        wj=w(j); if wj<=0,continue;end
        seeds=[O0(j),X0(j)];
        for ib=sb
            g=seeds(ib); if ~isfinite(g)||g<=0,continue;end
            N2i=complex(g);
            for it=1:12
                if ~isfinite(N2i)||real(N2i)<=0,N2i=NaN;break;end
                Ni=sqrt(N2i); kp=Ni*st*wj/cl; kr=Ni*ct*wj/cl;
                lam=(real(kr)*vth)^2/(2*wc^2);
                [K11,K12,K13,K22,K23,K33]=warm_tensor(wj,kp,lam,wp,wc,vth,Nh);
                Np2=N2i; Npar=Ni*st; Nperp=Ni*ct;
                L11=K11-Npar^2; L12=K12; L13=K13+Npar*Nperp;
                L22=K22-Np2; L23=K23; L33=K33-Nperp^2;
                Dv=L11*(L22*L33+L23^2)+L12*(L12*L33+L13*L23)+L13*(L12*L23-L13*L22);
                dN=max(1e-6,abs(N2i)*1e-5); N2p=N2i+dN;
                Np=sqrt(N2p); kp2=Np*st*wj/cl; kr2=Np*ct*wj/cl;
                lp=(real(kr2)*vth)^2/(2*wc^2);
                [P11,P12,P13,P22,P23,P33]=warm_tensor(wj,kp2,lp,wp,wc,vth,Nh);
                Np_par=Np*st; Np_perp=Np*ct;
                M11=P11-Np_par^2; M12=P12; M13=P13+Np_par*Np_perp;
                M22=P22-N2p; M23=P23; M33=P33-Np_perp^2;
                Dp=M11*(M22*M33+M23^2)+M12*(M12*M33+M13*M23)+M13*(M12*M23-M13*M22);
                dd=(Dp-Dv)/dN; if abs(dd)<1e-30,break;end
                del=-Dv/dd;
                al=1; if abs(del)>.5*abs(N2i), al=.5*abs(N2i)/abs(del); end
                N2n=N2i+al*del;
                if real(N2n)<=0, N2n=abs(real(N2i))*.5+1i*imag(N2n); end
                if abs(N2n-N2i)/max(1,abs(N2i))<1e-8, N2i=N2n; break; end
                N2i=N2n;
            end
            if isfinite(N2i)&&real(N2i)>0
                if real(N2i)>max(3*g,5), continue; end  % Bernstein guard
                if ib==1, N2_O(j)=N2i; else, N2_X(j)=N2i; end
            end
        end
    end
    N2_O=reshape(N2_O,w_size); N2_X=reshape(N2_X,w_size);
end


function [K11,K12,K13,K22,K23,K33] = warm_tensor(w,kp,lam,wp,wc,vth,Nh)
% Kinetic dielectric tensor, |n| <= Nh.
% Ref: Stix (1992) Ch. 10; Fitzpatrick (2015) Sec. 8.8.
    kps=kp; if abs(kps)<1e-10, kps=1e-10*sign(real(kps)+1e-30); end
    lam=max(real(lam),1e-20); eml=exp(-lam);
    den=kps*vth; pref=wp^2/(w*den);
    In=zeros(1,Nh+1); Inp=zeros(1,Nh+1);
    for nn=0:Nh
        In(nn+1)=besseli(nn,lam);
        if nn==0, Inp(nn+1)=besseli(1,lam);
        else, Inp(nn+1)=besseli(nn-1,lam)-(nn/lam)*In(nn+1); end
    end
    S11=0;S12=0;S13=0;S22=0;S23=0;S33=0;
    for n=-Nh:Nh
        an=abs(n); I_n=In(an+1); I_np=Inp(an+1);
        Gn=eml*I_n; dn=eml*(I_np-I_n);
        xi=(w-n*wc)/den; Z=1i*sqrt(pi)*faddeeva(xi); Zp=-2*(1+xi*Z);
        if n==0, t11=0; else, t11=n^2*Gn/lam*Z; end
        t12=1i*n*dn*Z;
        t22=(n^2*Gn/lam+2*lam*dn)*Z;
        t33=Gn*xi*Zp;
        if lam>1e-15, t13=n/sqrt(2*lam)*Gn*Zp; t23=1i*n*sqrt(lam/2)*dn*Zp;
        else, t13=0; t23=0; end
        S11=S11+t11; S12=S12+t12; S13=S13+t13;
        S22=S22+t22; S23=S23+t23; S33=S33+t33;
    end
    K11=1+pref*S11; K12=pref*S12; K13=pref*S13;
    K22=1+pref*S22; K23=pref*S23; K33=1+pref*S33;
end


%% ==================== POLARIZATION ==================== %%

function [eR, eL, ez] = branch_pol(N2, w, phi, wp, wc, vth, cl, model, Nh, nu)
%BRANCH_POL  Branch-consistent polarization eigenvector.
%
%  Builds the wave matrix Lambda from the SAME dielectric model that
%  produced the branch N^2, extracts the null eigenvector E, and
%  converts to normalized circular components (eR, eL, ez).
%
%  The wave matrix uses the exact structure from the branch determinant
%  solver: L21=-K12, L31=K13+NpNr, L32=-K23 (not Hermitian conjugates).
%  For the kinetic branch, N^2 is the full complex value.
%  For the Booker branch with nu>0, collision-regularized Stix parameters
%  are used (w_eff = w + i*nu), matching the Booker dispersion solve.
%
%    kinetic   -> kinetic Lambda (Harris/Stix tensor)
%    warmfluid -> warm-fluid Lambda (pressure tensor)
%    booker    -> cold-plasma Lambda (Stix S,D,P; regularized if nu>0)
%
%  Convention: exp(-iwt), eR = (Ex-iEy)/sqrt(2), eL = (Ex+iEy)/sqrt(2).
%  |eR|^2 + |eL|^2 + |ez|^2 = 1.
%
%  Ref: Stix, Waves in Plasmas (1992), Sec. 1-4;
%       Bornatici et al. (1983), Sec. 3;
%       Brambilla, Phys. Plasmas 2, 1094 (1995) [collision reg.].

    N = sqrt(N2);
    st = sin(phi);  ct = cos(phi);
    Npar = N*st;  Nperp = N*ct;

    % Build dielectric tensor K_ij for the selected model
    switch model
        case 'kinetic'
            kpar  = Npar*w/cl;
            kperp = Nperp*w/cl;
            lam   = (real(kperp)*vth)^2 / (2*wc^2);
            [K11,K12,K13,K22,K23,K33] = warm_tensor(w, kpar, lam, wp, wc, vth, Nh);

        case 'warmfluid'
            [K11,K12,K13,K22,K23,K33] = wf_tensor_elems(N2, w, wp, wc, vth, cl, st, ct);

        otherwise  % cold (booker), with collision regularization if nu>0
            we = w + 1i*nu;
            Lc = 1 - wp^2/(w*(we-wc));
            Rc = 1 - wp^2/(w*(we+wc));
            S  = .5*(Rc+Lc);  D = .5*(Rc-Lc);  P = 1-wp^2/(w*we);
            K11 = S;  K12 = -1i*D;  K13 = 0;
            K22 = S;  K23 = 0;      K33 = P;
    end

    % Wave matrix Lambda: SAME structure as the branch determinant solver.
    %   L21 = -K12,  L31 = K13+NpNr,  L32 = -K23
    % NOT Hermitian conjugates (which differ for complex kinetic K_ij).
    % Ref: Stix (1992) Sec. 1-4.
    L11 = K11 - Npar^2;
    L12 = K12;
    L13 = K13 + Npar*Nperp;
    L21 = -K12;
    L22 = K22 - N2;
    L23 = K23;
    L31 = K13 + Npar*Nperp;
    L32 = -K23;
    L33 = K33 - Nperp^2;

    % Null eigenvector via cross products of rows
    r1 = [L11, L12, L13];
    r2 = [L21, L22, L23];
    r3 = [L31, L32, L33];

    v1 = cross(r1, r2);
    v2 = cross(r1, r3);
    v3 = cross(r2, r3);

    n1 = abs(v1(1))^2+abs(v1(2))^2+abs(v1(3))^2;
    n2 = abs(v2(1))^2+abs(v2(2))^2+abs(v2(3))^2;
    n3 = abs(v3(1))^2+abs(v3(2))^2+abs(v3(3))^2;

    [~, idx] = max([n1, n2, n3]);
    vecs = {v1, v2, v3};
    E = vecs{idx};

    en = sqrt(abs(E(1))^2 + abs(E(2))^2 + abs(E(3))^2);
    if en < 1e-30
        eR = 1/sqrt(2);  eL = 1/sqrt(2);  ez = 0;
        return;
    end
    E = E / en;

    % Circular components: eR = (Ex-iEy)/sqrt(2), eL = (Ex+iEy)/sqrt(2)
    eR = (E(1) - 1i*E(2)) / sqrt(2);
    eL = (E(1) + 1i*E(2)) / sqrt(2);
    ez = E(3);
end


function [K11,K12,K13,K22,K23,K33] = wf_tensor_elems(N2, w, wp, wc, vth, cl, st, ct)
%WF_TENSOR_ELEMS  Warm-fluid dielectric tensor elements.
%
%  Same physics as wf_det, but returns individual elements for
%  eigenvector extraction.  Angle convention: st=sin(phi)=cos(theta_kB).
%
%  Ref: Stix (1992) Ch. 2; Swanson (2003) Ch. 2-4.

    wp2 = wp^2;  wc2 = wc^2;  vt2 = vth^2;
    w2 = w^2;  st2 = st^2;  ct2 = ct^2;
    km = sqrt(N2)*w/cl;  k2v = km^2 * vt2;

    Den = w2*(w2-k2v) - wc2*(w2-k2v*st2);
    if abs(Den) < 1e-30*w2^2
        K11=1; K12=0; K13=0; K22=1; K23=0; K33=1;
        return;
    end

    A = w2 - k2v*st2;
    B = k2v*st*ct;
    C = w2 - wc2 - k2v*ct2;

    K11 = 1 - wp2*A/Den;
    K22 = 1 - wp2*(w2-k2v)/Den;
    K12 = 1i*(wc/w)*wp2*A/Den;
    K13 = -1i*wp2*B/Den;
    K23 = -1i*(wc/w)*wp2*B/Den;
    K33 = 1 - wp2*C/Den;
end


%% ==================== UTILITY FUNCTIONS ==================== %%

function n_r = ne_pedestal(r, a, nc, ne, nped, rped, dc, dp)
    rho=r/a; rho0=min(max(rped/a,.85),.99);
    dc=max(dc,1e-8); dp=max(dp,1e-8);
    if rho<=rho0
        n_r = nped + (nc-nped)*(exp(-dc*rho/rho0)-exp(-dc))/(1-exp(-dc));
    else
        x=(rho-rho0)/max(1-rho0,eps);
        n_r = nped + (nped-ne)/(1-exp(-dp))*(exp(-dp*x)-1);
    end
    if ~isfinite(n_r)||n_r<=0, n_r=max(ne,eps); end
end


function printRow(f, W)
    for i=1:numel(f), fprintf('%s',cpad(f{i},W(i))); end; fprintf('\n');
end
function o = cpad(s, W)
    s=char(string(s)); L=length(s);
    if L>=W, o=s(1:W); return; end
    p=W-L; lf=floor(p/2); o=[repmat(' ',1,lf),s,repmat(' ',1,p-lf)];
end


function Z = plasma_Z(xi)
    Z = 1i*sqrt(pi)*faddeeva(xi);
end

function w = faddeeva(z)
    if numel(z)>1, w=zeros(size(z)); for i=1:numel(z), w(i)=faddeeva(z(i)); end; return; end
    x=real(z); y=imag(z);
    if y>4.29||(y>=.195*abs(x)-.176&&abs(x)<6.2)
        if abs(z)<.1, z2=z^2; w=1-(2/sqrt(pi))*z*(1-z2/1.5+z2^2/3.75-z2^3/8.75);
        else, w=cf_w(z); end
    elseif abs(x)<5.5, w=hum_w4(x,y);
    else, z2=z^2; w=(1i/sqrt(pi))*(1/z)*(1-.5/z2+.75/z2^2-1.875/z2^3); 
    end
    if y<0, w=2*exp(-z^2)-conj(w); end
end

function w = cf_w(z)
    tiny=1e-30; a1=1i/sqrt(pi); f=z; C=z; D=0;
    for n=1:100
        if n==1,an=a1;else,an=-(n-1)*.5;end
        D=z+an*D; if abs(D)<tiny,D=tiny;end
        C=z+an/C; if abs(C)<tiny,C=tiny;end
        D=1/D; d=C*D; f=f*d; if abs(d-1)<1e-12,break;end
    end; w=f;
end

function w = hum_w4(x, y)
    s=abs(x)+y; t=y-1i*x;
    if s>=15, w=t*.5641896/(.5+t^2);
    elseif s>=5.5, u=t*t; w=t*(1.410474+u*.5641896)/(.75+u*(3+u));
    else, u=t*t; v=u*u;
        w=exp(u)-t*(36183.31-v*(3321.99-v*(1540.787-v*(219.031-...
          v*(35.7668-v*(1.320522-v*.56419))))))/(32066.6-v*(24322.8-...
          v*(9022.23-v*(2186.18-v*(364.219-v*(61.5704-v*(1.84144-v)))))));
    end
end