function out = ct_simulate(p)
% CT_SIMULATE  Polychromatic CT simulation engine.
%   out = ct_simulate(p) where p is a parameter struct. Missing fields take
%   defaults from ct_defaults(). Every field is a candidate GUI control.
%
%   Returns struct with fields:
%     .recon      reconstructed image (1/cm)
%     .sinogram   effective sinogram used for reconstruction
%     .mu_centre .mu_edge .cupping   region metrics
%     .mu_detail .contrast           detail metrics (NaN if no detail)
%     .noise_sd                      background standard deviation
%     .resid                         SIRT residual history ([] for FBP)

    d = ct_defaults();
    fn = fieldnames(d);
    for k = 1:numel(fn)                       % fill unspecified fields
        if ~isfield(p, fn{k}), p.(fn{k}) = d.(fn{k}); end
    end
    pix_cm = p.obj_diam_cm / p.circ_pix;      % cm per pixel

    % ---- materials ---------------------------------------------------
    [Eo, mo] = load_material(p.obj_material);
    if isempty(p.det_material)
        Ed = Eo; md = mo;
    else
        [Ed, md] = load_material(p.det_material);
    end

    % ---- spectrum ----------------------------------------------------
    S = readmatrix(sprintf('spectrum_%dkvp.csv', p.kvp));
    S = S(~any(isnan(S),2), :);
    energies = S(:,1);  fluence = S(:,2);

    % ---- phantom -----------------------------------------------------
    [xx, yy] = meshgrid(0:p.grid-1, 0:p.grid-1);
    c = (p.grid-1)/2;
    obj_mask = ((xx-c).^2 + (yy-c).^2) <= (p.circ_pix/2)^2;

    if isempty(p.det_material) || p.det_diam_cm <= 0
        det_mask = false(size(obj_mask));  dc = [c c];
    else
        det_r = (p.det_diam_cm / p.obj_diam_cm) * p.circ_pix / 2;
        if strcmpi(p.det_position, 'centre'), dc = [c, c];
        else,                                 dc = [c + 0.7*p.circ_pix/2, c];
        end
        det_mask = (((xx-dc(1)).^2 + (yy-dc(2)).^2) <= det_r^2) & obj_mask;
    end
    obj_only = obj_mask & ~det_mask;

    % ---- geometry: path lengths --------------------------------------
    switch lower(p.geometry)
      case 'parallel'
        theta = 0 : 180/p.n_angles : 180-180/p.n_angles;
        path_o = radon(double(obj_only), theta) * pix_cm;
        path_d = radon(double(det_mask), theta) * pix_cm;
      case 'fan'
        fanargs = {'FanSensorGeometry','line', ...
                   'FanSensorSpacing', p.det_pixel_pix, ...
                   'FanRotationIncrement', p.rot_increment};
        path_o = fanbeam(double(obj_only), p.source_dist_pix, fanargs{:}) * pix_cm;
        path_d = fanbeam(double(det_mask), p.source_dist_pix, fanargs{:}) * pix_cm;
        theta = 0 : p.rot_increment : 360-p.rot_increment;   % for the parallel-approx path
      otherwise
        error('geometry must be ''parallel'' or ''fan''');
    end

    % ---- polychromatic forward model ---------------------------------
    sig = zeros(size(path_o));  sig0 = 0;
    for k = 1:numel(energies)
        E = energies(k);  N0 = fluence(k);
        mu_o = exp(interp1(log(Eo), log(mo), log(E)));
        mu_d = exp(interp1(log(Ed), log(md), log(E)));
        T = exp(-(mu_o*path_o + mu_d*path_d));
        switch lower(p.detector)
          case 'counting',    w = N0;
          case 'integrating', w = E*N0;
          otherwise, error('detector must be ''counting'' or ''integrating''');
        end
        sig = sig + w*T;   sig0 = sig0 + w;
    end

    % ---- incident statistics ------------------------------------------
    if isfinite(p.n0_per_pixel) && p.n0_per_pixel > 0
        s = p.n0_per_pixel / sig0;
        sig = poissrnd(sig * s) / s;
    end
    sino = -log(sig / sig0);

    % ---- reconstruction ------------------------------------------------
    resid = [];
    switch lower(p.recon_method)
      case 'fbp'
        if strcmpi(p.geometry, 'fan')
            rec = ifanbeam(sino, p.source_dist_pix, fanargs{:}, ...
                           'OutputSize', p.grid, 'Filter', p.recon_filter);
        else
            rec = iradon(sino, theta, p.recon_filter, p.grid);
        end
      case 'fbp_parallel_approx'   % deliberately wrong: fan data, parallel recon
        rec = iradon(sino, theta, p.recon_filter, p.grid);
      case 'sirt'
        [rec, resid] = sirt_recon(sino, theta, p.grid, p.n_iter, p.relax);
      otherwise
        error('recon_method must be fbp | fbp_parallel_approx | sirt');
    end
    rec = rec / pix_cm;

    % ---- metrics ---------------------------------------------------------
    r = sqrt((xx-c).^2 + (yy-c).^2);
    centre = r < 0.15*(p.circ_pix/2);
    edge   = (r > 0.80*(p.circ_pix/2)) & (r < 0.92*(p.circ_pix/2));
    bg     = r > 1.15*(p.circ_pix/2);

    out.recon    = rec;
    out.sinogram = sino;
    out.params   = p;
    out.resid    = resid;
    out.mu_centre = mean(rec(centre));
    out.mu_edge   = mean(rec(edge));
    out.cupping   = 100*(out.mu_edge - out.mu_centre)/out.mu_edge;
    out.noise_sd  = std(rec(bg));
    if any(det_mask(:))
        dz = ((xx-dc(1)).^2 + (yy-dc(2)).^2) <= (0.6*det_r)^2;
        bgz = obj_only & (r < 0.5*p.circ_pix/2);
        out.mu_detail = mean(rec(dz));
        out.contrast  = 100*(out.mu_detail - mean(rec(bgz)))/mean(rec(bgz));
    else
        out.mu_detail = NaN;  out.contrast = NaN;
    end
end

% ======================================================================
function d = ct_defaults()
    d.kvp             = 80;
    d.obj_material    = 'water';
    d.obj_diam_cm     = 10.0;
    d.det_material    = '';        % '' = no detail
    d.det_diam_cm     = 2.0;
    d.det_position    = 'extremity';
    d.geometry        = 'parallel';   % 'parallel' | 'fan'
    d.source_dist_pix = 500;
    d.det_pixel_pix   = 1.0;
    d.rot_increment   = 1.0;
    d.detector        = 'integrating';  % 'counting' | 'integrating'
    d.n_angles        = 180;
    d.n0_per_pixel    = Inf;       % Inf = noiseless
    d.recon_method    = 'fbp';     % 'fbp' | 'fbp_parallel_approx' | 'sirt'
    d.recon_filter    = 'Hann';
    d.n_iter          = 100;
    d.relax           = 1.0;
    d.grid            = 180;
    d.circ_pix        = 160;
end

% ======================================================================
function [rec, resid] = sirt_recon(sino, theta, grid, n_iter, relax)
    norm_img  = iradon(ones(size(sino)), theta, 'none', grid);
    norm_img(norm_img < 1e-6) = 1e-6;
    norm_sino = radon(ones(grid), theta);
    norm_sino(norm_sino < 1e-6) = 1e-6;
    rec = zeros(grid);  resid = zeros(1, n_iter);
    for it = 1:n_iter
        fwd = radon(rec, theta);
        rec = rec + relax * (iradon((sino - fwd)./norm_sino, theta, 'none', grid) ./ norm_img);
        rec(rec < 0) = 0;
        resid(it) = norm(sino - radon(rec, theta), 'fro');
    end
end

% ======================================================================
function [E_keV, mu_cm] = load_material(name)
    switch lower(name)
        case 'water', f = 'water_mu.txt'; rho = 1.00;
        case 'pmma',  f = 'pmma_mu.txt';  rho = 1.19;
        case 'al',    f = 'al_mu.txt';    rho = 2.70;
        otherwise, error('Unknown material: %s', name);
    end
    M = readmatrix(f, 'FileType','text');
    M = M(~any(isnan(M),2), :);
    E_keV = M(:,1)*1000;  mu_cm = M(:,2)*rho;
    keep = [true; diff(E_keV) > 0];
    E_keV = E_keV(keep);  mu_cm = mu_cm(keep);
end