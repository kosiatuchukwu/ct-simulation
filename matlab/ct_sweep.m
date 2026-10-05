% ct_sweep.m - controlled parameter sweep for Results section 4.2
% One parameter varied at a time from a fixed baseline.
clear; clc;

base = struct('kvp',80, 'obj_material','water', 'obj_diam_cm',10, ...
              'det_material','pmma', 'det_diam_cm',2, 'det_position','centre', ...
              'geometry','parallel', 'n_angles',180, 'detector','integrating', ...
              'n0_per_pixel',1e6, 'recon_method','fbp', 'recon_filter','Hann');

fprintf('BASELINE: 80 kVp, water + PMMA detail (centre), parallel,\n');
fprintf('180 angles, integrating, 1e6 photons/pixel, Hann filter\n\n');

% ---- arm 1: reconstruction filter ------------------------------------
fprintf('--- Filter ---\n');
fprintf('%-14s %10s %10s %10s\n','filter','mu_centre','cupping','noise_sd');
for f = {'Ram-Lak','Shepp-Logan','Cosine','Hamming','Hann'}
    p = base;  p.det_material = '';  p.recon_filter = f{1};
    r = ct_simulate(p);
    fprintf('%-14s %10.4f %9.2f%% %10.4f\n', f{1}, r.mu_centre, r.cupping, r.noise_sd);
end

% ---- arm 2: detail material -------------------------------------------
fprintf('\n--- Detail material ---\n');
fprintf('%-8s %10s %10s %10s %10s\n','material','mu_detail','contrast','cupping','noise_sd');
for m = {'','pmma','al'}
    p = base;  p.det_material = m{1};
    r = ct_simulate(p);
    nm = m{1}; if isempty(nm), nm = 'none'; end
    fprintf('%-8s %10.4f %9.1f%% %9.2f%% %10.4f\n', ...
            nm, r.mu_detail, r.contrast, r.cupping, r.noise_sd);
end

% ---- arm 3: incident statistics ---------------------------------------
fprintf('\n--- Incident statistics ---\n');
fprintf('%-10s %10s %10s %10s\n','photons','contrast','cupping','noise_sd');
for n = [1e4 1e5 1e6 1e7]
    p = base;  p.det_material = '';  p.n0_per_pixel = n;
    r = ct_simulate(p);
    fprintf('%-10.0e %9.1f%% %9.2f%% %10.4f\n', n, r.contrast, r.cupping, r.noise_sd);
end

% ---- fit the noise-vs-statistics power law -----------------------------
N   = [1e4 1e5 1e6 1e7];
sd  = zeros(size(N));
for i = 1:numel(N)
    p = base;  p.det_material = '';  p.n0_per_pixel = N(i);
    r = ct_simulate(p);
    sd(i) = r.noise_sd;
end
cf = polyfit(log(N), log(sd), 1);
fprintf('\nnoise ~ N^%.3f   (theory: -0.5)\n', cf(1));

figure;
loglog(N, sd, 'o-', 'LineWidth', 1.5); hold on;
loglog(N, exp(cf(2))*N.^cf(1), 'k--');
loglog(N, sd(1)*(N/N(1)).^-0.5, 'r:');
legend('measured', sprintf('fit: N^{%.2f}', cf(1)), 'theory: N^{-0.5}');
xlabel('incident photons per pixel'); ylabel('background noise (sd)');
grid on; title('Noise versus incident statistics');