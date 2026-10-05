% ct_check.m - does the consolidated engine reproduce the known benchmarks?
clear; clc;

fprintf('--- benchmark: parallel, 80 kVp, noiseless ---\n');
p = struct('kvp',80,'detector','integrating');
r = ct_simulate(p);
fprintf('cupping = %.2f%%   (expect 4.43)\n', r.cupping);

fprintf('\n--- fan-beam, D=500 ---\n');
p = struct('kvp',80,'geometry','fan','source_dist_pix',500);
r = ct_simulate(p);
fprintf('cupping = %.2f%%   (expect 4.43)\n', r.cupping);

fprintf('\n--- same fan data, parallel approximation ---\n');
p.recon_method = 'fbp_parallel_approx';
r = ct_simulate(p);
fprintf('cupping = %.2f%%   (expect 2.94)\n', r.cupping);

fprintf('\n--- SIRT, 60 angles, 1e5 photons ---\n');
p = struct('kvp',80,'n_angles',60,'n0_per_pixel',1e5,'recon_method','sirt');
r = ct_simulate(p);
fprintf('cupping = %.2f%%  noise = %.4f   (expect ~4.43, ~0.0012)\n', ...
        r.cupping, r.noise_sd);