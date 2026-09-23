% fig_sirt_sparse.m - Figure 4.4: FBP vs SIRT at 20 projection angles
% Same noisy sinogram reconstructed both ways.

p = struct('kvp',80, 'det_material','', 'n_angles',20, ...
           'n0_per_pixel',1e5, 'recon_filter','Hann');

rng(1);                              % fixed seed: reproducible noise
p.recon_method = 'fbp';
r_fbp  = ct_simulate(p);

rng(1);                              % same noise realisation
p.recon_method = 'sirt';
p.n_iter = 100;
r_sirt = ct_simulate(p);

fprintf('FBP  : cupping=%.2f%%  noise=%.4f  min mu=%.4f\n', ...
        r_fbp.cupping,  r_fbp.noise_sd,  min(r_fbp.recon(:)));
fprintf('SIRT : cupping=%.2f%%  noise=%.4f  min mu=%.4f\n', ...
        r_sirt.cupping, r_sirt.noise_sd, min(r_sirt.recon(:)));

figure('Position',[100 100 900 420]);
cl = [-0.05 0.30];                   % shared colour scale for fair comparison

subplot(1,2,1);
imagesc(r_fbp.recon, cl);  axis image off; colormap gray; colorbar;
title(sprintf('(a) FBP, 20 angles  (noise sd = %.4f)', r_fbp.noise_sd));

subplot(1,2,2);
imagesc(r_sirt.recon, cl); axis image off; colorbar;
title(sprintf('(b) SIRT, 100 iterations  (noise sd = %.4f)', r_sirt.noise_sd));

exportgraphics(gcf, 'fig_sirt_sparse.png', 'Resolution', 200);
disp('Saved fig_sirt_sparse.png');