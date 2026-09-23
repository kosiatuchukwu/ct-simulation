% ct_sirt_demo.m - hand-implemented SIRT, parallel-beam
% Compares FBP with iterative reconstruction on the same noisy sinogram.

KVP        = 80;
N_ITER     = 100;      % SIRT iterations
RELAX      = 1.0;      % relaxation factor (step size)
N_ANGLES   = 20;       % deliberately sparse: where iterative should shine
N0_PER_PIXEL = 1e5;    % incident photons per pixel (noise level)

GRID = 180; CIRC_PIX = 160; OBJ_DIAM_CM = 10.0;
PIXEL_SIZE_CM = OBJ_DIAM_CM / CIRC_PIX;

% ---- water mu(E) ----------------------------------------------------
M = readmatrix('water_mu.txt','FileType','text');
M = M(~any(isnan(M),2), :);
E_tab = M(:,1)*1000;  mu_tab = M(:,2)*1.0;

% ---- spectrum -------------------------------------------------------
S = readmatrix(sprintf('spectrum_%dkvp.csv', KVP));
S = S(~any(isnan(S),2), :);
energies = S(:,1);  fluence = S(:,2);

% ---- phantom --------------------------------------------------------
[xx, yy] = meshgrid(0:GRID-1, 0:GRID-1);
c = (GRID-1)/2;
mask = double(((xx-c).^2 + (yy-c).^2) <= (CIRC_PIX/2)^2);

theta = 0 : 180/N_ANGLES : 180-180/N_ANGLES;
path_cm = radon(mask, theta) * PIXEL_SIZE_CM;

% ---- polychromatic forward model, energy-integrating ----------------
sig = zeros(size(path_cm));  sig0 = 0;
for k = 1:numel(energies)
    E = energies(k);  N0 = fluence(k);
    mu = exp(interp1(log(E_tab), log(mu_tab), log(E)));
    sig  = sig  + E*N0*exp(-mu*path_cm);
    sig0 = sig0 + E*N0;
end
scale = N0_PER_PIXEL / sig0;
sig   = poissrnd(sig*scale) / scale;          % incident statistics
sino  = -log(sig / sig0);                     % measured sinogram (1/pixel units)

% ---- (a) FBP reference ----------------------------------------------
rec_fbp = iradon(sino, theta, 'Hann', GRID) / PIXEL_SIZE_CM;

% ---- (b) SIRT, hand-implemented -------------------------------------
% normalisation: how much each pixel is "hit" by the projector pair
ones_sino = ones(size(sino));
norm_img  = iradon(ones_sino, theta, 'none', GRID);
norm_img(norm_img < 1e-6) = 1e-6;             % avoid divide-by-zero

norm_sino = radon(ones(GRID), theta);
norm_sino(norm_sino < 1e-6) = 1e-6;

rec = zeros(GRID);                             % start from a blank image
resid_hist = zeros(1, N_ITER);
for it = 1:N_ITER
    fwd   = radon(rec, theta);                 % what my guess would measure
    resid = (sino - fwd) ./ norm_sino;         % error, per-ray normalised
    corr  = iradon(resid, theta, 'none', GRID) ./ norm_img;   % back-project it
    rec   = rec + RELAX * corr;                % nudge the image
    rec(rec < 0) = 0;                          % physical: mu cannot be negative
    resid_hist(it) = norm(sino - radon(rec, theta), 'fro');
end
rec_sirt = rec / PIXEL_SIZE_CM;

% ---- numbers ---------------------------------------------------------
r = sqrt((xx-c).^2 + (yy-c).^2);
centre = r < 0.15*(CIRC_PIX/2);
edge   = (r > 0.80*(CIRC_PIX/2)) & (r < 0.92*(CIRC_PIX/2));
bg     = r > 1.15*(CIRC_PIX/2);                % outside the object
fprintf('%d angles, %g photons/pixel, %d SIRT iterations\n', ...
        N_ANGLES, N0_PER_PIXEL, N_ITER);
for p = {{'FBP', rec_fbp}, {'SIRT', rec_sirt}}
    nm = p{1}{1}; q = p{1}{2};
    fprintf('%-6s mu_centre=%.4f  mu_edge=%.4f  cupping=%5.2f%%  bg noise(sd)=%.4f\n', ...
            nm, mean(q(centre)), mean(q(edge)), ...
            100*(mean(q(edge))-mean(q(centre)))/mean(q(edge)), std(q(bg)));
end

% ---- figure ----------------------------------------------------------
figure;
subplot(2,2,1); imagesc(rec_fbp);  axis image off; colormap gray; colorbar;
title(sprintf('FBP (%d angles)', N_ANGLES));
subplot(2,2,2); imagesc(rec_sirt); axis image off; colorbar;
title(sprintf('SIRT (%d iterations)', N_ITER));
mid = round(GRID/2);
subplot(2,2,3);
plot(rec_fbp(mid,:), 'r'); hold on; plot(rec_sirt(mid,:), 'b'); grid on;
legend('FBP','SIRT'); xlabel('pixel'); ylabel('\mu (1/cm)');
title('Centre profiles');
subplot(2,2,4);
semilogy(resid_hist, 'k', 'LineWidth', 1.2); grid on;
xlabel('iteration'); ylabel('||sinogram residual||');
title('SIRT convergence');