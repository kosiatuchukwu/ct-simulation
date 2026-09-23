% ct_engine_demo.m - MATLAB port of the polychromatic CT engine
% Benchmark to reproduce (Python, 80 kVp): counting 6.07%, integrating 4.42%

KVP = 120;
GRID = 180; CIRC_PIX = 160; OBJ_DIAM_CM = 10.0;
N_ANGLES = 180;
PIXEL_SIZE_CM = OBJ_DIAM_CM / CIRC_PIX;

% ---- water mu(E) from the XCOM file --------------------------------
M = readmatrix('water_mu.txt', 'FileType','text');
M = M(~any(isnan(M),2), :);          % drop header/blank lines
E_tab  = M(:,1) * 1000;              % MeV -> keV
mu_tab = M(:,2) * 1.0;               % x density of water -> 1/cm

% ---- spectrum from SpekPy (exported CSV) ---------------------------
S = readmatrix(sprintf('spectrum_%dkvp.csv', KVP));
S = S(~any(isnan(S),2), :);
energies = S(:,1);  fluence = S(:,2);

% ---- phantom mask and geometric path lengths -----------------------
[xx, yy] = meshgrid(0:GRID-1, 0:GRID-1);
c = (GRID-1)/2;
mask = double(((xx-c).^2 + (yy-c).^2) <= (CIRC_PIX/2)^2);
theta = 0 : 180/N_ANGLES : 180-180/N_ANGLES;      % 180 angles, like Python
path_cm = radon(mask, theta) * PIXEL_SIZE_CM;

% ---- per-energy Beer-Lambert, both detector models -----------------
sig_c = zeros(size(path_cm));  sig_i = zeros(size(path_cm));
sig0_c = 0;  sig0_i = 0;
for k = 1:numel(energies)
    E  = energies(k);  N0 = fluence(k);
    mu = exp(interp1(log(E_tab), log(mu_tab), log(E)));
    T  = exp(-mu * path_cm);
    sig_c = sig_c + N0  * T;    sig0_c = sig0_c + N0;
    sig_i = sig_i + E*N0 * T;   sig0_i = sig0_i + E*N0;
end

% ---- reconstruction and cupping ------------------------------------
r = sqrt((xx-c).^2 + (yy-c).^2);
centre = r < 0.15*(CIRC_PIX/2);
edge   = (r > 0.80*(CIRC_PIX/2)) & (r < 0.92*(CIRC_PIX/2));

names = {'counting','integrating'};
sigs  = {sig_c, sig_i};  sig0s = {sig0_c, sig0_i};
recons = cell(1,2);
for m = 1:2
    sino = -log(sigs{m} / sig0s{m});
    rec  = iradon(sino, theta, 'Hamming', GRID) / PIXEL_SIZE_CM;
    recons{m} = rec;
    mu_c = mean(rec(centre));  mu_e = mean(rec(edge));
    fprintf('%-12s mu_centre=%.4f  mu_edge=%.4f  cupping=%.2f%%\n', ...
            names{m}, mu_c, mu_e, 100*(mu_e-mu_c)/mu_e);
end

% ---- figure --------------------------------------------------------
figure;
subplot(2,2,1); imagesc(recons{1}); axis image off; colormap gray; colorbar;
title(sprintf('Photon-counting (%d kVp)', KVP));
subplot(2,2,2); imagesc(recons{2}); axis image off; colorbar;
title(sprintf('Energy-integrating (%d kVp)', KVP));
mid = round(GRID/2);
subplot(2,1,2);
plot(recons{1}(mid,:), 'r'); hold on; plot(recons{2}(mid,:), 'b');
legend('photon-counting','energy-integrating'); grid on;
xlabel('pixel'); ylabel('\mu (1/cm)'); title('Centre profiles');