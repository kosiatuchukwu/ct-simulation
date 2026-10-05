% ct_fanbeam_demo.m - fan-beam geometry, fixed parameters
% Fan forward projection + fan reconstruction, with a deliberately-wrong
% parallel reconstruction of the same data for comparison.

KVP            = 80;
SOURCE_DIST_PIX = 2000;    % source-to-rotation-centre distance (pixels)
DET_PIXEL_PIX   = 1.0;    % detector pixel size (pixels) = FanSensorSpacing
ROT_INCREMENT   = 1.0;    % degrees between views (360 views over a full turn)
RECON_FILTER    = 'Hann';

GRID = 180; CIRC_PIX = 160; OBJ_DIAM_CM = 10.0;
PIXEL_SIZE_CM = OBJ_DIAM_CM / CIRC_PIX;

% ---- water mu(E) ----------------------------------------------------
M = readmatrix('water_mu.txt', 'FileType','text');
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

% ---- FAN geometry: path lengths through the object ------------------
path_fan = fanbeam(mask, SOURCE_DIST_PIX, ...
                   'FanSensorGeometry','line', ...
                   'FanSensorSpacing', DET_PIXEL_PIX, ...
                   'FanRotationIncrement', ROT_INCREMENT) * PIXEL_SIZE_CM;

% ---- polychromatic, energy-integrating detector ---------------------
sig = zeros(size(path_fan));  sig0 = 0;
for k = 1:numel(energies)
    E = energies(k);  N0 = fluence(k);
    mu = exp(interp1(log(E_tab), log(mu_tab), log(E)));
    sig  = sig  + E*N0 * exp(-mu*path_fan);
    sig0 = sig0 + E*N0;
end
sino_fan = -log(sig / sig0);

% ---- reconstruct: (a) correctly, with ifanbeam ----------------------
rec_fan = ifanbeam(sino_fan, SOURCE_DIST_PIX, ...
                   'FanSensorGeometry','line', ...
                   'FanSensorSpacing', DET_PIXEL_PIX, ...
                   'FanRotationIncrement', ROT_INCREMENT, ...
                   'OutputSize', GRID, ...
                   'Filter', RECON_FILTER) / PIXEL_SIZE_CM;

% ---- reconstruct: (b) deliberately WRONG, as if parallel ------------
% treat each fan view as a parallel projection: valid only if the source
% is far away compared with the object size.
theta_par = 0 : ROT_INCREMENT : 360-ROT_INCREMENT;
rec_par = iradon(sino_fan, theta_par, RECON_FILTER, GRID) / PIXEL_SIZE_CM;

% ---- numbers --------------------------------------------------------
r = sqrt((xx-c).^2 + (yy-c).^2);
centre = r < 0.15*(CIRC_PIX/2);
edge   = (r > 0.80*(CIRC_PIX/2)) & (r < 0.92*(CIRC_PIX/2));
fprintf('Source distance = %g px  (object diameter = %d px, ratio %.1f)\n', ...
        SOURCE_DIST_PIX, CIRC_PIX, SOURCE_DIST_PIX/CIRC_PIX);
for p = {{'fan recon', rec_fan}, {'parallel recon of fan data', rec_par}}
    nm = p{1}{1}; rec = p{1}{2};
    mu_c = mean(rec(centre));  mu_e = mean(rec(edge));
    fprintf('%-28s mu_centre=%.4f  mu_edge=%.4f  cupping=%.2f%%\n', ...
            nm, mu_c, mu_e, 100*(mu_e-mu_c)/mu_e);
end

% ---- figure ---------------------------------------------------------
figure;
subplot(2,2,1); imagesc(sino_fan); colormap gray; colorbar;
title('Fan-beam sinogram'); xlabel('view'); ylabel('detector element');
subplot(2,2,2); imagesc(rec_fan); axis image off; colorbar;
title(sprintf('Fan reconstruction (D = %g px)', SOURCE_DIST_PIX));
subplot(2,2,3); imagesc(rec_par); axis image off; colorbar;
title('Parallel reconstruction of the SAME fan data');
mid = round(GRID/2);
subplot(2,2,4);
plot(rec_fan(mid,:), 'b'); hold on; plot(rec_par(mid,:), 'r'); grid on;
legend('fan recon','parallel recon'); xlabel('pixel'); ylabel('\mu (1/cm)');
title('Centre profiles');