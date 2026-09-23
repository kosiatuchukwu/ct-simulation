% fig_kvp_profiles.m - Figure 4.1: centre-line profiles at 50, 80, 120 kVp
% Homogeneous water phantom, noiseless, energy-integrating, Hann filter.

kvps = [50 80 120];
cols = {'#c0392b', '#2c3e50', '#27ae60'};
figure('Position',[100 100 700 450]); hold on;

for i = 1:numel(kvps)
    p = struct('kvp', kvps(i), 'detector','integrating', ...
               'det_material','', 'recon_filter','Hann');
    r = ct_simulate(p);
    mid = round(size(r.recon,1)/2);
    plot(r.recon(mid,:), 'Color', cols{i}, 'LineWidth', 1.4, ...
         'DisplayName', sprintf('%d kVp (cupping %.2f%%)', kvps(i), r.cupping));
    fprintf('%3d kVp: mu_centre=%.4f  mu_edge=%.4f  cupping=%.2f%%\n', ...
            kvps(i), r.mu_centre, r.mu_edge, r.cupping);
end

grid on; box on;
xlabel('Pixel position across the phantom');
ylabel('\mu (cm^{-1})');
title('Centre-line attenuation profiles');
legend('Location','south');

exportgraphics(gcf, 'fig_kvp_profiles.png', 'Resolution', 200);
disp('Saved fig_kvp_profiles.png');