% fig_phantom_masks.m - Figure 3.1: phantom construction
GRID = 180; CIRC_PIX = 160;
OBJ_DIAM_CM = 10.0; DET_DIAM_CM = 2.0;

[xx, yy] = meshgrid(0:GRID-1, 0:GRID-1);
c = (GRID-1)/2;

obj_mask = ((xx-c).^2 + (yy-c).^2) <= (CIRC_PIX/2)^2;

det_r = (DET_DIAM_CM / OBJ_DIAM_CM) * CIRC_PIX / 2;
dc    = [c + 0.7*CIRC_PIX/2, c];
det_mask = (((xx-dc(1)).^2 + (yy-dc(2)).^2) <= det_r^2) & obj_mask;

obj_only = obj_mask & ~det_mask;

figure('Position',[100 100 900 320]);
subplot(1,3,1); imagesc(obj_mask);  axis image off; colormap(gray);
title('(a) Object mask');
subplot(1,3,2); imagesc(det_mask);  axis image off;
title('(b) Detail mask');
subplot(1,3,3); imagesc(obj_only);  axis image off;
title('(c) Object excluding detail');

exportgraphics(gcf, 'fig_phantom_masks.png', 'Resolution', 200);
disp('Saved fig_phantom_masks.png');