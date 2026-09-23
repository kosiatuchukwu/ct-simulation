% fig_geometry.m - Figure 3.x: acquisition geometries
% Two panels: parallel-beam and fan-beam. Schematic only.

figure('Position',[100 100 950 420]);

% ================= PANEL (a): PARALLEL =================
subplot(1,2,1); hold on; axis equal off;

R = 1.0;                                   % object radius
th = linspace(0,2*pi,200);
fill(R*cos(th), R*sin(th), [0.85 0.90 0.95], 'EdgeColor',[0.3 0.3 0.4], 'LineWidth',1.2);

% parallel rays
for y = -1.5:0.5:1.5
    plot([-2.6 2.6],[y y], 'Color',[0.75 0.30 0.25], 'LineWidth',1.0);
end

% detector
plot([2.6 2.6],[-1.8 1.8],'k','LineWidth',3);
text(2.75,0,'detector','Rotation',90,'HorizontalAlignment','center','FontSize',9);

% detector coordinate t
plot([2.45 2.45],[0 1.0],'k','LineWidth',0.8);
plot(2.45,1.0,'k^','MarkerFaceColor','k','MarkerSize',4);
text(2.25,0.55,'$t$','Interpreter','latex','FontSize',11);

% rotation angle
plot(0,0,'k+','MarkerSize',8,'LineWidth',1.2);
ang = linspace(0,0.9,40);
plot(0.55*cos(ang),0.55*sin(ang),'k','LineWidth',0.8);
text(0.68,0.28,'$\theta$','Interpreter','latex','FontSize',11);
plot([0 1.6],[0 0],'k--','LineWidth',0.6);
plot([0 1.6*cos(0.9)],[0 1.6*sin(0.9)],'k--','LineWidth',0.6);

title('(a) Parallel-beam','FontSize',11);
xlim([-3.2 3.4]); ylim([-2.2 2.2]);

% ================= PANEL (b): FAN =================
subplot(1,2,2); hold on; axis equal off;

D  = 3.2;                                  % source to rotation centre
Dd = 1.6;                                  % centre to detector
fill(R*cos(th), R*sin(th), [0.85 0.90 0.95], 'EdgeColor',[0.3 0.3 0.4], 'LineWidth',1.2);

% source
plot(-D,0,'ko','MarkerFaceColor',[0.75 0.30 0.25],'MarkerSize',8);
text(-D,-0.42,'source','HorizontalAlignment','center','FontSize',9);

% diverging rays
det_half = 1.8;
for yd = linspace(-det_half, det_half, 9)
    plot([-D Dd],[0 yd],'Color',[0.75 0.30 0.25],'LineWidth',1.0);
end

% detector
plot([Dd Dd],[-det_half-0.15 det_half+0.15],'k','LineWidth',3);
text(Dd+0.22,0,'detector','Rotation',90,'HorizontalAlignment','center','FontSize',9);

% rotation centre
plot(0,0,'k+','MarkerSize',8,'LineWidth',1.2);

% D annotation
plot([-D 0],[-1.35 -1.35],'k','LineWidth',0.8);
plot([-D -D],[-1.28 -1.42],'k','LineWidth',0.8);
plot([0 0],[-1.28 -1.42],'k','LineWidth',0.8);
text(-D/2,-1.62,'$D$','Interpreter','latex','FontSize',11,'HorizontalAlignment','center');

% fan half-angle subtended by the object
a = atan2(R, D);
plot([-D -D+ (D)*1.0],[0 0],'k:','LineWidth',0.6);
ang2 = linspace(0,a,30);
plot(-D + 1.5*cos(ang2), 1.5*sin(ang2), 'k','LineWidth',0.8);
text(-D+1.62, 0.30, '$\gamma$','Interpreter','latex','FontSize',11);

% detector element spacing
plot([Dd+0.05 Dd+0.05],[0.45 0.90],'k','LineWidth',0.8);
text(Dd+0.14,0.68,'$\Delta_{\mathrm{det}}$','Interpreter','latex','FontSize',10);

title('(b) Fan-beam','FontSize',11);
xlim([-4.0 2.6]); ylim([-2.2 2.2]);

exportgraphics(gcf,'fig_geometry.png','Resolution',220);
disp('Saved fig_geometry.png');