clear all
% close all
clc
%% two sets
% Polyhedral set B
B1 = [0.2, 1;];
B2 = [-0.2, 1;];
b2 = 1;
gap = b2 - B2(1);
b1 = B1(1) + gap;

C1 = [-B1(2)/B1(1), 1];
c1 = - B1(1) + C1(1) + b1;
C2 = [-B2(2)/B2(1), 1];
c2 = - B2(1) + C2(1) + b2;

ProjX = @(x) Proj2X(x, B1,B2,b1,b2);

% Affine set A
A = [-0.19, 1];
a = 1/2;
% a = A(1) + gap;

ProjY = @(x) x - (A')/ (A*A') *(A*x - a);

para.maxits = 1e4;
para.tol = 1e-12;

outputType = 'pdf';
%% plot the two sets
t = 5e-2;

Xw = -9.0;
Xe = 11.0;

Ys = -2.0;
Yn = 14.0;

axesFontSize = 8;
labelFontSize = 10;
legendFontSize = 8;

resolution = 300; % output resolution
output_size = 300 *[10, 8]; % output size

figure(100), clf;
set(0,'DefaultAxesFontSize', axesFontSize);
set(gcf,'paperunits','centimeters','paperposition',[-0.5 -0.5 output_size/resolution]);
set(gcf,'papersize',output_size/resolution-[1.1 1.1]);

% plot the axis
plot([Xw,Xe], [0, 0], 'k--', 'LineWidth', 1);
hold on;
plot([0, 0], [Ys,Yn], 'k--', 'LineWidth', 1);

% Set A
plot([Xw, Xe], -A(1)*[Xw, Xe]+a, 'k', 'LineWidth', 1.0);

plot([Xw, 1], -B1(1)*[Xw, 1]+b1, 'r', 'LineWidth', 1.0);
plot([1, Xe], -B2(1)*[1, Xe]+b2, 'r', 'LineWidth', 1.0);

% plot([Xw, Xe], -C1(1)*[Xw, Xe]+c1, 'b', 'LineWidth', 1.5);
% plot([Xw, Xe], -C2(1)*[Xw, Xe]+c2, 'b', 'LineWidth', 1.5);

% % Set B
patch([Xw, 1, Xe, Xe, Xw], [-B1(1)*Xw+b1, -B1(1)*1+b1, -B2(1)*Xe+b2, Yn,Yn],...
    'r','FaceAlpha',0.3, 'EdgeAlpha',0.0, 'EdgeColor','r');

text(-8,-1, '$\mathcal{X}$', 'FontSize',10, 'Interpreter', 'latex')
text(-8,3, '$\mathcal{Y}$', 'FontSize',10, 'Interpreter', 'latex')

 
% axis properties
grid on;
axis equal;
axis([Xw Xe Ys Yn]);

set(gca, 'yTick', [Ys:2:Yn]);
set(gca, 'xTick', [Xw:2:Xe]);

% print('dist-of-2sets.pdf', '-dpdf');
%% Original FISTA-BT
fprintf(sprintf('performing FISTA-BT...\n'));

r = 4;

p = 1;
q = 1;
[x,y, its, ek] = func_FISTA_Mod(p,q,r, para, ProjX, ProjY);

fprintf('\n');
%% FISTA-Mod
fprintf(sprintf('performing FISTA-Mod...\n'));

p = 1/5;
q = 1/1;

[x_m1,y_m1, its_m1, ek_m1] = func_FISTA_Mod(p,q,r, para, ProjX, ProjY);
% fprintf('\n');

%

p = 1/50;
q = 1/1e1;

[x_m2,y_m2, its_m2, ek_m2] = func_FISTA_Mod(p,q,r, para, ProjX, ProjY);

fprintf('\n');
%% plot ||x_{k}-x_{k-1}||
linewidth = 1;

axesFontSize = 8;
labelFontSize = 8;
legendFontSize = 8;

resolution = 300; % output resolution
output_size = 300 *[10, 8]; % output size

%%%%%% relative error

figure(101), clf;
set(0,'DefaultAxesFontSize', axesFontSize);
set(gcf,'paperunits','centimeters','paperposition',[-0.1 -0.1 output_size/resolution]);
set(gcf,'papersize',output_size/resolution-[0.85 0.5]);

grey1 = [0.65,0.65,0.65];
p1 = semilogy(ek, 'Color',grey1, 'LineWidth',linewidth);


hold on,


% blue1 = [0.12,0.48,1.0];
% p1 = semilogy(phik_m1-phisol, 'Color',blue1, 'LineWidth',linewidth);
blue2 = [0.0,0.0,0.85];
p2 = semilogy(ek_m2, 'Color',blue2, 'LineWidth',linewidth);

grid on;
ax = gca;
ax.GridLineStyle = '--';

% v = axis;
axis([1 length(ek)/1 1e-14 1e4]);
ytick = [1e-12, 1e-8, 1e-4, 1e0, 1e4];
set(gca, 'yTick', ytick);

ylb = ylabel({'$\|x_{k}-x_{k-1}\|$'}, 'FontSize', labelFontSize,...
    'FontAngle', 'normal', 'Interpreter', 'latex');
set(ylb, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
xlb = xlabel({'\vspace{-1.0mm}';'$k$'}, 'FontSize', labelFontSize,...
    'FontAngle', 'normal', 'Interpreter', 'latex');
set(xlb, 'Units', 'Normalized', 'Position', [1/2, -0.055, 0]);

lg = legend([p1, p2],...
    'FISTA-BT',...
    'FISTA-Mod, $p = \frac{1}{50}, q = \frac{1}{10}$');
set(lg,'FontSize', legendFontSize);
set(lg, 'Interpreter', 'latex');
%
pos = get(lg, 'Position');
set(lg, 'Position', [pos(1)-0.125, pos(2)-0.075, pos(3:4)]);
pos_ = get(lg, 'Position');
legend('boxoff');


epsname = sprintf('cmp_fista_ek_bp.%s', outputType);
if strcmp(outputType, 'png')
    print(epsname, '-dpng');
else
    print(epsname, '-dpdf');
end