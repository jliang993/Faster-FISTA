clear all
% close all
clc
%% new examples
% Line Y
A = [-0.2, 1];
a = 1/2;

ProjY = @(x) x - (A')/ (A*A') *(A*x - a);

% Line X
B = [-0.15, 1;];
b = 1/2;

ProjX = @(x) x - (B')/ (B*B') *(B*x - b);

para.maxits = 1e4;
para.tol = 1e-12;

outputType = 'pdf';
%% plot the two sets
% t = 5e-2;
% 
% Xw = -2.5;
% Xe = 2.5;
% 
% Ys = -1.5;
% Yn = 2.5;
%
% 
% resolution = 300; % output resolution
% output_size = 300 *[8, 8]; % output size
% 
% figure(100), clf;
% set(0,'DefaultAxesFontSize', axesFontSize);
% set(gcf,'paperunits','centimeters','paperposition',[-0.8 -1.5 output_size/resolution]);
% set(gcf,'papersize',output_size/resolution-[1.3 2.7]);
% 
% % % plot the axis
% % plot([Xw,Xe], [0, 0], 'k--', 'LineWidth', 1);
% % hold on;
% % plot([0, 0], [Ys,Yn], 'k--', 'LineWidth', 1);
% 
% % Set A
% plot([Xw, Xe], -A(1)*[Xw, Xe]+a, 'k', 'LineWidth', 1.25);
% hold on;
% 
% plot([Xw, Xe], -B(1)*[Xw, Xe]+b, 'r', 'LineWidth', 1.25);
% 
% % axis properties
% grid on;
% axis equal;
% axis([Xw Xe Ys Yn]);
% 
% set(gca, 'yTick', [-2:1:2]);
% set(gca, 'xTick', [-2:1:2]);
% %
% set(gca, 'yTickLabel', {' ', ' ', ' ', ' ', ' '});
% set(gca, 'xTickLabel', {' ', ' ', ' ', ' ', ' '});
% 
% grid on;
% ax = gca;
% ax.GridLineStyle = '--';
% 
% print('dist-of-2lines.pdf', '-dpdf');
%%
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

p = 1/30;
q = 1/1e1;

[x_m2,y_m2, its_m2, ek_m2] = func_FISTA_Mod(p,q,r, para, ProjX, ProjY);

fprintf('\n');
%% Adaptive + Restarting
fprintf(sprintf('performing restarting AdaFISTA...\n'));

r = 4;
p = 1/2;
q = p^2;

[x_ra,y_ra, its_ra, ek_ra] = func_RAdaFISTA(p,q,r, para, ProjX, ProjY);

fprintf('\n');
%% Restarting FISTA
fprintf(sprintf('performing restarting FISTA...\n'));

r = 4;

p = 1;
q = 1;

[x_r,y_r, its_r, ek_r] = func_FISTA_Restart(p,q,r, para, ProjX, ProjY);

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

grey1 = [0.15,0.15,0.15];
p1 = semilogy(ek, 'Color',grey1, 'LineWidth',linewidth);


hold on,


% blue1 = [0.12,0.48,1.0];
% p1 = semilogy(phik_m1-phisol, 'Color',blue1, 'LineWidth',linewidth);
blue2 = [0.9,0.0,0.0];
p2 = semilogy(ek_m2, 'Color',blue2, 'LineWidth',linewidth);

pasr = semilogy(ek_ra, 'Color',[0.99,0.01,0.99], 'LineWidth',linewidth);

pr = semilogy(ek_r, '-.', 'Color',[0.33,0.33,0.33], 'LineWidth',linewidth);

grid on;
ax = gca;
ax.GridLineStyle = '--';

% v = axis;
axis([1 length(ek)/3 1e-10 1e2]);
ytick = [1e-10, 1e-6, 1e-2, 1e2];
set(gca, 'yTick', ytick);

ylb = ylabel({'$\|x_{k}-x_{k-1}\|$'}, 'FontSize', labelFontSize,...
    'FontAngle', 'normal', 'Interpreter', 'latex');
set(ylb, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
xlb = xlabel({'\vspace{-1.0mm}';'$k$'}, 'FontSize', labelFontSize,...
    'FontAngle', 'normal', 'Interpreter', 'latex');
set(xlb, 'Units', 'Normalized', 'Position', [1/2, -0.055, 0]);

lg = legend([p1, p2, pasr, pr], 'FISTA-BT',...
    'FISTA-Mod, $p = \frac{1}{30}, q = \frac{1}{10}$',...
    '$\alpha$-RAda-FISTA',...
    'Restarting FISTA');
set(lg,'FontSize', legendFontSize);
set(lg, 'Interpreter', 'latex');
%
pos = get(lg, 'Position');
set(lg, 'Position', [pos(1)-0.125, pos(2)-0.075, pos(3:4)]);
pos_ = get(lg, 'Position');
legend('boxoff');


epsname = sprintf('cmp_fista_ek_2lines.%s', outputType);
if strcmp(outputType, 'png')
    print(epsname, '-dpng');
else
    print(epsname, '-dpdf');
end