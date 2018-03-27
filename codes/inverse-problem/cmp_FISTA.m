clear all
close all
clc
%% problem set up
J = 'lasso';
% J = 'glasso';
% J = 'infty';

[para, gradF,proxJ, objPhi] = problem_setup(J);

para.tol = 1e-11;
para.maxits = 1e5;

outputType = 'png';
%% Original FISTA-BT
fprintf(sprintf('performing FISTA-BT...\n'));

r = 4;

p = 1;
q = 1;
[x, its, ek, phik] = func_FISTA_Mod(para, gradF, proxJ, objPhi, J, p,q,r);

fprintf('\n');
%% FISTA-Mod
fprintf(sprintf('performing FISTA-Mod...\n'));

p = 1/5;
q = 1/1;

[x_m1, its_m1, ek_m1, phik_m1] = func_FISTA_Mod(para, gradF, proxJ, objPhi, J, p,q,r);
% fprintf('\n');

%

p = 1/50;
q = 1/1e1;

[x_m2, its_m2, ek_m2, phik_m2] = func_FISTA_Mod(para, gradF, proxJ, objPhi, J, p,q,r);

fprintf('\n');
%% Adaptive-FISTA, AdaFISTA
fprintf(sprintf('performing Ada-FISTA...\n'));
r = 4;
p = 1/1;
q = p^2;

[x_a, its_a, ek_a, phik_a, r_a, Rk_a, Vk_a] = func_AdaFISTA_s1(para, gradF, proxJ, objPhi, J, p,q,r);

fprintf('\n');
%% Adaptive + Restarting
r = 4;
p = 1/2;
q = p^2;

[x_ar, its_ar, ek_ar, phik_ar, r_ar, Rk, Vk,Wk] = func_AdaFISTA_sR(para, gradF, proxJ, objPhi, J, p,q,r);

fprintf('\n');
%% FISTA-AC
fprintf(sprintf('performing FISTA-CD...\n'));

d = 2;
[x_c1, its_c1, ek_c1, phik_c1] = func_FISTA_CD(para, gradF, proxJ, objPhi, J, d);

%
d = 75;
[x_c2, its_c2, ek_c2, phik_c2] = func_FISTA_CD(para, gradF, proxJ, objPhi, J, d);
fprintf('\n\n');
%% Restarting FISTA
fprintf(sprintf('performing restarting FISTA...\n'));

r = 4;

p = 1;
q = 1;

[x_r, its_r, ek_r, phik_r] = func_FISTA_Restart(para, gradF, proxJ, objPhi, J, p,q,r);

fprintf('\n');

[its_a, its_ar, its_r]
%% plot ||x_{k} - x_{k-1}||
linewidth = 1;

axesFontSize = 8;
labelFontSize = 8;
legendFontSize = 8;

resolution = 300; % output resolution
output_size = 300 *[10, 8]; % output size

%%%%%% relative error

figure(100), clf;
set(0,'DefaultAxesFontSize', axesFontSize);
set(gcf,'paperunits','centimeters','paperposition',[-0.1 -0.1 output_size/resolution]);
set(gcf,'papersize',output_size/resolution-[0.85 0.5]);

grey1 = [0.3,0.3,0.3];
p1 = semilogy(ek, 'Color',grey1, 'LineWidth',linewidth);
hold on,

% blue1 = [0.12,0.48,1.0];
% p1 = semilogy(ek_m1, 'Color',blue1, 'LineWidth',linewidth);
blue2 = [0.0,0.0,0.85];
p2 = semilogy(ek_m2, 'Color',blue2, 'LineWidth',linewidth);

pas1 = semilogy(ek_a, 'Color',[0.0,0.5,0.0], 'LineWidth',linewidth);
pasr = semilogy(ek_ar, 'Color',[0.99,0.01,0.99], 'LineWidth',linewidth);

pr = semilogy(ek_r, '-.', 'Color',[0.33,0.33,0.33], 'LineWidth',linewidth);

%%%%%
red1 = [1.0,0.4,0.4];
p3 = semilogy(ek_c1, 'Color',red1, 'LineWidth',linewidth);
red2 = [0.85,0.0,0.0];
p4 = semilogy(ek_c2, 'Color',red2, 'LineWidth',linewidth);

grid on;
ax = gca;
ax.GridLineStyle = '--';

% v = axis;
axis([1 length(ek)/1 1e-10 1e2]);
if strcmp(J, 'infty'); axis([1 length(ek)/6 1e-10 1e2]); end
ytick = [1e-10, 1e-6, 1e-2, 1e2];
set(gca, 'yTick', ytick);

ylb = ylabel({'$\|x_{k}-x_{k-1}\|$'}, 'FontSize', labelFontSize,...
    'FontAngle', 'normal', 'Interpreter', 'latex');
set(ylb, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
xlb = xlabel({'\vspace{-1.0mm}';'$k$'}, 'FontSize', labelFontSize,...
    'FontAngle', 'normal', 'Interpreter', 'latex');
set(xlb, 'Units', 'Normalized', 'Position', [1/2, -0.055, 0]);

lg = legend([p1, p2, p3,p4, pas1,pasr, pr], 'FISTA-BT',...
    'FISTA-Mod, $p = \frac{1}{50}, q = \frac{1}{10}$',...
    'FISTA-CD, $d = 2$', 'FISTA-CD, $d = 75$',...
    'Ada-FISTA, $p=1, q=1$', 'Ada-FISTA+Restart, $p=\frac{1}{2}, q=1$',...
    'Restarting FISTA');
set(lg,'FontSize', legendFontSize);
set(lg, 'Interpreter', 'latex');
%
pos = get(lg, 'Position');
set(lg, 'Position', [pos(1)-0.125, pos(2)-0.075, pos(3:4)]);
pos_ = get(lg, 'Position');
legend('boxoff');


epsname = sprintf('cmp_fista_ek_%s.%s', J, outputType);
if strcmp(outputType, 'png')
    print(epsname, '-dpng');
else
    print(epsname, '-dpdf');
end
%% plot Phi(x_{k}) - Phi(x*) 
phisol = min(phik);

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

grey1 = [0.3,0.3,0.3];
p1 = semilogy(phik-phisol, 'Color',grey1, 'LineWidth',linewidth);
hold on,

% blue1 = [0.12,0.48,1.0];
% p1 = semilogy(phik_m1-phisol, 'Color',blue1, 'LineWidth',linewidth);
blue2 = [0.0,0.0,0.85];
p2 = semilogy(phik_m2-phisol, 'Color',blue2, 'LineWidth',linewidth);

pas1 = semilogy(phik_a-phisol, 'Color',[0.0,0.5,0.0], 'LineWidth',linewidth);
pasr = semilogy(phik_ar-phisol, 'Color',[0.99,0.01,0.99], 'LineWidth',linewidth);

pr = semilogy(phik_r-phisol, '-.', 'Color',[0.33,0.33,0.33], 'LineWidth',linewidth);

%%%%%
red1 = [1.0,0.4,0.4];
p3 = semilogy(phik_c1-phisol, 'Color',red1, 'LineWidth',linewidth);
red2 = [0.85,0.0,0.0];
p4 = semilogy(phik_c2-phisol, 'Color',red2, 'LineWidth',linewidth);

grid on;
ax = gca;
ax.GridLineStyle = '--';

% v = axis;
axis([1 length(phik)/1 1e-14 1e4]);
if strcmp(J, 'infty'); axis([1 length(phik)/6 1e-14 1e2]); end
ytick = [1e-12, 1e-8, 1e-4, 1e0, 1e4];
set(gca, 'yTick', ytick);

ylb = ylabel({'$\Phi(x_{k})-\Phi(x^\star)$'}, 'FontSize', labelFontSize,...
    'FontAngle', 'normal', 'Interpreter', 'latex');
set(ylb, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
xlb = xlabel({'\vspace{-1.0mm}';'$k$'}, 'FontSize', labelFontSize,...
    'FontAngle', 'normal', 'Interpreter', 'latex');
set(xlb, 'Units', 'Normalized', 'Position', [1/2, -0.055, 0]);

lg = legend([p1, p2, p3,p4, pas1,pasr, pr], 'FISTA-BT',...
    'FISTA-Mod, $p = \frac{1}{50}, q = \frac{1}{10}$',...
    'FISTA-CD, $d = 2$', 'FISTA-CD, $d = 75$',...
    'Ada-FISTA, $p=1, q=1$', 'Ada-FISTA+Restart, $p=\frac{1}{2}, q=1$',...
    'Restarting FISTA');
set(lg,'FontSize', legendFontSize);
set(lg, 'Interpreter', 'latex');
%
pos = get(lg, 'Position');
set(lg, 'Position', [pos(1)-0.125, pos(2)-0.075, pos(3:4)]);
pos_ = get(lg, 'Position');
legend('boxoff');


epsname = sprintf('cmp_fista_phik_%s.%s', J, outputType);
if strcmp(outputType, 'png')
    print(epsname, '-dpng');
else
    print(epsname, '-dpdf');
end