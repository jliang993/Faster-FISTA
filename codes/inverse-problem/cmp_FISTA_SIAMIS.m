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

outputType = 'pdf';
%% Forward--Backward
fprintf(sprintf('performing Forward-Backward...\n'));

[x0, its0, ek0, phik0] = func_FBS(para, gradF, proxJ, objPhi, J);

fprintf('\n');
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

p = 1/30;
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
%% Restarting FISTA
fprintf(sprintf('performing restarting FISTA...\n'));

r = 4;

p = 1;
q = 1;

[x_r, its_r, ek_r, phik_r] = func_FISTA_Restart(para, gradF, proxJ, objPhi, J, p,q,r);

fprintf('\n');
%% Adaptive + Restarting
fprintf(sprintf('performing RAda-FISTA...\n'));
r = 4;
p = 1.0/1.05;
q = p^2;

[x_ar, its_ar, ek_ar, phik_ar, r_ar, Rk, Vk,Wk] = func_RAdaFISTA(para, gradF, proxJ, objPhi, J, p,q,r);

fprintf('\n');

[its_a, its_ar, its_r]

fprintf('\n');
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

grey1 = [0.15,0.15,0.15];
p1 = semilogy(phik-phisol, 'Color',grey1, 'LineWidth',linewidth);


hold on,


% blue1 = [0.12,0.48,1.0];
% p1 = semilogy(phik_m1-phisol, 'Color',blue1, 'LineWidth',linewidth);
blue2 = [0.95,0.0,0.0];
p2 = semilogy(phik_m2-phisol, 'Color',blue2, 'LineWidth',linewidth);

pas1 = semilogy(phik_a-phisol, 'Color',[0.0,0.5,0.0], 'LineWidth',linewidth);
pasr = semilogy(phik_ar-phisol, 'Color',[0.99,0.01,0.99], 'LineWidth',linewidth);

pr = semilogy(phik_r-phisol, '-.', 'Color',[0.33,0.33,0.33], 'LineWidth',linewidth);


grid on;
ax = gca;
ax.GridLineStyle = '--';

% v = axis;
axis([1 length(phik)/1 1e-14 1e4]);
if strcmp(J, 'infty'); axis([1 length(phik)/10 1e-14 1e2]); end
ytick = [1e-12, 1e-8, 1e-4, 1e0, 1e4];
set(gca, 'yTick', ytick);

ylb = ylabel({'$\Phi(x_{k})-\Phi(x^\star)$'}, 'FontSize', labelFontSize,...
    'FontAngle', 'normal', 'Interpreter', 'latex');
set(ylb, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
xlb = xlabel({'\vspace{-1.0mm}';'$k$'}, 'FontSize', labelFontSize,...
    'FontAngle', 'normal', 'Interpreter', 'latex');
set(xlb, 'Units', 'Normalized', 'Position', [1/2, -0.055, 0]);

lg = legend([p1, p2, pas1,pasr, pr],...
    'FISTA-BT',...
    'FISTA-Mod, $p = \frac{1}{30}, q = \frac{1}{10}$',...
    'Ada-FISTA', '$\alpha$-RAda-FISTA',...
    'Restarting FISTA');
set(lg,'FontSize', 10);
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
