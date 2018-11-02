clear all
close all
clc

addpath toolbox
set(groot,'defaultLineLineWidth',1.5);
%% problem set up
J = 'lasso';
J = 'glasso';
% J = 'infty';
% J = 'tv';

[para, gradF,proxJ, objPhi] = problem_FB(J);

para.J = J;
para.tol = 1e-15;
para.maxits = 5e4 + 1;
para.gamma = para.beta;

para.x0 = zeros(para.n, 1);

para.verbose = 1;
%% Lazy-start FISTA-Mod     
fprintf(sprintf('performing Lazy-start FISTA...\n'));

p = 1/20;
q = 1/2;
r = 4;

[xsol, ~, ~, ~, ~] = func_FISTA_Mod(p,q,r, para, proxJ,gradF, objPhi, 0);

[x2, its2, dk2, ek2, fk2] = func_FISTA_Mod(p,q,r, para, proxJ,gradF, objPhi, xsol);

fprintf('\n');
%% FISTA, original 
fprintf(sprintf('performing original FISTA...\n'));

p = 1;
q = 1;
r = 4;

[x1, its1, dk1, ek1, fk1] = func_FISTA_Mod(p,q,r, para, proxJ,gradF, objPhi, xsol);

fprintf('\n');
%% Restarting FISTA     
fprintf(sprintf('performing restarting FISTA...\n'));

p = 1;
q = 1;
r = 4; 

[x3, its3, dk3, ek3, fk3] = func_Restart_FISTA(p,q,r, para, proxJ,gradF, objPhi, xsol);

fprintf('\n');
%% Rada-FISTA       
fprintf(sprintf('performing Rada-FISTA...\n'));

p = 1;
q = 1;
r = 4;

[x4, its4, dk4, ek4, fk4] = func_Rada_FISTA(p,q,r, para, proxJ,gradF, objPhi, xsol);

fprintf('\n');

disp([its1, its2, its3, its4])

fprintf('\n');
%% Greedy FISTA 
fprintf(sprintf('performing Greedy FISTA...\n'));

para.c_gamma = 1.3;
para.a = @(k) 1.0; %max(2/(1+k/5), 1.0);


[x5, its5, dk5, ek5, fk5] = func_Greedy_FISTA(para, proxJ,gradF, objPhi, xsol);

fprintf('\n');

disp([its1, its2, its3, its4, its5])
%% color map   
hh = parula;
%% distance error ||x_{k}-x^\star||  
linewidth = 1;

axesFontSize = 8;
labelFontSize = 8;
legendFontSize = 8;

resolution = 300; % output resolution
output_size = 300 *[10, 8]; % output size

%%%%%% relative error

figure(101), clf;
set(0,'DefaultAxesFontSize', axesFontSize);
set(gcf,'paperunits','centimeters','paperposition',[-0.1 -0.0 output_size/resolution]);
set(gcf,'papersize',output_size/resolution-[0.85 0.4]);

p1d = semilogy(dk1, 'color', hh(46,:), 'LineWidth',linewidth);
hold on,

p2d = semilogy(dk2, 'color', hh(31,:), 'LineWidth',linewidth);

p3d = semilogy(dk3, 'color', hh(16,:), 'LineWidth',linewidth);

p4d = semilogy(dk4, 'color', hh(1,:), 'LineWidth',linewidth);

p5d = semilogy(dk5, 'color', [1,0,0], 'LineWidth',linewidth);


uistack(p1d, 'bottom');

grid on;
ax = gca;
ax.GridLineStyle = '--';

axis([1, 4*its3, 1e-10, 2*max(dk1)]);
ytick = [1e-8, 1e-4, 1e-0, 1e4];
set(gca, 'yTick', ytick);

ylb = ylabel({'$\|x_{k}-x^\star\|$'}, 'FontSize', labelFontSize,...
    'FontAngle', 'normal', 'Interpreter', 'latex');
set(ylb, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
xlb = xlabel({'\vspace{-1.0mm}';'$k$'}, 'FontSize', labelFontSize,...
    'FontAngle', 'normal', 'Interpreter', 'latex');
set(xlb, 'Units', 'Normalized', 'Position', [1/2, -0.075, 0]);


lg = legend([p1d, p2d, p3d, p4d, p5d], ...
    'FISTA-BT',...
    'FISTA-Mod, $p = \frac{1}{20}, q = \frac{1}{2}$',...
    'Restarging FISTA', 'Rada-FISTA', 'Greedy FISTA');
set(lg,'FontSize', legendFontSize);
set(lg, 'Interpreter', 'latex');
% set(lg, 'Location', 'best');
legend('boxoff');


epsname = sprintf('cmp-ip-%s-dk.pdf', J);
print(epsname, '-dpdf');
epsname = sprintf('cmp-ip-%s-dk.png', J);
print(epsname, '-dpng');
%% plot Phi(x_{k}) - Phi(x*)    
min_f = min([min(fk1), min(fk1), min(fk2), min(fk3), min(fk4)]);
linewidth = 1;

axesFontSize = 8;
labelFontSize = 8;
legendFontSize = 8;

resolution = 300; % output resolution
output_size = 300 *[10, 8]; % output size

%%%%%% relative error

figure(103), clf;
set(0,'DefaultAxesFontSize', axesFontSize);
set(gcf,'paperunits','centimeters','paperposition',[-0.1 -0.0 output_size/resolution]);
set(gcf,'papersize',output_size/resolution-[0.85 0.4]);

p1f = semilogy(fk1-min_f, 'color', hh(46,:), 'LineWidth',linewidth);
hold on,

p2f = semilogy(fk2-min_f, 'color', hh(31,:), 'LineWidth',linewidth);

p3f = semilogy(fk3-min_f, 'color', hh(16,:), 'LineWidth',linewidth);

p4f = semilogy(fk4-min_f, 'color', hh(1,:), 'LineWidth',linewidth);

p5f = semilogy(fk5-min_f, 'color', [1,0,0], 'LineWidth',linewidth);


uistack(p1f, 'bottom');

grid on;
ax = gca;
ax.GridLineStyle = '--';

axis([1, 4*its3, 1e-14, 1*max(fk1)]);
ytick = [1e-14, 1e-10, 1e-6, 1e-2, 1e2, 1e6];
set(gca, 'yTick', ytick);

ylb = ylabel({'$\Phi(x_{k})-\Phi(x^\star)$'}, 'FontSize', labelFontSize,...
    'FontAngle', 'normal', 'Interpreter', 'latex');
set(ylb, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
xlb = xlabel({'\vspace{-1.0mm}';'$k$'}, 'FontSize', labelFontSize,...
    'FontAngle', 'normal', 'Interpreter', 'latex');
set(xlb, 'Units', 'Normalized', 'Position', [1/2, -0.075, 0]);


lg = legend([p1f, p2f, p3f, p4f, p5f], ...
    'FISTA-BT',...
    'FISTA-Mod, $p = \frac{1}{20}, q = \frac{1}{2}$',...
    'Restarging FISTA', 'Rada-FISTA', 'Greedy FISTA');
set(lg,'FontSize', legendFontSize);
% set(lg, 'Location', 'best');
set(lg, 'Interpreter', 'latex');
legend('boxoff');

epsname = sprintf('cmp-ip-%s-fk.pdf', J);
print(epsname, '-dpdf');
epsname = sprintf('cmp-ip-%s-fk.png', J);
print(epsname, '-dpng');