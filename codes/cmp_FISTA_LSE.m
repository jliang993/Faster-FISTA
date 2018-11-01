clear all
close all
clc

addpath toolbox
set(groot,'defaultLineLineWidth',1.5);
%% set-up        
n = 2e2 + 1;

A = 2*eye(n) - diag(ones(n-1,1), -1) - diag(ones(n-1,1), 1);

x_ob = randn(n, 1);
b_ob = A *x_ob;
b = b_ob + 1e-1*randn(n, 1);
b = 0*b;

% para.beta = 1/norm(A)^2;
para.gamma = 1.0 /norm(A)^2;
para.maxits = 1e6 + 1;
para.tol = 1e-16;
para.n = n;

para.mu = 0;

proxJ = @(x, t) x;
gradF = @(x) (A')*(A*x - b);
objF = @(x) norm(A*x-b)^2 /2;

para.x0 = 1e4*ones(n, 1);

para.verbose = 1;
%% compute strong convexity     
gamma = para.gamma;

v = svd((A')*A);
alpha = min(v);

eta = 1 - gamma*alpha;
a_opt = (1-sqrt(1-eta))/(1+sqrt(1-eta));
%% Optimal scheme 
fprintf(sprintf('performing optimal scheme...\n'));

p = 1;
q = 1;
r = 4*(1-p) ...
    + 4*p*(1-sqrt(gamma*alpha))/(1+sqrt(gamma*alpha)) ...
    + 4*gamma*alpha*(p^2-q)/((1+sqrt(gamma*alpha))^2);

[xsol, ~, ~, ~, ~] = func_FISTA_Mod(p,q,r, para, proxJ,gradF, objF, 0);
xsol = 0*xsol;

[x2, its2, dk2, ek2, fk2] = func_FISTA_Mod(p,q,r, para, proxJ,gradF, objF, xsol);

fprintf('\n');
%% Gradient Descent 
fprintf(sprintf('performing Gradient Descent...\n'));

p = 1;
q = 1;
r = 0;

[x1, its1, dk1, ek1, fk1] = func_FISTA_Mod(p,q,r, para, proxJ,gradF, objF, xsol);

fprintf('\n');
%% FISTA, original 
fprintf(sprintf('performing original FISTA...\n'));

p = 1;
q = 1;
r = 4;

[x3, its3, dk3, ek3, fk3] = func_FISTA_Mod(p,q,r, para, proxJ,gradF, objF, xsol);

fprintf('\n');
%% Lazy-start FISTA-Mod 
fprintf(sprintf('performing Lazy-FISTA...\n'));

p = 1/20;
q = 1/2;
r = 4;

[x4, its4, dk4, ek4, fk4] = func_FISTA_Mod(p,q,r, para, proxJ,gradF, objF, xsol);

fprintf('\n');
%% Restarting FISTA     
fprintf(sprintf('performing restarting FISTA...\n'));

p = 1;
q = 1;
r = 4; 

[x5, its5, dk5, ek5, fk5] = func_Restart_FISTA(p,q,r, para, proxJ,gradF, objF, xsol);

fprintf('\n');
%% Rada-FISTA       
fprintf(sprintf('performing Rada-FISTA...\n'));

p = 1;
q = 1;
r = 4;

[x6, its6, dk6, ek6, fk6] = func_Rada_FISTA(p,q,r, para, proxJ,gradF, objF, xsol);

fprintf('\n');
%% Greedy FISTA 
fprintf(sprintf('performing Greedy FISTA...\n'));

para.c_gamma = 1.3;
para.a = @(k) 1.0; %max(2/(1+k/5), 1.0);


[x7, its7, dk7, ek7, fk7] = func_Greedy_FISTA(para, proxJ,gradF, objF, xsol);

fprintf('\n');

disp([its3, its4, its5, its6, its7]);
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

p1d = semilogy(dk1, 'color', hh(51,:), 'LineWidth',linewidth);
hold on,

p3d = semilogy(dk3, 'color', hh(41,:), 'LineWidth',linewidth);

p4d = semilogy(dk4, 'color', hh(31,:), 'LineWidth',linewidth);

p2d = semilogy(dk2, 'color', hh(21,:), 'LineWidth',linewidth);
p5d = semilogy(dk5, 'color', hh(11,:), 'LineWidth',linewidth);

p6d = semilogy(dk6, 'color', .8*hh(1,:), 'LineWidth',linewidth);

p7d = semilogy(dk7, 'color', [1,0,0], 'LineWidth',linewidth);

uistack(p3d, 'bottom');

delete(p1d);

grid on;
ax = gca;
ax.GridLineStyle = '--';

axis([1, min([3*find(dk2<1e-10, 1), para.maxits]), 1e-6, 2*max(dk1)]);
ytick = [1e-10, 1e-6, 1e-2, 1e2];
set(gca, 'yTick', ytick);

ylb = ylabel({'$\|x_{k}-x^\star\|$'}, 'FontSize', labelFontSize,...
    'FontAngle', 'normal', 'Interpreter', 'latex');
set(ylb, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
xlb = xlabel({'\vspace{-1.0mm}';'$k$'}, 'FontSize', labelFontSize,...
    'FontAngle', 'normal', 'Interpreter', 'latex');
set(xlb, 'Units', 'Normalized', 'Position', [1/2, -0.075, 0]);


% lg = legend([p1d, p3d, p4d, p2d, p5d, p6d, p7d], ...
%     'Gradient Descent', 'FISTA-BT',...
%     'FISTA-Mod, $p = \frac{1}{20}, q = \frac{1}{2}$', 'Optimal scheme',...
%     'Restarging FISTA', 'Rada-FISTA', 'Greedy FISTA');
lg = legend([p3d, p4d, p2d, p5d, p6d, p7d], ...
    'FISTA-BT',...
    'FISTA-Mod, $p = \frac{1}{20}, q = \frac{1}{2}$', '$\alpha$-FISTA',...
    'Restarging FISTA', 'Rada-FISTA', 'Greedy FISTA');
set(lg,'FontSize', legendFontSize);
set(lg, 'Interpreter', 'latex');
% set(lg, 'Location', 'southeast');
legend('boxoff');

filename = ['results', filesep, sprintf('cmp-lse-dk.pdf')];
print(filename, '-dpdf');
filename = ['results', filesep, sprintf('cmp-lse-dk.png')];
print(filename, '-dpng');
%% plot Phi(x_{k}) - Phi(x*)    
min_f = min([min(fk1), min(fk2), min(fk3), min(fk4), min(fk5), min(fk6)]);
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

p1f = semilogy(fk1-min_f, 'color', hh(51,:), 'LineWidth',linewidth);
hold on,

p3f = semilogy(fk3-min_f, 'color', hh(41,:), 'LineWidth',linewidth);

p4f = semilogy(fk4-min_f, 'color', hh(31,:), 'LineWidth',linewidth);

p2f = semilogy(fk2-min_f, 'color', hh(21,:), 'LineWidth',linewidth);
p5f = semilogy(fk5-min_f, 'color', hh(11,:), 'LineWidth',linewidth);

p6f = semilogy(fk6-min_f, 'color', .8*hh(1,:), 'LineWidth',linewidth);

p7f = semilogy(fk7-min_f, 'color', [1,0,0], 'LineWidth',linewidth);

uistack(p3f, 'bottom');

delete(p1f);

grid on;
ax = gca;
ax.GridLineStyle = '--';

axis([1, min([3*find(fk2-min_f<1e-12, 1), para.maxits]), 1e-14, 2*max(fk1)]);
ytick = [1e-14, 1e-10, 1e-6, 1e-2, 1e2, 1e6];
set(gca, 'yTick', ytick);

ylb = ylabel({'$\Phi(x_{k})-\Phi(x^\star)$'}, 'FontSize', labelFontSize,...
    'FontAngle', 'normal', 'Interpreter', 'latex');
set(ylb, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
xlb = xlabel({'\vspace{-1.0mm}';'$k$'}, 'FontSize', labelFontSize,...
    'FontAngle', 'normal', 'Interpreter', 'latex');
set(xlb, 'Units', 'Normalized', 'Position', [1/2, -0.075, 0]);


% lg = legend([p1f, p3f, p4f, p2f, p5f, p6f, p7f], ...
%     'Gradient Descent', 'FISTA-BT',...
%     'FISTA-Mod, $p = \frac{1}{20}, q = \frac{1}{2}$', 'Optimal scheme',...
%     'Restarging FISTA', 'Rada-FISTA', 'Greedy FISTA');
lg = legend([p3f, p4f, p2f, p5f, p6f, p7f], ...
    'FISTA-BT',...
    'FISTA-Mod, $p = \frac{1}{20}, q = \frac{1}{2}$', '$\alpha$-FISTA',...
    'Restarging FISTA', 'Rada-FISTA', 'Greedy FISTA');
set(lg,'FontSize', legendFontSize);
% set(lg, 'Location', 'southeast');
set(lg, 'Interpreter', 'latex');
legend('boxoff');

filename = ['results', filesep, sprintf('cmp-lse-fk.pdf')];
print(filename, '-dpdf');
filename = ['results', filesep, sprintf('cmp-lse-fk.png')];
print(filename, '-dpng');


