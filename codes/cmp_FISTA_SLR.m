clear all
% close all
clc

addpath data
addpath toolbox
set(groot,'defaultLineLineWidth',1.5);
%% load and scale data
load('australian_label.mat');
load('australian_sample.mat');

h = full(h);
% rescale the data
for j=1:size(h,2)
    h(:,j) = rescale(h(:,j), -1, 1);
end
%% parameters
[m, n] = size(h);

para.m = m;
para.n = n;

para.W = h;
para.y = l;

para.mu = 1e-2;

Li = zeros(m, 1);
for i=1:m
    Wi = para.W(i,:);
    Li(i) = norm(Wi)^2 /4;
end
para.beta_fi = 1 /max(Li);

para.tol = 1e-15; % stopping criterion
para.maxits = 1e5; % max # of iteration

para.beta = 4*m/norm(para.W)^2;
para.gamma = para.beta;

gradF = @(x) grad_logistic(x, para.W, para.y) /m;% + para.mu* x;
% ProxJ = @(x, t) [wthresh(x(1:end-1), 's', t); x(end)];
proxJ = @(x, t) wthresh(x, 's', t);
objPhi = @(x) 1;

para.x0 = zeros(n, 1);

para.verbose = 1;
%% Lazy-start Para-FISTA     
fprintf(sprintf('performing Lazy-start FISTA...\n'));

p = 1/20;
q = 1/2;
r = 4;

[xsol] = func_FISTA_Mod(p,q,r, para, proxJ,gradF, objPhi, 0);
[x2, its2, dk2, ek2] = func_FISTA_Mod(p,q,r, para, proxJ,gradF, objPhi, xsol);

fprintf('\n');
%% FISTA, original   
fprintf(sprintf('performing original FISTA...\n'));

p = 1;
q = 1;
r = 4;

[x1, its1, dk1, ek1] = func_FISTA_Mod(p,q,r, para, proxJ,gradF, objPhi, xsol);

fprintf('\n');
%% Restarting FISTA         
fprintf(sprintf('performing restarting FISTA...\n'));

p = 1;
q = 1;
r = 4;

[x3, its3, dk3, ek3] = func_Restart_FISTA(p,q,r, para, proxJ,gradF, objPhi, xsol);

fprintf('\n');
%% Rada-FISTA       
fprintf(sprintf('performing Rada-FISTA...\n'));

p = 1;
q = 1;
r = 4;

[x4, its4, dk4, ek4] = func_Rada_FISTA(p,q,r, para, proxJ,gradF, objPhi, xsol);

fprintf('\n');

disp([its1, its2, its3, its4])
%% Greedy FISTA     
fprintf(sprintf('performing Greedy FISTA...\n'));

para.c_gamma = 1.25;
para.a = @(k) max(2/(1+k/12), 1);

[x5, its5, dk5, ek5] = func_Greedy_FISTA(para, proxJ,gradF, objPhi, xsol);

disp([its1, its2, its3, its4, its5])

fprintf('\n');
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

p1d = semilogy(dk1, 'color', [1/1.75,1/1/1.75,1/1/1.75], 'LineWidth',linewidth);
hold on,

p2d = semilogy(dk2, 'color', [1/6,1/6,1/6], 'LineWidth',linewidth);

p3d = semilogy(dk3, 'color', [1/6,1/1.4,1/6], 'LineWidth',linewidth);

p4d = semilogy(dk4, 'color', [1/6,1/4,1/1.2], 'LineWidth',linewidth);

p5d = semilogy(dk5, 'color', [1,0,0], 'LineWidth',linewidth);

uistack(p1d, 'bottom');

grid on;
ax = gca;
ax.GridLineStyle = '--';

axis([1, 3*its3, 1e-10, 2*max(dk1)]);
ytick = [1e-8, 1e-4, 1e-0, 1e4];
set(gca, 'yTick', ytick);

ylb = ylabel({'$\|x_{k}-x^\star\|$'}, 'FontSize', labelFontSize,...
    'FontAngle', 'normal', 'Interpreter', 'latex');
set(ylb, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
xlb = xlabel({'\vspace{-1.0mm}';'$k$'}, 'FontSize', labelFontSize,...
    'FontAngle', 'normal', 'Interpreter', 'latex');
set(xlb, 'Units', 'Normalized', 'Position', [1/2, -0.075, 0]);


lg = legend([p1d, p2d, p3d, p4d], ...
    'FISTA-BT',...
    'Para-FISTA, $p = \frac{1}{20}, q = \frac{1}{2}$',...
    'Restarging FISTA', 'Rada-FISTA');
set(lg,'FontSize', legendFontSize);
set(lg, 'Interpreter', 'latex');
% set(lg, 'Location', 'best');
legend('boxoff');

filename = ['results', filesep, sprintf('cmp-slr-dk.pdf')];
print(filename, '-dpdf');
filename = ['results', filesep, sprintf('cmp-slr-dk.png')];
print(filename, '-dpng');
