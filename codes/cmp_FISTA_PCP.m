clear all
close all
clc

addpath data
set(groot,'defaultLineLineWidth',1.5);
%%
type = 'mtx';
type = 'video';

if strcmp(type, 'mtx')
    n = [32, 32] *2;
    
    r = 2; % rank
    xl0 = rand(n(1),r)*diag(r:-1:1)*rand(r,n(2));
    xl0 = xl0/max(abs(xl0(:))) *50;
    
    xs0 = rand(n)*rand(n);
    xs0 = xs0/max(abs(xs0(:))) *50;
    
    ratio = 0.2; % sparsity
    mask = proj_mask(xs0, ratio, 'p');
    xs0 = xs0 .* mask;
    
    f0 = xs0 + xl0;
    sigma = 1e-2*std(f0(:));
    
    f = xs0 + xl0 + sigma*randn(n);
else
    load lobby.mat;
    
    mov = output.mov(1:1:end, 1:1:160, :);
    clear output;
    
    [n1,n2,n3] = size(mov);
    
    f = zeros(n1*n2, n3);
    for i=1:n3
        frame = mov(:,:, i);
        f(:, i) = frame(:);
    end
    f = rescale(f, 0,1);
    
    n = size(f);
end

para.x0 = 0* f;

para.verbose = 1;
%% parameters
para.mu = 1 /sqrt(max(n)); % weight for the sparse part
para.nu = 2; % weight for the low-rank part

% para.beta = 1; % cocoercivity of the gradient
para.gamma = 1;

para.f = f;
para.n = n;

para.tol = 1e-16; % stopping criterion
para.maxits = 1e3; % max # of iteration

proxJ = @(x, t) wthresh(x, 's', t);
gradF = @(x) - ((f-x) - svt(f-x, para.nu));
objPhi = @(x) 1;
%% Lazy-start FISTA-Mod
fprintf(sprintf('performing Lazy-start FISTA...\n'));

p = 1/20;
q = 1/2;
r = 4;

[xs_sol, ~, ~, ~, ~] = func_FISTA_Mod(p,q,r, para, proxJ,gradF, objPhi, 0);

[xs2, its2, dk2, ek2, fk2] = func_FISTA_Mod(p,q,r, para, proxJ,gradF, objPhi, xs_sol);

xl2 = svt(f-xs2, para.nu);

fprintf('\n');
%% FISTA, original
fprintf(sprintf('performing original FISTA...\n'));

p = 1;
q = 1;
r = 4;

[xs1, its1, dk1, ek1, fk1] = func_FISTA_Mod(p,q,r, para, proxJ,gradF, objPhi, xs_sol);

xl1 = svt(f-xs1, para.nu);

fprintf('\n');
%% Restarting FISTA
fprintf(sprintf('performing restarting FISTA...\n'));

p = 1;
q = 1;
r = 4;

[xs3, its3, dk3, ek3, fk3] = func_Restart_FISTA(p,q,r, para, proxJ,gradF, objPhi, xs_sol);

xl3 = svt(f-xs3, para.nu);

fprintf('\n');
%% Rada-FISTA
fprintf(sprintf('performing Rada-FISTA...\n'));

p = 1;
q = 1;
r = 4;

[xs4, its4, dk4, ek4, fk4] = func_Rada_FISTA(p,q,r, para, proxJ,gradF, objPhi, xs_sol);

xl4 = svt(f-xs4, para.nu);

fprintf('\n');
%% Greedy restarting
fprintf(sprintf('performing greedy restarting...\n'));

para.c_gamma = 1.1;
para.a = @(k) 1; %max(2/(1+k/12), 1);

[xs5, its5, dk5, ek5, fk5] = func_Greedy_FISTA(para, proxJ,gradF, objPhi, xs_sol);

xl5 = svt(f-xs5, para.nu);

fprintf('\n');
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

figure(100), clf;
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

axis([1, its1, 1e-8, 2*max(dk1)]);
ytick = [1e-8, 1e-4, 1e-0, 1e4];
set(gca, 'yTick', ytick);

ylb = ylabel({'$\|x_{l,k}-x_{l}^\star\|$'}, 'FontSize', labelFontSize,...
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

filename = ['results', filesep, sprintf('cmp-pcp-dk-%s.pdf', type)];
print(filename, '-dpdf');
filename = ['results', filesep, sprintf('cmp-pcp-dk-%s.png', type)];
print(filename, '-dpng');
%% print images
kk = 30;

if strcmp(type, 'video') && 1
    resolution = 300; % output resolution
    output_size = 300 *[8, 8]; % output size
    
    figure(101), clf;
    set(0,'DefaultAxesFontSize', axesFontSize);
    set(gcf,'paperunits','centimeters','paperposition',[-1.015 -1.13 output_size/resolution]);
    set(gcf,'papersize',output_size/resolution-[1.74 2.49]);
    
    imgsc(reshape(f(:,kk), n1,n2));
    
    filename = ['results', filesep, sprintf('original-frame.pdf')];
    print(filename, '-dpdf');
    filename = ['results', filesep, sprintf('original-frame.png')];
    print(filename, '-dpng');
    
    
    figure(102), clf;
    set(0,'DefaultAxesFontSize', axesFontSize);
    set(gcf,'paperunits','centimeters','paperposition',[-1.015 -1.13 output_size/resolution]);
    set(gcf,'papersize',output_size/resolution-[1.74 2.49]);
    
    imgsc(reshape(xs3(:,kk), n1,n2));
    
    filename = ['results', filesep, sprintf('sparse-component.pdf')];
    print(filename, '-dpdf');
    filename = ['results', filesep, sprintf('sparse-component.png')];
    print(filename, '-dpng');
    
    figure(103), clf;
    set(0,'DefaultAxesFontSize', axesFontSize);
    set(gcf,'paperunits','centimeters','paperposition',[-1.015 -1.13 output_size/resolution]);
    set(gcf,'papersize',output_size/resolution-[1.74 2.49]);
    
    imgsc(reshape(xl3(:,kk), n1,n2));
    
    filename = ['results', filesep, sprintf('lowrank-component.pdf')];
    print(filename, '-dpdf');
    filename = ['results', filesep, sprintf('lowrank-component.png')];
    print(filename, '-dpng');
end