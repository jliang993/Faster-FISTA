clear all
close all
clc
%%
n = [32, 32];

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
%% parameters
para.mu1 = 1 /sqrt(max(n)); % weight for the sparse part
para.mu2 = 2; % weight for the low-rank part

para.beta = 1; % cocoercivity of the gradient

para.f = f;
para.n = n;

para.tol = 1e-10; % stopping criterion
para.maxits = 1e5; % max # of iteration

GradF = @(x) - ((f-x) - svt(f-x, para.mu2));
ProxJ = @(x, t) wthresh(x, 's', t);
% %
% [xs, xl, its, ek] = func_tFB_pcp(para, GradF, ProxJ);

objF = @(l,s) para.mu1*sum(abs(s(:))) + para.mu2*sum(svd(l)) + 1/2*norm(f(:)-s(:)-l(:))^2;
%% FISTA-Mod
fprintf(sprintf('performing FISTA-Mod...\n'));

r = 4;

p = 1/1;
q = 1/1;

[xs1,xl1, its1, ek1, fk1] = func_FISTA_Mod(para, GradF, ProxJ, p,q,r, objF);

p = 1/5;
q = 1/1;

[xs2,xl2, its2, ek2, fk2] = func_FISTA_Mod(para, GradF, ProxJ, p,q,r, objF);


p = 1/30;
q = 1/10;

[xs3,xl3, its3, ek3, fk3] = func_FISTA_Mod(para, GradF, ProxJ, p,q,r, objF);

fprintf('\n');
%% RAd-FISTA
fprintf(sprintf('performing restarting RAda-FISTA...\n'));
r = 4;
p = 1/1.5;
q = p^2;

[xs_ra,xl_ra, its_ra, ek_ra, fk_ra] = func_RAdaFISTA(para, GradF, ProxJ, p,q,r, objF);

fprintf('\n');
%% Restarting FISTA
fprintf(sprintf('performing restarting FISTA...\n'));

r = 4;

p = 1;
q = 1;

[xsr,xlr, its_r, ek_r, fk_r] = func_FISTA_Restart(para, GradF, ProxJ, p,q,r, objF);

fprintf('\n');
%% plot Phi(\xk)-\Phi(\xsol)
fsol = min( [min(fk1), min(fk2), min(fk3), min(fk_ra), min(fk_r)] );
linewidth = 1;

axesFontSize = 8;
labelFontSize = 8;
legendFontSize = 8;

resolution = 300; % output resolution
output_size = 300 *[10, 8]; % output size

%%%%%% relative error

figure(100), clf;
set(0,'DefaultAxesFontSize', axesFontSize);
set(gcf,'paperunits','centimeters','paperposition',[-0.1 -0.17 output_size/resolution]);
set(gcf,'papersize',output_size/resolution-[0.85 0.5]);

grey1 = [0.3,0.3,0.3];
p1 = semilogy(fk1 - fsol, 'Color',grey1, 'LineWidth',linewidth);
hold on,

% blue1 = [0.12,0.48,1.0];
% p2 = semilogy(fk2 - fsol, 'Color',blue1, 'LineWidth',linewidth);
blue2 = [0.9,0.0,0.0];
p3 = semilogy(fk3 - fsol, 'Color',blue2, 'LineWidth',linewidth);

pasr = semilogy(fk_ra - fsol, 'Color',[0.99,0.01,0.99], 'LineWidth',linewidth);
pr = semilogy(fk_r - fsol, '-.', 'Color',[0.33,0.33,0.33], 'LineWidth',linewidth);

grid on;
ax = gca;
ax.GridLineStyle = '--';

% v = axis;
axis([1 length(ek1)/2 1e-10 1e4]);
ytick = [1e-10, 1e-6, 1e-2, 1e2];
set(gca, 'yTick', ytick);

ylb = ylabel({'$\Phi(x_{k})-\Phi(x_{k-1})$'}, 'FontSize', labelFontSize,...
    'FontAngle', 'normal', 'Interpreter', 'latex');
set(ylb, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
xlb = xlabel({'\vspace{-1.0mm}';'$k$'}, 'FontSize', labelFontSize,...
    'FontAngle', 'normal', 'Interpreter', 'latex');
set(xlb, 'Units', 'Normalized', 'Position', [1/2, -0.055, 0]);

lg = legend([p1, p3, pasr, pr], 'FISTA-BT',...
    'FISTA-Mod, $p = \frac{1}{30}, q = \frac{1}{10}$',...
    'Ada-FISTA+Restart, $p=\frac{1}{1.5}, q=1$',...
    'Restarting FISTA');
set(lg,'FontSize', legendFontSize);
set(lg, 'Interpreter', 'latex');
%
pos = get(lg, 'Position');
% set(lg, 'Position', [pos(1)-0.125, pos(2)-0.075, pos(3:4)]);
% pos_ = get(lg, 'Position');
legend('boxoff');


epsname = sprintf('cmp_fista_pcp_mtx_objf.png');
print(epsname, '-dpng');
%% plot ek
linewidth = 1;

axesFontSize = 8;
labelFontSize = 8;
legendFontSize = 8;

resolution = 300; % output resolution
output_size = 300 *[10, 8]; % output size

%%%%%% relative error

figure(101), clf;
set(0,'DefaultAxesFontSize', axesFontSize);
set(gcf,'paperunits','centimeters','paperposition',[-0.1 -0.17 output_size/resolution]);
set(gcf,'papersize',output_size/resolution-[0.85 0.5]);

grey1 = [0.3,0.3,0.3];
p1 = semilogy(ek1, 'Color',grey1, 'LineWidth',linewidth);
hold on,

blue1 = [0.12,0.48,1.0];
p2 = semilogy(ek2, 'Color',blue1, 'LineWidth',linewidth);
blue2 = [0.9,0.0,0.0];
p3 = semilogy(ek3, 'Color',blue2, 'LineWidth',linewidth);

pasr = semilogy(ek_ra, 'Color',[0.99,0.01,0.99], 'LineWidth',linewidth);
pr = semilogy(ek_r, '-.', 'Color',[0.33,0.33,0.33], 'LineWidth',linewidth);

grid on;
ax = gca;
ax.GridLineStyle = '--';

% v = axis;
axis([1 length(ek1)/2 1e-10 1e2]);
ytick = [1e-10, 1e-6, 1e-2, 1e2];
set(gca, 'yTick', ytick);

ylb = ylabel({'$\|x_{k}-x_{k-1}\|$'}, 'FontSize', labelFontSize,...
    'FontAngle', 'normal', 'Interpreter', 'latex');
set(ylb, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
xlb = xlabel({'\vspace{-1.0mm}';'$k$'}, 'FontSize', labelFontSize,...
    'FontAngle', 'normal', 'Interpreter', 'latex');
set(xlb, 'Units', 'Normalized', 'Position', [1/2, -0.055, 0]);

lg = legend([p1, p2, p3, pasr, pr], 'FISTA-BT',...
    'FISTA-Mod, $p = \frac{1}{5}, q = {1}$',...
    'FISTA-Mod, $p = \frac{1}{30}, q = \frac{1}{10}$',...
    'Ada-FISTA+Restart, $p=\frac{1}{1.5}, q=1$',...
    'Restarting FISTA');
set(lg,'FontSize', legendFontSize);
set(lg, 'Interpreter', 'latex');
%
pos = get(lg, 'Position');
% set(lg, 'Position', [pos(1)-0.125, pos(2)-0.075, pos(3:4)]);
% pos_ = get(lg, 'Position');
legend('boxoff');


epsname = sprintf('cmp_fista_pcp_mtx_ek.png');
print(epsname, '-dpng');
%% print images
% resolution = 300; % output resolution
% output_size = 300 *[8, 8]; % output size
% 
% figure(101), clf;
% set(0,'DefaultAxesFontSize', axesFontSize);
% set(gcf,'paperunits','centimeters','paperposition',[-1.015 -1.01 output_size/resolution]);
% set(gcf,'papersize',output_size/resolution-[1.74 1.75]);
% 
% imgsc(f);
% 
% epsname = sprintf('observation.png');
% print(epsname, '-dpng');
% 
% 
% figure(102), clf;
% set(0,'DefaultAxesFontSize', axesFontSize);
% set(gcf,'paperunits','centimeters','paperposition',[-1.015 -1.01 output_size/resolution]);
% set(gcf,'papersize',output_size/resolution-[1.74 1.75]);
% 
% imgsc(xs3);
% 
% epsname = sprintf('sparse-mtx.png');
% print(epsname, '-dpng');
% 
% figure(103), clf;
% set(0,'DefaultAxesFontSize', axesFontSize);
% set(gcf,'paperunits','centimeters','paperposition',[-1.015 -1.01 output_size/resolution]);
% set(gcf,'papersize',output_size/resolution-[1.74 1.75]);
% 
% imgsc(xl3);
% 
% epsname = sprintf('lowrank-mtx.png');
% print(epsname, '-dpng');