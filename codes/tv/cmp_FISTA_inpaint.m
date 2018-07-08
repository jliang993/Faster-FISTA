clear all
close all
clc
%%
name = 'cameraman.png';

u = double(imread(name));
u = u(1:2:end, 1:2:end);
[m, n] = size(u);

A = proj_mask(u, 0.5, 'p');

% f0 = imfilter(u, h_, 'circular');
f0 = A.* u;

var = 1e0;
f = f0 + var * randn(m,n);

mu = 8;

para.f = f;
para.mu = mu;
para.beta = 1;

para.tol = 1e-7;
para.maxits = 1e4;

gradF = @(x) A.*( A.*x - f);
proxJ = @(x, t) perform_prox_tv1D(x, t);

dxf = @(x) [diff(x,1,2), zeros(m,1)];
dyf = @(x) [diff(x,1,1); zeros(1,n)];

objF = @(x) mu*( sum(sum( abs(dxf(x)) )) + sum(sum( abs(dyf(x)) )) )...
    + 1/2*norm(A(:).*x(:) - f(:), 'fro')^2;

% objF = @(x) 1/2*norm(A(:).*x(:) - f(:), 'fro')^2;
%% FISTA-BT
fprintf(sprintf('performing FISTA-BT...\n'));

r = 4;

p = 1;
q = 1;

[x1, its1, ek1, fk1, sk1] = func_FISTA_Mod(gradF, proxJ, objF, para, p,q,r);
fprintf('\n');
%% FISTA-Mod
fprintf(sprintf('performing FISTA-Mod...\n'));

p = 1/30;
q = 1/10;

[x2, its2, ek2, fk2, sk2] = func_FISTA_Mod(gradF, proxJ, objF, para, p,q,r);

fprintf('\n');
%% Adaptive + Restarting
fprintf(sprintf('performing restarting AdaFISTA...\n'));

r = 4;
p = 1/2;
q = p^2;

[x_ra, its_ra, ek_ra, fk_ra] = func_RAdaFISTA(gradF, proxJ, objF, para, p,q,r);

fprintf('\n');
%% Restarting FISTA
fprintf(sprintf('performing restarting FISTA...\n'));

r = 4;

p = 1;
q = 1;

[x_r, its_r, ek_r, fk_r] = func_FISTA_Restart(gradF, proxJ, objF, para, p,q,r);

fprintf('\n');
%% plot \Phi(\xk)-\Phi(\xsol)
fsol = min( [min(fk1), min(fk2), min(fk_ra), min(fk_r)] );
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
p1 = semilogy(fk1-fsol, 'Color',grey1, 'LineWidth',linewidth);
hold on,

blue1 = [0.12,0.48,1.0];
p2 = semilogy(fk2-fsol, 'Color',blue1, 'LineWidth',linewidth);

pasr = semilogy(fk_ra-fsol, 'Color',[0.99,0.01,0.99], 'LineWidth',linewidth);

pr = semilogy(fk_r-fsol, '-.', 'Color',[0.33,0.33,0.33], 'LineWidth',linewidth);

grid on;
ax = gca;
ax.GridLineStyle = '--';

% v = axis;
axis([1 length(ek1)/4 1e-10 1e2]);
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


epsname = sprintf('cmp_fista_tvinpaint_objf.pdf');
print(epsname, '-dpdf');
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

pasr = semilogy(ek_ra, 'Color',[0.99,0.01,0.99], 'LineWidth',linewidth);

pr = semilogy(ek_r, '-.', 'Color',[0.33,0.33,0.33], 'LineWidth',linewidth);

grid on;
ax = gca;
ax.GridLineStyle = '--';

% v = axis;
axis([1 length(ek1)/4 1e-6 1e2]);
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


epsname = sprintf('cmp_fista_tvinpaint_ek.pdf');
print(epsname, '-dpdf');
%% print images
% close all;
% 
% resolution = 300; % output resolution
% output_size = 300 *[8, 8]; % output size
% 
% figure(101), clf;
% set(0,'DefaultAxesFontSize', axesFontSize);
% set(gcf,'paperunits','centimeters','paperposition',[-1.015 -1.01 output_size/resolution]);
% set(gcf,'papersize',output_size/resolution-[1.74 1.75]);
% 
% imgsc(u);
% 
% epsname = sprintf('original-img.pdf');
% print(epsname, '-dpdf');
% 
% 
% figure(102), clf;
% set(0,'DefaultAxesFontSize', axesFontSize);
% set(gcf,'paperunits','centimeters','paperposition',[-1.015 -1.01 output_size/resolution]);
% set(gcf,'papersize',output_size/resolution-[1.74 1.75]);
% 
% imgsc(f);
% 
% epsname = sprintf('original-miss.pdf');
% print(epsname, '-dpdf');
% 
% figure(103), clf;
% set(0,'DefaultAxesFontSize', axesFontSize);
% set(gcf,'paperunits','centimeters','paperposition',[-1.015 -1.01 output_size/resolution]);
% set(gcf,'papersize',output_size/resolution-[1.74 1.75]);
% 
% imgsc(x2);
% 
% epsname = sprintf('original-inpaint.pdf');
% print(epsname, '-dpdf');