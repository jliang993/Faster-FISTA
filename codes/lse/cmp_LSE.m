clear all;
close all;
clc;
%%
n = 2e1 + 1;

A = 2*eye(n) - diag(ones(n-1,1), -1) - diag(ones(n-1,1), 1);
% A = randn(n-10, n) /sqrt(n);

x_ob = randn(n, 1);
b_ob = A *x_ob;
b = b_ob + 1e-1*randn(n, 1);

para.beta = 1/norm(A)^2;

para.maxits = 1e7 + 1;
para.tol = 1e-10;
para.n = n;

GradF = @(x) (A')*(A*x - b);
ObjF = @(x) norm(A*x-b)^2 /2;

x0 = 1e4*randn(n, 1);

outputType = 'pdf';
%% Gradient Descent
fprintf(sprintf('performing GD...\n'));

[x0, ek0, fk0, its0] = func_GD(x0, para, GradF, ObjF);

fprintf('\n');
%% compute strong convexity
gamma = para.beta;

v = svd((A')*A);
alpha = min(v);

eta = 1 - gamma*alpha;
%% Heavy-ball Method, optimal a
fprintf(sprintf('performing Heavy Ball...\n'));

a_opt = (1-sqrt(1-eta))^2/eta;
a = a_opt - 1e-8;

[x1, ek1, fk1, its1] = func_Heavyball(x0, a, para, GradF, ObjF);

fprintf('\n');
%% FISTA, original
fprintf(sprintf('performing FISTA...\n'));

r = 4;
p = 1;
q = 1;

[x2, ek2, fk2, its2] = func_FISTA_Mod(x0, p,q,r, para, GradF, ObjF);

fprintf('\n');
%% Lazy FISTA-Mod
fprintf(sprintf('performing Lazy-FISTA...\n'));

r = 4;
p = 1/12;
q = 1/2;

[x3, ek3, fk3, its3] = func_FISTA_Mod(x0, p,q,r, para, GradF, ObjF);

fprintf('\n');
%% Restarting FISTA
fprintf(sprintf('performing restarting FISTA...\n'));

r = 4; 
p = 1;
q = p^2;

[x5, ek5, fk5, its5] = func_FISTA_Restart(x0, p,q,r, para, GradF, ObjF);

fprintf('\n');
%% Adaptive FISTA
fprintf(sprintf('performing restarting AdaFISTA...\n'));

r = 4;%*(1 - sqrt(alpha*gamma))^2/(1-alpha*gamma);
p = 1/1.1;
q = p^2;

[x4, ek4, fk4, its4] = func_RAdaFISTA(x0, p,q,r, para, GradF, ObjF);

fprintf('\n');

[its3, its4, its5]

fprintf('\n');
%% relative error ||x_{k}-x_{k-1}|| 
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

p0e = semilogy(ek0, 'color',[0.75,0.0,0.0], 'LineWidth',linewidth);
hold on,
p1e = semilogy(ek1, 'r', 'LineWidth',linewidth);

p2e = semilogy(ek2, 'color',[0.4,0.4,0.4], 'LineWidth',linewidth);
p3e = semilogy(ek3, 'b', 'LineWidth',linewidth);

p4e = semilogy(ek4, 'color',[0.99,0.0,0.99], 'LineWidth',linewidth);
p5e = semilogy(ek5, 'k-.', 'LineWidth',linewidth);

uistack(p2e, 'bottom');

grid on;
ax = gca;
ax.GridLineStyle = '--';

axis([1, length(ek2)/1, 1e-10, 1e2]);
ytick = [1e-10, 1e-6, 1e-2, 1e2];
set(gca, 'yTick', ytick);

ylb = ylabel({'$\|x_{k}-x_{k-1}\|$'}, 'FontSize', labelFontSize,...
    'FontAngle', 'normal', 'Interpreter', 'latex');
set(ylb, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
xlb = xlabel({'\vspace{-1.0mm}';'$k$'}, 'FontSize', labelFontSize,...
    'FontAngle', 'normal', 'Interpreter', 'latex');
set(xlb, 'Units', 'Normalized', 'Position', [1/2, -0.075, 0]);


lg = legend([p0e, p1e, p2e, p3e, p4e, p5e], ...
    'Gradient Descent', 'Optimal Heavyball',...
    'FISTA-BT', 'Lazy FISTA-Mod, $p = \frac{1}{20}, q = \frac{1}{1}$',...
    'Ada-FISTA+Restart',...
    'Restarting FISTA');
set(lg,'FontSize', legendFontSize);
set(lg, 'Interpreter', 'latex');
legend('boxoff');

epsname = sprintf('cmp_lse_ek.%s', outputType);
if strcmp(outputType, 'png')
    print(epsname, '-dpng');
else
    print(epsname, '-dpdf');
end
%% plot Phi(x_{k}) - Phi(x*)
phistar = min([min(fk0), min(fk1), min(fk2), min(fk3), min(fk4)]);
linewidth = 1;

axesFontSize = 8;
labelFontSize = 8;
legendFontSize = 8;

resolution = 300; % output resolution
output_size = 300 *[10, 8]; % output size

%%%%%% relative error

figure(102), clf;
set(0,'DefaultAxesFontSize', axesFontSize);
set(gcf,'paperunits','centimeters','paperposition',[-0.1 -0.0 output_size/resolution]);
set(gcf,'papersize',output_size/resolution-[0.85 0.4]);

p0e = semilogy(fk0-phistar, 'color',[0.75,0.0,0.0], 'LineWidth',linewidth);
hold on,
p1e = semilogy(fk1-phistar, 'r', 'LineWidth',linewidth);

p2e = semilogy(fk2-phistar, 'color',[0.4,0.4,0.4], 'LineWidth',linewidth);
p3e = semilogy(fk3-phistar, 'b', 'LineWidth',linewidth);

p4e = semilogy(fk4-phistar, 'color',[0.99,0.0,0.99], 'LineWidth',linewidth);
p5e = semilogy(fk5-phistar, 'k-.', 'LineWidth',linewidth);



uistack(p2e, 'bottom');

grid on;
ax = gca;
ax.GridLineStyle = '--';

axis([1, length(ek2)/1, 1e-14, 1e6]);
ytick = [1e-14, 1e-10, 1e-6, 1e-2, 1e2, 1e6];
set(gca, 'yTick', ytick);

ylb = ylabel({'$\Phi(x_{k})-\Phi(x^\star)$'}, 'FontSize', labelFontSize,...
    'FontAngle', 'normal', 'Interpreter', 'latex');
set(ylb, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
xlb = xlabel({'\vspace{-1.0mm}';'$k$'}, 'FontSize', labelFontSize,...
    'FontAngle', 'normal', 'Interpreter', 'latex');
set(xlb, 'Units', 'Normalized', 'Position', [1/2, -0.075, 0]);


lg = legend([p0e, p1e, p2e, p3e, p4e, p5e], ...
    'Gradient Descent', 'Optimal Heavyball',...
    'FISTA-BT', 'Lazy FISTA-Mod, $p = \frac{1}{20}, q = \frac{1}{1}$',...
    'Ada-FISTA+Restart',...
    'Restarting FISTA');
set(lg,'FontSize', legendFontSize);
set(lg, 'Interpreter', 'latex');
legend('boxoff');

epsname = sprintf('cmp_lse_fk.%s', outputType);
if strcmp(outputType, 'png')
    print(epsname, '-dpng');
else
    print(epsname, '-dpdf');
end
%% plot Phi(x_{k}) - Phi(x*), EPSRC
phistar = min([min(fk0), min(fk1), min(fk2)]);
linewidth = 1;

axesFontSize = 10;
labelFontSize = 12;
legendFontSize = 10;

resolution = 300; % output resolution
output_size = 300 *[12, 8]; % output size

%%%%%% relative error

figure(102), clf;
set(0,'DefaultAxesFontSize', axesFontSize);
set(gcf,'paperunits','centimeters','paperposition',[-0.1 -0.0 output_size/resolution]);
set(gcf,'papersize',output_size/resolution-[0.85 0.4]);

gap = 5;
N = numel(fk0);
p0e = semilogy(1:gap:N, (fk0(1:gap:N)-phistar)/1e2, 'color',[0.99,0.0,0.0], 'LineWidth',linewidth);
hold on,

p2e = semilogy(1:gap:N, fk2(1:gap:N)-phistar, 'color',[0.1,0.1,0.1], 'LineWidth',linewidth);

p1 = semilogy(1:gap:N, 5e6./(1:gap:N), '--', 'color',[0.1,0.1,0.99], 'LineWidth',linewidth);
p2 = semilogy(1:gap:N, 1e8./(1:gap:N).^2, '-.', 'color',[0.1,0.1,0.99], 'LineWidth',linewidth);

uistack(p2e, 'bottom');

grid on;
ax = gca;
ax.GridLineStyle = '--';

axis([1, para.maxits/gap, 1e-9, 1e4]);
ytick = [1e-14, 1e-10, 1e-6, 1e-2, 1e2, 1e6];
set(gca, 'yTick', ytick);

ylb = ylabel({'$F(x_{k})-F(x^\star)$'}, 'FontSize', labelFontSize,...
    'FontAngle', 'normal', 'Interpreter', 'latex');
set(ylb, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
xlb = xlabel({'\vspace{-0.7mm}';'$k$'}, 'FontSize', labelFontSize,...
    'FontAngle', 'normal', 'Interpreter', 'latex');
set(xlb, 'Units', 'Normalized', 'Position', [1/2, -0.075, 0]);


lg = legend([p1, p0e, p2, p2e], ...
    '$O(1/k)$', 'Gradient descent',...
    '$O(1/k^2)$', 'Nesterov''s scheme');
set(lg,'FontSize', legendFontSize);
set(lg, 'Interpreter', 'latex');
% set(lg, 'Location', 'Best');
legend('boxoff');

pos = get(lg, 'Position');
set(lg, 'Position', [pos(1)-0.05, pos(2)-0.25, pos(3:4)]);
pos_ = get(lg, 'Position');
legend('boxoff');


epsname = sprintf('epsrc_lse_fk.%s', outputType);
if strcmp(outputType, 'png')
    print(epsname, '-dpng');
else
    print(epsname, '-dpdf');
end
















