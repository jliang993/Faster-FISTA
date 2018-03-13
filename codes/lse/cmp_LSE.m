clear all;
close all;
clc;
%%
n = 2e2 + 1;

A = 2*eye(n) - diag(ones(n-1,1), -1) - diag(ones(n-1,1), 1);

x_ob = randn(n, 1);
b_ob = A *x_ob;
b = b_ob + 1e-1*randn(n, 1);

para.beta = 1/norm(A)^2;

para.maxits = 1e6 + 1;
para.tol = 1e-10;
para.n = n;

GradF = @(x) (A')*(A*x - b);
ObjF = @(x) norm(A*x-b)^2 /2;

outputType = 'png';
%% Gradient Descent
fprintf(sprintf('performing GD...\n'));

[x0, ek0, fk0, its0] = func_GD(para, GradF, ObjF);

fprintf('\n');
%% compute strong convexity
gamma = para.beta;

v = svd((A')*A);
alpha = min(v);

eta = 1 - gamma*alpha;
%% Heavy-ball Method, optimal a
fprintf(sprintf('performing Heavy Ball...\n'));

a_opt = (1-sqrt(1-eta))^2/eta;
a = a_opt - 1e-6;

[x1, ek1, fk1, its1] = func_Heavyball(a, para, GradF, ObjF);

fprintf('\n');
%% FISTA, original
fprintf(sprintf('performing FISTA...\n'));

r = 4;
p = 1;
q = 1;

[x2, ek2, fk2, its2] = func_FISTA_Mod(p,q,r, para, GradF, ObjF);

fprintf('\n');
%% Lazy FISTA-Mod
fprintf(sprintf('performing Lazy-FISTA...\n'));

r = 4;
p = 1/12;
q = 1/2;

[x3, ek3, fk3, its3] = func_FISTA_Mod(p,q,r, para, GradF, ObjF);

fprintf('\n');
%%
d = 20;
[x3c, ek3c, fk3c, its3c] = func_FISTA_CD(d, para, GradF, ObjF);

fprintf('\n');
%% Adaptive FISTA
fprintf(sprintf('performing adaptive FISTA...\n'));

r = 4*(1 - sqrt(alpha*gamma))^2/(1-alpha*gamma);
p = 1;
q = 1;

[x4, ek4, fk4, its4] = func_AdaFISTA(p,q,r, para, GradF, ObjF);

fprintf('\n');
%% relative error ||x_{k}-x_{k-1}|| 
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

p0e = semilogy(ek0, 'm', 'LineWidth',linewidth);
hold on,
p1e = semilogy(ek1, 'r', 'LineWidth',linewidth);

p2e = semilogy(ek2, 'b', 'LineWidth',linewidth);
p3e = semilogy(ek3, 'k', 'LineWidth',linewidth);

p4e = semilogy(ek4, 'k-.', 'LineWidth',linewidth);

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


lg = legend([p0e, p1e, p2e, p3e, p4e], ...
    'Gradient Descent', 'Optimal Heavyball',...
    'FISTA-BT', 'Lazy FISTA-Mod, $p = \frac{1}{20}, q = \frac{1}{1}$',...
    'Ada-FISTA, $p=1, q=1$');
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

figure(101), clf;
set(0,'DefaultAxesFontSize', axesFontSize);
set(gcf,'paperunits','centimeters','paperposition',[-0.1 -0.0 output_size/resolution]);
set(gcf,'papersize',output_size/resolution-[0.85 0.4]);

p0e = semilogy(fk0-phistar, 'm', 'LineWidth',linewidth);
hold on,
p1e = semilogy(fk1-phistar, 'r', 'LineWidth',linewidth);

p2e = semilogy(fk2-phistar, 'b', 'LineWidth',linewidth);
p3e = semilogy(fk3-phistar, 'k', 'LineWidth',linewidth);

p4e = semilogy(fk4-phistar, 'k-.', 'LineWidth',linewidth);

uistack(p2e, 'bottom');

grid on;
ax = gca;
ax.GridLineStyle = '--';

axis([1, length(ek2)/1, 1e-14, 1e2]);
ytick = [1e-14, 1e-10, 1e-6, 1e-2, 1e2];
set(gca, 'yTick', ytick);

ylb = ylabel({'$\Phi(x_{k})-\Phi(x^\star)$'}, 'FontSize', labelFontSize,...
    'FontAngle', 'normal', 'Interpreter', 'latex');
set(ylb, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
xlb = xlabel({'\vspace{-1.0mm}';'$k$'}, 'FontSize', labelFontSize,...
    'FontAngle', 'normal', 'Interpreter', 'latex');
set(xlb, 'Units', 'Normalized', 'Position', [1/2, -0.075, 0]);


lg = legend([p0e, p1e, p2e, p3e, p4e], ...
    'Gradient Descent', 'Optimal Heavyball',...
    'FISTA-BT', 'Lazy FISTA-Mod, $p = \frac{1}{20}, q = \frac{1}{1}$',...
    'Ada-FISTA, $p=1, q=1$');
set(lg,'FontSize', legendFontSize);
set(lg, 'Interpreter', 'latex');
legend('boxoff');

epsname = sprintf('cmp_lse_fk.%s', outputType);
if strcmp(outputType, 'png')
    print(epsname, '-dpng');
else
    print(epsname, '-dpdf');
end















