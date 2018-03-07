clear all
close all
clc
%%
name = 'cameraman.png';

u = double(imread(name));
u = u(1:2:end, 1:2:end);
[m, n] = size(u);

type = 'gaussian'; hsize = 11; sigma = 1;
% type = 'motion'; hsize = 5; sigma = 1;
h_ = fspecial(type, hsize, sigma);

h = zeros(m,n);
h(1:size(h_,1),1:size(h_,2)) = h_;
h = circshift(h,-floor(size(h_)/2));

% f0 = imfilter(u, h_, 'circular');
f0 = ifft2( fft2(u) .* fft2(h,m,n) );

var = 1e0;
f = f0 + var * randn(m,n);

mu = 1;
%% FISTA-Mod
fprintf(sprintf('performing FISTA-Mod...\n'));

r = 4;

p = 1/1;
q = 1/1;

[x1, its1, ek1] = func_FISTA_Mod(f,h, mu, p,q,r);


p = 1/5;
q = 1/1;

[x2, its2, ek2] = func_FISTA_Mod(f,h, mu, p,q,r);
%%
p = 1/30;
q = 1/10;

[x3, its3, ek3] = func_FISTA_Mod(f,h, mu, p,q,r);

fprintf('\n');
%% plot
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
p1 = semilogy(ek1, 'Color',grey1, 'LineWidth',linewidth);
hold on,

blue1 = [0.12,0.48,1.0];
p2 = semilogy(ek2, 'Color',blue1, 'LineWidth',linewidth);
blue2 = [0.9,0.0,0.0];
p3 = semilogy(ek3, 'Color',blue2, 'LineWidth',linewidth);

grid on;
ax = gca;
ax.GridLineStyle = '--';

% v = axis;
axis([1 length(ek1)+0 1e-10 1e2]);
ytick = [1e-10, 1e-6, 1e-2, 1e2];
set(gca, 'yTick', ytick);

ylb = ylabel({'$\|x_{k}-x_{k-1}\|$'}, 'FontSize', labelFontSize,...
    'FontAngle', 'normal', 'Interpreter', 'latex');
set(ylb, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
xlb = xlabel({'\vspace{-1.0mm}';'$k$'}, 'FontSize', labelFontSize,...
    'FontAngle', 'normal', 'Interpreter', 'latex');
set(xlb, 'Units', 'Normalized', 'Position', [1/2, -0.055, 0]);

lg = legend([p1, p2, p3], 'FISTA-BT',...
    'FISTA-Mod, $p = \frac{1}{5}, q = {1}$',...
    'FISTA-Mod, $p = \frac{1}{30}, q = \frac{1}{10}$');
set(lg,'FontSize', 10);
set(lg, 'Interpreter', 'latex');
%
pos = get(lg, 'Position');
set(lg, 'Position', [pos(1)-0.125, pos(2)-0.075, pos(3:4)]);
pos_ = get(lg, 'Position');
legend('boxoff');


epsname = sprintf('cmp_fista_tvdeblur.png');
print(epsname, '-dpng');
%% print images
close all;

resolution = 300; % output resolution
output_size = 300 *[8, 8]; % output size

figure(101), clf;
set(0,'DefaultAxesFontSize', axesFontSize);
set(gcf,'paperunits','centimeters','paperposition',[-1.015 -1.01 output_size/resolution]);
set(gcf,'papersize',output_size/resolution-[1.74 1.75]);

imgsc(u);

epsname = sprintf('original-img.png');
print(epsname, '-dpng');


figure(102), clf;
set(0,'DefaultAxesFontSize', axesFontSize);
set(gcf,'paperunits','centimeters','paperposition',[-1.015 -1.01 output_size/resolution]);
set(gcf,'papersize',output_size/resolution-[1.74 1.75]);

imgsc(f);

epsname = sprintf('original-blur.png');
print(epsname, '-dpng');

figure(103), clf;
set(0,'DefaultAxesFontSize', axesFontSize);
set(gcf,'paperunits','centimeters','paperposition',[-1.015 -1.01 output_size/resolution]);
set(gcf,'papersize',output_size/resolution-[1.74 1.75]);

imgsc(x3);

epsname = sprintf('original-deblur.png');
print(epsname, '-dpng');