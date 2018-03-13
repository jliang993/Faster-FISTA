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

outputType = 'png';
%% Original FISTA-BT
fprintf(sprintf('performing FISTA-BT...\n'));

r = 4;

p = 1;
q = 1;
[x, its, ek, phik] = func_FISTA_Mod(para, gradF, proxJ, objPhi, J, p,q,r);

fprintf('\n');
%% FISTA-Mod
K = 20;

p = 1/50;
q = 1/1e1;

tic;
for i=1:K
    [x_m2, its_m2, ek_m2, phik_m2] = func_FISTA_Mod(para, gradF, proxJ, objPhi, J, p,q,r);
end
toc/K

fprintf('\n');
%% Adaptive-FISTA, AdaFISTA
fprintf(sprintf('performing Ada-FISTA...\n'));
r = 4;
p = 1;
q = p^2;

tic
for i=1:K
    [x_a, its_a, ek_a, phik_a] = func_AdaFISTA_s1(para, gradF, proxJ, objPhi, J, p,q,r);
end
toc/K

fprintf('\n');
%% Restarting FISTA
fprintf(sprintf('performing restarting FISTA...\n'));
r = 4;

p = 1;
q = 1;

tic
for i=1:K
    [x_r, its_r, ek_r, phik_r] = func_FISTA_Restart(para, gradF, proxJ, objPhi, J, p,q,r);
end
toc/K

fprintf('\n');