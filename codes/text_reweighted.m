clear all
close all
clc

addpath toolbox
set(groot,'defaultLineLineWidth',1.5);
%% problem set up
J = 'lasso';

[para, gradF,proxJ, objPhi] = problem_FB(J);

para.J = J;
para.tol = 1e-15;
para.maxits = 5e4 + 1;
para.gamma = para.beta;

para.x0 = zeros(para.n, 1);

para.verbose = 1;
%% Lazy-start FISTA-Mod     
fprintf(sprintf('performing Lazy-start FISTA...\n'));

[xsol, ~, ~, ~] = func_FB(para, proxJ,gradF, 0);

[x2, its2, dk2, ek2] = func_FB(para, proxJ,gradF, xsol);

fprintf('\n');