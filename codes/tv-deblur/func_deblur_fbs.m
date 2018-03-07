function [x, its, ek,ekT, ek_bar,ekT_bar]...
    = func_deblur_fbs(f,h, nu)
%% set up
tol = 1e-8;
maxits = 1e5;

%%%

beta = 1;
gamma = 1.0*beta;

alpha = 2*beta / (4*beta - gamma);
% lambda = 1 /(1.01 *alpha); % lambda \in ]0, 2[
lambda = 1;

tau = nu*gamma;

[m,n] = size(f);
F = fft2(f, m,n);
H = fft2(h, m,n);
conjH = conj(H);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gradF = @(x) ifft2( conjH.*( H.*fft2(x) - F) );
proxJ = @(x, t) perform_prox_tv1D(x, t);

T_fbs = @(x, gamma, tau) ...
    (1-lambda)* x + lambda* proxJ(x-gamma*gradF(x), tau);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = f;
x0 = x;

ek = zeros(1, maxits);
ek_bar = zeros(1,maxits);

lambdak = 0;
lambdaK = zeros(1,maxits);

its = 1;
normE = 1;
while((its<maxits)&&(normE>tol))||(normE>1e10)
    x_old = x;
    
    x = T_fbs(x, gamma, tau);
    %%%%%%% stop?
    normE = norm(x_old - x, 'fro');
    str = sprintf('step %06d: norm(ek) = %.3e', its, normE);
    if its==1 fprintf(str); end
    fprintf('%s%s', char(8*ones(1,length(str))), str); 
    
    ek(its) = normE;

    lambdak = lambdak + lambda;
    lambdaK(its) = lambdak;
    ek_bar(its) = norm( x0-x, 'fro') /lambdak;
    
    its = its + 1;
end
fprintf('\n');

ek = ek(1:its);
ek_bar = ek_bar(1:its);

d0 = norm(x0-x, 'fro');

lambda_ = min(lambda);
ekT = d0 ./ sqrt( (1:its) .* lambda_*(1/alpha-lambda_) );
ekT_bar = 2*d0 ./ lambdaK(1:its);

%%%% EOF

