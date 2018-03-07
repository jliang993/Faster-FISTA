function [x, its, ek,ekT, ek_bar,ekT_bar] = func_deblur_gfb(f,h, nu)
%%%%
tol = 1e-10;
maxits = 1e4;

w1 = 1/2;
w2 = 1 - w1;

beta = 1;
gamma = 1.0*beta;

alpha = 2*beta / (4*beta - gamma);
lambda = ones(1,maxits);

tau = nu*gamma/w1;

[m,n] = size(f);
F = fft2(f, m,n);
H = fft2(h, m,n);
conjH = conj(H);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prox1 = @(x, t) perform_prox_tv1D(x, t); 
prox2 = @(x) max(min(x,255), 0); % box constraint

grad = @(x) ifft2( conjH.*( H.*fft2(x) - F) ); % gradient

proj = @(z1,z2) [w1*z1 + w2*z2]; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z1 = zeros(size(f));
z2 = z1;

x = proj(z1, z2);

z10 = z1;
z20 = z2;

ek = zeros(1, maxits);
ek_bar = zeros(1, maxits);

lambdak = 0;
lambdaK = zeros(1, maxits);

its = 1;
while(its<maxits)
        
    z1_old = z1;
    z2_old = z2;
    
    % forward step
    gamma_g = gamma* grad(x);
    
    % J_A1, SVT
    u1 = prox1(2*x-z1-gamma_g, tau);
    z1 = z1 + lambda(its)*( u1 - x );
    
    % J_A2, Non-negative constraint
    u2 = prox2(2*x-z2-gamma_g);
    z2 = z2 + lambda(its)*( u2 - x );
    
    % Projection
    x = proj(z1, z2);
    % stop?
    normE = norm(([z1_old; z2_old]-[z1; z2])/lambda(its), 'fro');
    str = sprintf('step %05d: ||ek|| = %.3e', its, normE);
    if its==1 fprintf(str); end
    fprintf('%s%s', char(8*ones(1,length(str))), str); % 8 is tha ascii name of \b
    
    ek(its) = normE;
    
    lambdak = lambdak + lambda(its);
    lambdaK(its) = lambdak;
    ek_bar(its) = norm( ([z10; z20] - [z1; z2]), 'fro') /lambdak;
    
    if normE<tol
        break;
    end
    
    its = its + 1;
    
end
fprintf('\n');

ek = ek(1:its);
ek_bar = ek_bar(1:its);

d0 = norm([z10; z20]-[z1; z2], 'fro');

lambda_ = min(lambda);
ekT = d0 ./ sqrt( (1:its) .* lambda_*(1/alpha-lambda_) );
ekT_bar = 2*d0 ./ lambdaK(1:its);

% EOF
