function [x, its, ek] = func_deblur_inertial_gfb(f,h, nu, a,b)
%%%%
tol = 1e-8;
maxits = 1e4;

w1 = 1/2;
w2 = 1 - w1;

beta = 1;
gamma = 1.0*beta;

tau = nu*gamma/w1;

[m,n] = size(f);
F = fft2(f, m,n);
H = fft2(h, m,n);
conjH = conj(H);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prox1 = @(x, t) perform_prox_tv1D(x, t); 
prox2 = @(x) max(min(x,255), 0); % box constraint

grad = @(x) ifft2( conjH.*( H.*fft2(x) - F) ); % gradient

proj = @(z1,z2) w1*z1 + w2*z2; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z1 = zeros(size(f));
z2 = z1;

x = proj(z1, z2);

ya1 = z1;
ya2 = z2;

% yb1 = z1;
% yb2 = z2;

xa = w1*ya1 + w2*ya2;
% xb = w1*yb1 + w2*yb2;
xb = x;

ek = zeros(1, maxits);

its = 1;
while(its<maxits)
        
    x_old = x;
    
    z1_old = z1;
    z2_old = z2;
    
    % forward step
    gamma_g = gamma* grad(xb);
    
    % J_A1, SVT
    z1 = ya1 + ( prox1(2*xa-ya1-gamma_g, tau) - xa );
    
    % J_A2, Non-negative constraint
    z2 = ya2 + ( prox2(2*xa-ya2-gamma_g) - xa );
    
    % Projection
    x = proj(z1, z2);
    
    %%%%%%% inertial step
    ya1 = z1 + a*(z1 - z1_old);
    ya2 = z2 + a*(z2 - z2_old);
    
    % yb1 = z1 + b*(z1 - z1_old);
    % yb2 = z2 + b*(z2 - z2_old);
    
    xa = w1*ya1 + w2*ya2;
    xb = x + b*(x - x_old);
    
    % stop?
    normE = norm(([z1_old; z2_old]-[z1; z2]), 'fro');
    str = sprintf('step %05d: ||ek|| = %.3e', its, normE);
    if its==1 fprintf(str); end
    fprintf('%s%s', char(8*ones(1,length(str))), str); % 8 is tha ascii name of \b
    
    ek(its) = normE;
    
    if normE<tol break; end
    
    its = its + 1;
    
end
fprintf('\n');

x = proj(z1, z2);

ek = ek(1:its);

% EOF
