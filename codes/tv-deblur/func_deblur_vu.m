function [u, its, ek,ekT, ek_bar,ekT_bar] = func_deblur_vu(f, h, w)
% ANIsotropic TV
[m, n] = size(f);

tol = 1e-10;
maxits = 1e5;

F = fft2(f, m,n);
H = fft2(h, m,n);
conjH = conj(H);

L = sqrt(8); % ||\nabla||_2

mu = 1;
beta = mu;

rho = 1.50 /(2*beta);

sigma = 1 /( rho + L );
tau = sigma;

delta = max(1/tau, 1/sigma);

alpha = (2*beta*rho) /(4*beta*rho - 1);

lambda = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dx = @(x) [ diff(x,1,2), zeros(n,1) ];       % forward difference
dy = @(x) [ diff(x,1,1); zeros(1,m) ];       % discrete y-derivative
dxt = @(x) [ -x(:,1), -diff(x(:,1:m-1),1,2), x(:,m-1) ];      % transpose x-derivative
dyt = @(x) [ -x(1,:); -diff(x(1:n-1,:),1,1); x(n-1,:) ];      % transpose y-derivative

Proj = @(x) max(min(x,255), 0); % box constraint
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u = f;

v = f;
vx = dx(v);
vy = dy(v);

z0 = [u;vx;vy];

ek = zeros(1, maxits);
ek_bar = zeros(1,maxits);

lambdak = 0;
lambdaK = zeros(1,maxits);

its = 1;
while(its < maxits)
    
    z_old = [u;vx;vy];
    
    % primal update
    grad = mu*ifft2( conjH.*( H.*fft2(u) - F) );
    p = Proj( u - tau*(dxt(vx)+dyt(vy) + grad) );
    
    y = 2*p - u;
    u = u + lambda*(p - u);
    
    % dual update
    tmp_x = vx/sigma + dx(y);
    tmp_y = vy/sigma + dy(y);
    
    qx = sigma*( tmp_x - wthresh(tmp_x, 's', w/sigma) );
    qy = sigma*( tmp_y - wthresh(tmp_y, 's', w/sigma) );
    
    vx = vx + lambda*(qx - vx);
    vy = vy + lambda*(qy - vy);
    
    z = [u;vx;vy];
    % % stop?
    normE = norm(z_old-z, 'fro');
    % stem(itr, normE, '.');
    str = sprintf('step %06d: norm(ek) = %.3e...', its, normE);
    if its==1 fprintf(str); end
    fprintf(char(8*ones(1,length(str)))); % 8 is tha ascii name of \b
    fprintf(str);
    
    if (normE<tol)||(normE>1e20); break; end
    
    ek(its) = normE;
    lambdak = lambdak + lambda;
    lambdaK(its) = lambdak;
    ek_bar(its) = norm( (z0 - z), 'fro') / lambdak;
    
    its = its + 1;
end
fprintf('\n');

ek = ek(1:its);

d0 = norm(z0-z, 'fro');

ekT = 2*delta*d0/rho ./ sqrt( (1:its) .* lambda*(1/alpha-lambda) );
ekT_bar = 2*delta*d0/rho ./ lambdaK(1:its); % make gap smaller