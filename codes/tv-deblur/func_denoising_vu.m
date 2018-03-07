function [u, its, ek] = func_denoising_vu(f, w)
% ANIsotropic TV
tol = 1e-10;
maxits = 1e5;

beta = 1;
rho = 1.5 /(2*beta);

L = sqrt(8); % ||\nabla||_2
sigma = 1 /( rho + L );
tau = sigma;

alpha = (2*beta*rho) /(4*beta*rho - 1);
lambda = 1/(1.01*alpha);
% lambda = 1;

[m, n] = size(f);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dx = @(x) [ diff(x,1,2), zeros(n,1) ];       % forward difference
dy = @(x) [ diff(x,1,1); zeros(1,m) ];       % discrete y-derivative
dxt = @(x) [ -x(:,1), -diff(x(:,1:m-1),1,2), x(:,m-1) ];      % transpose x-derivative
dyt = @(x) [ -x(1,:); -diff(x(1:n-1,:),1,1); x(n-1,:) ];      % transpose y-derivative
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u = f;

v = zeros(m,n);
vx = dx(v);
vy = dy(v);

ek = zeros(1, maxits);

its = 1;
while(its < maxits)
    
    z_old = [u;vx;vy];
    
    % primal update
    p = u - tau*( dxt(vx)+dyt(vy) + (u - f) );
    
    y = 2*p - u;
    u = u + lambda*(p - u);
    
    % dual update
    tmp_x = vx/sigma + dx(y);
    tmp_y = vy/sigma + dy(y);
    
    qx = sigma*( tmp_x - wthresh(tmp_x, 's', w/sigma) );
    qy = sigma*( tmp_y - wthresh(tmp_y, 's', w/sigma) );
    
    vx = vx + lambda*(qx - vx);
    vy = vy + lambda*(qy - vy);
    
    % % stop?
    normE = norm(z_old-[u;vx;vy], 'fro');
    % stem(itr, normE, '.');
    str = sprintf('step %06d: norm(ek) = %.3e...', its, normE);
    if its==1 fprintf(str); end
    fprintf(char(8*ones(1,length(str)))); % 8 is tha ascii name of \b
    fprintf(str);
    
    if normE<tol
        break;
    end
    
    ek(its) = normE;
    its = its + 1;
end
fprintf('\n');

ek = ek(1:its);