function [x, its, ek, fk] = func_FISTA_Restart(gradF, proxJ, objF, para, p,q,r)
itsprint(sprintf('        step %08d: norm(ek) = %.3e', 1,1), 1);

beta = para.beta; % convolution kernel
mu = para.mu;
gamma = 1.0*beta;

tau = mu*gamma;

f = para.f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FBS = @(x, gamma, tau) proxJ(x-gamma*gradF(x), tau);

% Forward difference, dx for j, dy for i
% dxf = @(x) [diff(x,1,2), zeros(m,1)];
% dyf = @(x) [diff(x,1,1); zeros(1,n)];
% Backward difference
% dxb = @(x) [x(:,1), diff(x,1,2)];
% dyb = @(x) [x(1,:); diff(x,1,1);];

% objF = @(x) mu*( sum(sum(abs(diff(x,1,1)))) + sum(sum(abs(diff(x,1,2)))) )...
%     + 1/2*norm(ifft2( fft2(x) .* fft2(h,m,n) )-f, 'fro')^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = ones(size(f));
y = x;

tol = para.tol;
maxits = para.maxits;

ek = zeros(1, maxits);
fk = zeros(1, maxits);

t = 1;

its = 1;
while(its<maxits)
    
    fk(its) = objF(x);
    
    x_old = x;
    y_old = y;
    
    x = FBS(y, gamma, tau);
    
    t_old = t;
    t = (p + sqrt(q+r*t_old^2)) /2;
    a = (t_old-1) /t;
    y = x + a*(x-x_old);
    
    %%% update r_k
    vk = (y_old(:)-x(:))'*(x(:)-x_old(:));
    if vk >= 0; t = 1; y = x; end
    
    %%%%%%% stop?
    normE = norm(x_old - x, 'fro');
    if mod(its,1e2)==0
        itsprint(sprintf('        step %08d: norm(ek) = %.3e', its,normE), its);
    end
    
    ek(its) = normE;
    if (normE<tol)||(normE>1e10); break; end
    
    its = its + 1;
    
end
fprintf('\n');

its = its - 1;
ek = ek(1:its);
fk = fk(1:its);

% EoF