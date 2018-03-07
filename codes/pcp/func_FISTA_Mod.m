function [xs,xl, its, ek] = func_FISTA_Mod(para, GradF, ProxJ, p,q,r)
itsprint(sprintf('        step %08d: norm(ek) = %.3e', 1,1), 1);

% parameter initialization
beta = para.beta;
mu1 = para.mu1;
mu2 = para.mu2;

gamma = 1.0 *beta;
tau = mu1*gamma;
n = para.n;
f = para.f;

% Forward-Backward Step
FBS = @(x, tau) ProxJ(x-gamma*GradF(x), tau);

x0 = zeros(n);

% max number of iterations
maxits = 1e5;
ek = zeros(1, maxits);

its = 1;
ToL = 1e-14;

t = 1;

x = x0;
y = x0;

while(its<maxits)
    
    x_old = x;
    x = FBS(y, tau);
    
    t_old = t;
    t = (p + sqrt(q+r*t_old^2)) /2;
    a = (t_old-1) /t;
    y = x + a*(x-x_old);
    
    %%% stop?
    normE = norm(x(:)-x_old(:), 'fro');
    if mod(its,1e2)==0
        itsprint(sprintf('        step %08d: norm(ek) = %.3e', its,normE), its);
    end
    
    ek(its) = normE;
    if ((normE/prod(n))<ToL)||(normE>1e10) break; end
    
    its = its + 1;
    
end
fprintf('\n');

ek = ek(1:its-1);

xs = x;
xl = svt(f-xs, mu2);