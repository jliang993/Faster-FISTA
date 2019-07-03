function [x, its, dk, ek] = func_FB(para, ProxJ,GradF, xsol)
% MODified FISTA
verbose = para.verbose;

if verbose
    itsprint(sprintf('        step %08d: residual = %.3e', 1,1), 1);
end

n = para.n;
mu = para.mu;
gamma = para.gamma;
tol = para.tol;
maxits = para.maxits  + 1;

tau = mu*gamma;

% Forward--Backward Operator
FBO = @(x) ProxJ(x-gamma*GradF(x), tau ./(1e-4+abs(x)));
% FBO = @(x) ProxJ(x-gamma*GradF(x), taus);
%% FBS iteration
x0 = ones(n, 1);
x = x0;

dk = zeros(maxits, 1);
ek = zeros(maxits, 1);

its = 1;
while(its<maxits)
    
    if verbose
        dk(its) = norm(x-xsol, 'fro');
    end
    
    x_old = x;
    
    x = FBO(x);
    
    %%%%%%% stop?
    res = norm(x_old-x, 'fro');
    
    if verbose&&mod(its, 1e2)==0
        itsprint(sprintf('        step %08d: residual = %.3e', its,res), its);
    end
    
    ek(its) = res;
    if (res/prod(n)<tol)||(res>1e10); break; end
    
    its = its + 1;
    
end
fprintf('\n');

dk = dk(1:its-1);
ek = ek(1:its-1);

% EoF