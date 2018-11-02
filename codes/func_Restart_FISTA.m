function [x, its, dk, ek, fk] = func_Restart_FISTA(p,q,r, para, ProxJ,GradF, ObjPhi, xsol)
% Restarting FISTA
verbose = para.verbose;

if verbose
    itsprint(sprintf('        step %08d: residual = %.3e', 1,1), 1);
end

n = para.n;
% J = para.J;
mu = para.mu;
gamma = para.gamma;
tol = para.tol;
maxits = para.maxits  + 1;

tau = mu*gamma;

% Forward--Backward Operator
FBO = @(y) ProxJ(y-gamma*GradF(y), tau);
%% FBS iteration
x0 = para.x0;
x = x0;
y = x0;

dk = zeros(maxits, 1);
ek = zeros(maxits, 1);
fk = zeros(maxits, 1);

% t = t0;
t = 1;
tor = 0;

its = 1;
while(its<maxits)
    
    if verbose
        fk(its) = ObjPhi(x);
        dk(its) = norm(x-xsol, 'fro');
    end
    
    x_old = x;
    y_old = y;
    
    x = FBO(y);
    
    t_old = t;
    t = (p + sqrt(q+r*t_old^2)) /2;
    a = min(1, (t_old-1) /t);
    
    y = x + a*(x-x_old);
    
    %%% gradient criteria
    norm_ = norm(y_old(:)-x(:)) * norm(x(:)-x_old(:));
    vk = (y_old(:)-x(:))'*(x(:)-x_old(:));
    if vk >= tor* norm_
        t = 1;
        y = x; 
    end
    
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
fk = fk(1:its-1);

% EoF