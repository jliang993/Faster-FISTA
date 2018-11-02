function [x, its, dk, ek, fk] = func_Greedy_FISTA(para, ProxJ,GradF, ObjPhi, xsol)
% Greddy FISTA
verbose = para.verbose;

if verbose
    itsprint(sprintf('        step %08d: residual = %.3e', 1,1), 1);
end

n = para.n;
% J = para.J;
mu = para.mu;
gamma0 = para.gamma;
gamma = para.c_gamma* gamma0;
tol = para.tol;
maxits = para.maxits  + 1;

% tau = mu*gamma;

% Forward--Backward Operator
FBO = @(y, gamma) ProxJ(y-gamma*GradF(y), mu*gamma);
%% FBS iteration
x0 = para.x0;
x = x0;
y = x0;

dk = zeros(maxits, 1);
ek = zeros(maxits, 1);
fk = zeros(maxits, 1);

a = para.a;

tor = 0;
S = 1;
xi = 0.96;
% first = 1;
% e0 = 1e5;

its = 1;
while(its<maxits)
    
    if verbose
        fk(its) = ObjPhi(x);
        dk(its) = norm(x-xsol, 'fro');
    end
    
    x_old = x;
    y_old = y;
    
    x = FBO(y, gamma);
    
    y = x + a(its)*(x-x_old);
    
    %%% gradient criteria
    norm_ = norm(y_old(:)-x(:)) * norm(x(:)-x_old(:));
    vk = (y_old(:)-x(:))'*(x(:)-x_old(:));
    if vk >= tor* norm_
        y = x;
        % if first; e0 = max(ek(1:its)); first = 0; end
    end
    
    %%%%%%% stop?
    res = norm(x_old-x, 'fro');
    
    if verbose&&mod(its, 1e2)==0
        itsprint(sprintf('        step %08d: residual = %.3e', its,res), its);
    end
    
    ek(its) = res;
    if (res/prod(n)<tol)||(res>1e10); break; end
    
    %%% safeguard
     if res>S*ek(1); gamma = max(gamma0, gamma*xi); end % x = x_old; y = x_old; end
    
    its = its + 1;
    
end
fprintf('\n');

% if verbose; fprintf('\n'); disp(gamma/para.gamma); fprintf('\n'); end

dk = dk(1:its-1);
ek = ek(1:its-1);
fk = fk(1:its-1);

% EoF