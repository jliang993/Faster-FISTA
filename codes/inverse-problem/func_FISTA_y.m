function [x, its, ek, phik] = func_FISTA_y(para, GradF, ProxJ, ObjPhi, J, d)
% FISTA with y
itsprint(sprintf('        step %08d: norm(ek) = %.3e', 1,1), 1);
w = 10;
if strcmp(J, 'infty'); w = 1e3; end

% set up
beta = para.beta;
mu = para.mu;
n = para.n;
% f = para.f;

gamma = 1.0 *beta;
tau = mu*gamma;

if strcmp(J, 'mc')
    FBS = @(y, gamma, tau) ProxJ(y-gamma*GradF(y), tau, n);
else
    FBS = @(y, gamma, tau) ProxJ(y-gamma*GradF(y), tau);
end
%% FBS iteration
x0 = zeros(prod(n), 1);

x = x0;
y = x0;

tol = 1e-10;
maxits = 1e5;

ek = zeros(1, maxits);
phik = zeros(1, maxits);

t = 1;

its = 1;
while(its<maxits)
    
    x_old = x;
    y_old = y;
    
    x = FBS(y, gamma, tau);
    
    t_old = t;
    t = (its+d-1) /d;
    a = (t_old-1) /t;
    y = x + 1.0*a*(x-y_old);
    
    % d = min(d*1.002^its, 75);
    
    %%%%%%% stop?
    normE = norm(x_old-x, 'fro');
    
    if mod(its,w)==0
        itsprint(sprintf('        step %08d: norm(ek) = %.3e', its,normE), its);
    end
    
    ek(its) = normE;
    phik(its) = ObjPhi(x);
    if (normE<tol)||(normE>1e10); break; end
    
    its = its + 1;
    
end
fprintf('\n');

ek = ek(1:its-1);
phik = phik(1:its-1);

% EoF