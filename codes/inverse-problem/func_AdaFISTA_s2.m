function [x, its, ek, r, Ak,Tk] = func_AdaFISTA_s2(para, GradF,ProxJ,J, p,q,r)
% The MODified FISTA
itsprint(sprintf('        step %08d: norm(ek) = %.3e', 1,1), 1);
w1 = 10;
w2 = 30;
if strcmp(J, 'infty'); w1 = 1e3; w2 = 3e2; end

% set up
beta = para.beta;
mu = para.mu;
n = para.n;
f = para.f;

A = para.A;
alpha = 0; %strong convexity

gamma = 1 *beta;
tau = mu*gamma;

if strcmp(J, 'mc')
    FBS = @(y, gamma, tau) ProxJ(y-gamma*GradF(y, f), mu*gamma, n);
else
    FBS = @(y, gamma, tau) ProxJ(y-gamma*GradF(y, f), mu*gamma);
end
%% FBS iteration
x0 = zeros(prod(n), 1);

x = x0;
y = x0;

tol = 1e-10;
maxits = 1e5;

Ak = zeros(maxits, 1);
Tk = zeros(maxits, 1);
ek = zeros(maxits, 1);

t = 1;

its = 1;
while(its<maxits)
    
    x_old = x;
    x = FBS(y, gamma, tau);
    
    t_old = t;
    t = (p + sqrt(q+r*t_old^2)) /2;
    a = (t_old-1) /t;
    y = x + a*(x-x_old);
    
    if mod(its,w2)==0
        [r, gamma] = alpha_gamma_est(p,gamma, x, A, J);
        % if t>4; t = t/2; end
    end
    
    %%%%%%% stop?
    normE = norm(x_old-x, 'fro');
    
    if mod(its,w1)==0
        itsprint(sprintf('        step %08d: norm(ek) = %.3e', its,normE), its);
    end
    
    Ak(its) = a;
    Tk(its) = t;
    ek(its) = normE;
    % Ek(its) = ObjF(x);
    if (normE<tol)||(normE>1e10); break; end
    
    its = its + 1;
    
end
fprintf('\n');

ek = ek(1:its-1);
Ak = Ak(1:its-1);
Tk = Tk(1:its-1);

% EoF