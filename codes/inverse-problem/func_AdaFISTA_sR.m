function [x, its, ek, phik, r,Rk, Vk,Wk] = func_AdaFISTA_sR(para, GradF,ProxJ, ObjPhi, J, p,q,r)
% The Adaptive FISTA
% For this scheme, the update of r_k is combined with the restarting FISTA stragety

itsprint(sprintf('        step %08d: norm(ek) = %.3e', 1,1), 1);
w = 10;
if strcmp(J, 'infty'); w = 5e2; end

% set up
beta = para.beta;
mu = para.mu;
n = para.n;
% f = para.f;

A = para.A;
% alpha = 0; %strong convexity

gamma = 1 *beta;
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

ek = zeros(maxits, 1);
phik = zeros(maxits, 1);

theta = 0.99;

Rk = zeros(maxits, 1);
cR = 1;

first = 1;

Vk = zeros(maxits, 1);
Wk = zeros(maxits, 1);

t = 1;

its = 1;
while(its<maxits)
    
    x_old = x;
    y_old = y;
    
    x = FBS(y, gamma, tau);
    
    t_old = t;
    t = (p + sqrt(q+r*t_old^2)) /2;
    a = (t_old-1) /t;
    y = x + a*(x-x_old);
    
    %%% update r_k
    vk = (y_old(:)-x(:))'*(x(:)-x_old(:));
    Vk(its) = vk;
    if vk >= 0
        
        if first
            if its>1e4
                theta = 0.999999;
            elseif its>1e3
                theta = 0.99999;
            elseif its>5e2
                theta = 0.9995;
            elseif its>1e2
                theta = 0.995;
            elseif its>50
                theta = 0.985;
            else
                theta = 0.95;
            end
            
            first = 0;
        end
        
        Wk(cR) = its;
        Rk(cR) = r; cR = cR + 1;
        r = r * theta;
        if t >= 4*p/(4 - r)
            t = 4*p/(4 - r)/1;
            y = x;
        end
    end
    
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

Vk = Vk(1:its-1);

Rk = Rk(1:cR-1);
Wk = Wk(1:cR-1);

% EoF