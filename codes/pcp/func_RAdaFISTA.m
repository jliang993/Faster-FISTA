function [xs,xl, its, ek, fk] = func_RAdaFISTA(para, GradF, ProxJ, p,q,r, objF)
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
ek = zeros(maxits, 1);
fk = zeros(maxits, 1);

its = 1;
ToL = 1e-14;

t = 1;
first = 1;

x = x0;
y = x0;

while(its<maxits)
    
    fk(its) = objF(svt(f-x, mu2), x);
    
    x_old = x;
    y_old = y;
    
    x = FBS(y, tau);
    
    t_old = t;
    t = (p + sqrt(q+r*t_old^2)) /2;
    a = (t_old-1) /t;
    y = x + a*(x-x_old);
    
    
    %%% update r_k
    vk = (y_old(:)-x(:))'*(x(:)-x_old(:));
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
        
        r = r * theta;
        if t >= 4*p/(4 - r)
            t = 4*p/(4 - r)/1;
            y = x;
        end
    end
    
    %%% stop?
    normE = norm(x(:)-x_old(:), 'fro');
    if mod(its,1e2)==0
        itsprint(sprintf('        step %08d: norm(ek) = %.3e', its,normE), its);
    end
    
    ek(its) = normE;
    if ((normE/prod(n))<ToL)||(normE>1e10); break; end
    
    its = its + 1;
    
end
fprintf('\n');

ek = ek(1:its-1);
fk = fk(1:its-1);

xs = x;
xl = svt(f-xs, mu2);