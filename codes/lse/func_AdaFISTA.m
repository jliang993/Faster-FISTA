function [x, ek, fk, its] = func_AdaFISTA(p,q,r, para, GradF, ObjF)

itsprint(sprintf('      step %09d: norm(ek) = %.3e', 1,1), 1);

beta = para.beta;
gamma = 1.0* beta;

n = para.n;
tol = para.tol;
maxits = para.maxits;

ek = zeros(1, maxits);
fk = zeros(1, maxits);

x0 = zeros(n, 1);

x = x0;
y = x0;

t = 1e2;
% t = 4*p/(4 - r);

its = 1;
while(its<maxits)
    
    x_old = x;
    x = y - gamma*GradF(y);
    
    t_old = t;
    t = (p + sqrt(q+r*t_old^2)) /2;
    a = (t_old-1) /t;
    y = x + a*(x - x_old);
    
    %%% stop?
    normE = norm(x-x_old, 'fro');
    if mod(its, 1e4)==0
        itsprint(sprintf('      step %09d: norm(ek) = %.3e', its,normE), its);
    end
    
    ek(its) = normE;
    if (normE<tol)||(normE>1e10); break; end
    
    fk(its) = ObjF(x);
    its = its + 1;
    
end
fprintf('\n');

ek = ek(1:its-1);
fk = fk(1:its-1);
