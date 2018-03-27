function [x, ek, fk, its] = func_AdaFISTA_sR(p,q,r, para, GradF, ObjF)

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

t = 1;
first = 1;
% t = 4*p/(4 - r);

its = 1;
while(its<maxits)
    
    x_old = x;
    y_old = y;
    
    x = y - gamma*GradF(y);
    
    t_old = t;
    t = (p + sqrt(q+r*t_old^2)) /2;
    a = (t_old-1) /t;
    y = x + a*(x - x_old);
    
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
            kpos = its;
            first = 0;
        end
        
        r = r * theta;
        if t >= 4*p/(4 - r)
            t = 4*p/(4 - r)/1;
            y = x;
        end
    end
    
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

r
kpos

ek = ek(1:its-1);
fk = fk(1:its-1);
