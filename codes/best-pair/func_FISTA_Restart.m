function [x,y, its, ek] = func_FISTA_Restart(p,q,r, para, projX, projY);

itsprint(sprintf('        step %08d: norm(ek) = %.3e', 1,1), 1);

gamma = .9;

tol = para.tol;
maxits = para.maxits;

ek = zeros(maxits, 1);

x0 = 1e2*randn(2, 1);

x = x0;
y = x;

t = 1;

its = 1;
while(its<maxits)
    
    x_old = x;
    y_old = y;
    
    g = y - projY(y);
    x = projX(y - gamma*g);
    
    t_old = t;
    t = (p + sqrt(q+r*t_old^2)) /2;
    a = (t_old-1) /t;
    y = x + a*(x-x_old);
    
    %%% obective function value criteria
    % if its>3;  if phik(its) > phik(its-1); t = 1; end end
    
    %%% gradient criteria
    % v_check = (GradF(y_old))'*(x-x_old);
    v_check = (y_old-x)'*(x-x_old);
    if v_check > 0; t = 1; end
        
    normE = norm(x-x_old);
    
    if mod(its,10)==0; itsprint(sprintf('        step %08d: norm(ek) = %.3e', its,normE), its); end
    
    ek(its) = normE;
    if normE<=tol; break; end
    
    its = its + 1;
    
end
fprintf('\n');

y = projY(x);

ek = ek(1:its-1);

