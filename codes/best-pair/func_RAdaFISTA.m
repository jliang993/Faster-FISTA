function [x,y, its, ek] = func_RAdaFISTA(p,q,r, para, projX, projY)
% The RAdaptive FISTA
% For this scheme, the update of r_k is combined with the restarting FISTA stragety

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
    
    %%% update r_k
    vk = (y_old(:)-x(:))'*(x(:)-x_old(:));
    if vk >= 0
        
            if its>1e4
                theta = 0.99999;
                d = 8;
                p = sqrt(6.0);
            elseif its>1e3
                theta = 0.9999;
                d = 4;
                p = sqrt(4.5);
            elseif its>1e2
                theta = 0.999;
                d = 2;
                p = sqrt(3.0);
            else
                theta = 0.99;
                d = 1;
                p = sqrt(1.5);
            end

        r = r * theta;
            t = t/d;
            y = x;
            
%             p = 1.414;
            q = 1.414;
    end
        
    normE = norm(x-x_old);
    
    if mod(its,10)==0; itsprint(sprintf('        step %08d: norm(ek) = %.3e', its,normE), its); end
    
    ek(its) = normE;
    if normE<=tol; break; end
    
    its = its + 1;
    
end
fprintf('\n');

y = projY(x);

ek = ek(1:its-1);