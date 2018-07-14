function [x, ek, fk, its] = func_RAdaFISTA(x0, p,q,r, para, GradF, ObjF)

itsprint(sprintf('      step %09d: norm(ek) = %.3e', 1,1), 1);

beta = para.beta;
gamma = 1.0* beta;

n = para.n;
tol = para.tol;
maxits = para.maxits;

ek = zeros(maxits, 1);
fk = zeros(maxits, 1);

% x0 = zeros(n, 1);
% x0 = 1e4*ones(n, 1);

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
%             if its>1e4
%                 theta = 0.999999;
%             elseif its>1e3
%                 theta = 0.99999;
%             elseif its>5e2
%                 theta = 0.9995;
%             elseif its>1e2
%                 theta = 0.995;
%             elseif its>50
%                 theta = 0.985;
%             else
%                 theta = 0.95;
%             end
            kpos = its;
            first = 0;
        end
        
%         r = r * theta;
%         if t >= 4*p/(4 - r)
%             t = 4*p/(4 - r)/1;
%             y = x;
%         end

            if its>1e4
                theta = 0.99999;
                d = 8;
            elseif its>1e3
                theta = 0.9999;
                d = 4;
            elseif its>1e2
                theta = 0.999;
                d = 2;
            else
                theta = 0.99;
                d = 1;
            end

        
        r = r * theta;
            t = t/d;
            y = x;
            
            % p = 1.618;
            p = 1.618;
            q = 1.0;
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
