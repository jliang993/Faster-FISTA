function [x, ek, fk, its] = func_GD(x0, para, GradF, ObjF)

itsprint(sprintf('      step %09d: norm(ek) = %.3e', 1,1), 1);

beta = para.beta;
gamma = 1.5* beta;

n = para.n;
tol = para.tol;
maxits = para.maxits;

ek = zeros(maxits, 1);
fk = zeros(maxits, 1);

% x0 = zeros(n, 1);
% x0 = 1e4*ones(n, 1);

x = x0;


its = 1;
while(its<maxits)
    
    fk(its) = ObjF(x);
    
    x_old = x;
    x = x - gamma*GradF(x);
    
%     if mod(its, 10)==0
%         plot(x(1), x(2), 'b.', 'markersize', 12);
%     end
    
    %%% stop?
    normE = norm(x-x_old, 'fro');
    if mod(its, 1e4)==0
        itsprint(sprintf('      step %09d: norm(ek) = %.3e', its,normE), its);
    end
    
    ek(its) = normE;
    if (normE<tol)||(normE>1e10); break; end
    
    its = its + 1;
    
end
fprintf('\n');

ek = ek(1:its-1);
fk = fk(1:its-1);