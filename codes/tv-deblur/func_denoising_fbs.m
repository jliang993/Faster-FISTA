function [x, its, ek]...
    = tv_denoising_fbs(f, w)
%% set up
tol = 1e-10;
maxits = 1e5;

beta = 1;
gamma = 1.5 *beta;

% alpha = 2*beta / (4*beta - gamma);
% lambda = 1 /(1.01 *alpha);
lambda = 1.0;

tau = w*gamma;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
grad = @(x) (x - f);
prox = @(x, t) perform_prox_tv1D(x, t);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = f;

ek = zeros(1, maxits);

its = 1;
while(its<maxits)
    x_old = x;
    
    g = grad(x);
    u = prox(x-gamma*g, tau);
    x = x + lambda*( u - x );
    
    %%%%%%% stop?
    normE = norm(x_old - x, 'fro');
    str = sprintf('step %06d: norm(ek) = %.3e', its, normE);
    if its==1 fprintf(str); end
    fprintf('%s%s', char(8*ones(1,length(str))), str); 
    
    ek(its) = normE;
    
    if (normE<tol)||(normE>1e10)
        break;
    end
    
    its = its + 1;
end
fprintf('\n');

ek = ek(1:its);

%%%% EOF