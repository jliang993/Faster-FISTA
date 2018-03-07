function [x, k,its, ek, dk, chT] = func_FB_SpikeDecon(f,A, nu)
%% set up
tol = 1e-10;
maxits = 1e5;

beta = 1/norm(A)^2;
gamma = 1.5*beta;

tau = nu*gamma;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gradF = @(x) A'*(A*x - f);
proxJ = @(x, t) wthresh(x, 's', t);

T_fbs = @(x, gamma, tau) proxJ(x-gamma*gradF(x), tau);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = f;

ek.P = zeros(1, maxits);

K = 1e3;
flagT = 1;
chT = zeros(1,maxits);

ToL = 1e-15;

its = 1;
while(its<maxits)
    
    x_old = x;
    x = T_fbs(x, gamma, tau);
    
    %%%%%%% stop?
    normE = norm(x_old - x, 'fro');
    itsprint(sprintf('step %06d: norm(ek) = %.3e', its,normE), its);
    
    ek.P(its) = normE;
    if (normE<ToL)||(normE>1e10) break; end
    
    %%%%%%% T check
    chT(its) = numel(find(abs(x)>tol));  
    if (its>K)&&flagT
        cht = sum(abs(diff(chT(its-K:its))));
        if cht==0
            flagT = 0;
            k = its - K + 5;
        end
    end
    
    its = its + 1;
    
end
fprintf('\n');

its = its - 1;
ek.P = ek.P(1:its);
chT = chT(1:its);

eta = eta_L1(x, gamma, A);
D = 5*nextpow10(ek.P(k));
ek.T = D * eta.^(0:its-k);
%% ||x^k - x^star|| 
xstar = x;

x = f;

dk.P = zeros(1, maxits);

its = 1;
while(its<maxits)
    
    % x_old = x;
    x = T_fbs(x, gamma, tau);
    
    %%%%%%% stop?
    normE = norm(x - xstar, 'fro');
    itsprint(sprintf('step %06d: norm(ek) = %.3e', its,normE), its);
    
    dk.P(its) = normE;
    if (normE<tol)||(normE>1e10) break; end
    
    its = its + 1;
    
end
fprintf('\n');

its = its - 1;
dk.P = dk.P(1:its);

D = 5*nextpow10(dk.P(k));
dk.T = D * eta.^(0:its-k);

end
% EoF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function eta = eta_L1(x, gamma, A)

tol = 1e-8;

I = abs(x) >= tol;
AT = A(:, I);

s = svd(AT' * AT);

eta = max(abs(1-gamma*s(1)), abs(1-gamma*s(end)));
end