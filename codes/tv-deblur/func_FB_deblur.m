function [x, k,its, ek, dk,dkA, chT] = func_FB_deblur(f,h, nu, h_, acc)
%% set up
tol = 1e-10;
maxits = 1e5;

%%%

beta = 1;
gamma = 1.0*beta;

% alpha = 2*beta / (4*beta - gamma);
% lambda = 1 /(1.01 *alpha); % lambda \in ]0, 2[
lambda = 1;

tau = nu*gamma;

[m,n] = size(f);
F = fft2(f, m,n);
H = fft2(h, m,n);
conjH = conj(H);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gradF = @(x) ifft2( conjH.*( H.*fft2(x) - F) );
proxJ = @(x, t) perform_prox_tv1D(x, t);

T_fbs = @(x, gamma, tau) ...
    (1-lambda)* x + lambda* proxJ(x-gamma*gradF(x), tau);

% Forward difference, dx for j, dy for i
dxf = @(x) [diff(x,1,2), zeros(m,1)];
dyf = @(x) [diff(x,1,1); zeros(1,n)];
% Backward difference
% dxb = @(x) [x(:,1), diff(x,1,2)];
% dyb = @(x) [x(1,:); diff(x,1,1);];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = f;

ek.P = zeros(1, maxits);

K = 100;
flagT = 1;
chT = zeros(1,maxits);

ToL = 5e-13;

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
    gx = dxf(x);
    gy = dyf(x);
    
    Tx = numel(find(abs(gx)>tol));
    Ty = numel(find(abs(gy)>tol));
    
    chT(its) = Tx + Ty;
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

eta = eta_TV(x, gamma, m, h_);
D = 5*nextpow10(ek.P(k));
ek.T = D * eta.^(0:its-k);
%% ||x^k - x^star|| 
xstar = x;

x = f;

dk.P = zeros(1, maxits);

its = 1;
while(its<maxits)
    
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
%% Acceleration
if acc
    
    gammaA = gamma;
    tauA = tau;
    
    x = f;
    dkA.P = zeros(1, maxits);
    
    its = 1;
    while(its<maxits)
        
        if (its==k+1)
            [~, gammaA] = eta_TV(x, gamma, m, h_);
            tauA = nu*gammaA;
        end
        x = T_fbs(x, gammaA, tauA);
        %%%%%%% stop?
        normE = norm(x-xstar, 'fro');
        itsprint(sprintf('step %06d: norm(ek) = %.3e', its,normE), its);
        
        dkA.P(its) = normE;
        if (normE<tol)||(normE>1e10) break; end
        
        its = its + 1;
        
    end
    fprintf('\n');
    
    dkA.P = dkA.P(1:its);
    
    [eta, ~] = eta_TV(x, gammaA, m, h_);
    dkA.T = D * eta.^(0:its-k);
    
else
    
    dkA = 0;
    
end
end
% EoF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [eta, gammaout] = eta_TV(x, gamma, N, h)

[Dx, Dy] = tvMtx(N);
A_ = con2tpzmtx(h, [N,N]);

D = [Dx; Dy];
A = A_;

g = D *x(:);

I = find(abs(g) >= 1e-10);
Ic= setdiff(1:2*N^2, I);

PT = null(D(Ic,:));

AT = A*PT;

s = svd(AT' * AT);

eta = max(abs(1-gamma*s(1)), abs(1-gamma*s(end)));
gammaout = 2.0/ (s(1)+s(end));
end