function [x, its, ek, phik, Ak, r,Rk, Vk,Wk] = func_AdaFISTA_sR_Bisection(para, GradF,ProxJ, ObjPhi, J, p,q,r)
% The Adaptive FISTA
% For this scheme, the update of r_k is combined with the restarting FISTA stragety

% itsprint(sprintf('        step %08d: norm(ek) = %.3e', 1,1), 1);
w = 10;
if strcmp(J, 'infty'); w = 1e2; end

% set up
beta = para.beta;
mu = para.mu;
n = para.n;
% f = para.f;

A = para.A;
% alpha = 0; %strong convexity

gamma = 1 *beta;
tau = mu*gamma;

if strcmp(J, 'mc')
    FBS = @(y, gamma, tau) ProxJ(y-gamma*GradF(y), tau, n);
else
    FBS = @(y, gamma, tau) ProxJ(y-gamma*GradF(y), tau);
end
%% FBS iteration
x0 = zeros(prod(n), 1);

x = x0;
y = x0;

tol = 1e-14;
maxits = 1e5;

ek = zeros(1, maxits);
phik = zeros(1, maxits);

theta = 0.99;

Rk = zeros(maxits, 1);
cR = 1;

Ak = zeros(maxits, 1);
Vk = zeros(maxits, 1);
Wk = zeros(maxits, 1);

r_left = 1;
r_right = 4;

first = 1;
f_left = 1;
% gap = 1;
% k_old = 1;

t = 1;

its = 1;
while(its<maxits)
    
    x_old = x;
    y_old = y;
    
    x = FBS(y, gamma, tau);
    
    t_old = t;
    t = (p + sqrt(q+r*t_old^2)) /2;
    a = (t_old-1) /t;
    y = x + a*(x-x_old);
    
    Ak(its) = a;
    
    %%% update r_k
    vk = (y_old-x)'*(x-x_old);
    Vk(its) = vk;
    
    if first % fist oscillation
        
        if vk >= 0
            r = (r_left + r_right)/2;
            
            t = 4*p/(4 - r);
            % if t >= 4*p/(4 - r); t = 4*p/(4 - r); y = x; end
            
            a = (p-1)/p + r/(4*p);
            y = x + a*(x-x_old);
            
            kpos = its;
            gap = its;
            
            if gap>1e2; gap = 1e2; end
            
            [1/2, its/1e2, r, r_left, r_right]
            
            first = 0;
        end
        
        % f_right = 1;
        
    else
        
        if (vk <= 0)&&(its-kpos>=gap)
            
            sum_k = sum( diff( sign( Vk(kpos+1:its) ) ) );
            
            if sum_k==0
                r_left = r;
                r = (r_left + r_right)/2;
                
                t = 4*p/(4 - r);
                % if t >= 4*p/(4 - r); t = 4*p/(4 - r)/1.1; y = x; end
                
                a = (p-1)/p + r/(4*p);
                y = x + a*(x-x_old);
                
                kpos = its;
                
                [1, its/1e2, r, r_left, r_right]
                
                f_left = 1;
            end
            
        end
        
        if (vk >= 0)&&(f_left)
            r_right = r;
            r = (r_left + r_right)/2;
            
            t = 4*p/(4 - r);
            
            a = (p-1)/p + r/(4*p);
            y = x + a*(x-x_old);
            
            % if t >= 4*p/(4 - r); t = 4*p/(4 - r)/2; y = x; end
            
            kpos = its;
            
            [2, its/1e2, r, r_left, r_right]
            
            f_left = 0;
        end
        
    end
    
    
    %%%%%%% stop?
    normE = norm(x_old-x, 'fro');
    
    % if mod(its,w)==0; itsprint(sprintf('        step %08d: norm(ek) = %.3e', its,normE), its); end
    
    ek(its) = normE;
    phik(its) = ObjPhi(x);
    if (normE<tol)||(normE>1e10); break; end
    
    its = its + 1;
    
end
fprintf('\n');

% kpos

ek = ek(1:its-1);
phik = phik(1:its-1);

Ak = Ak(1:its-1);
Vk = Vk(1:its-1);

% EoF