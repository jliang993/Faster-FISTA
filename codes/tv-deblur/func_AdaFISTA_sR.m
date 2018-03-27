function [x, its, ek] = func_AdaFISTA_sR(f,h, nu, p,q,r)
itsprint(sprintf('        step %08d: norm(ek) = %.3e', 1,1), 1);

beta = 1; % convolution kernel
gamma = 1.0*beta;


tau = nu*gamma;

[m,n] = size(f);
F = fft2(f, m,n);
H = fft2(h, m,n);
conjH = conj(H);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gradF = @(x) ifft2( conjH.*( H.*fft2(x) - F) );
proxJ = @(x, t) perform_prox_tv1D(x, t);

FBS = @(x, gamma, tau) proxJ(x-gamma*gradF(x), tau);

% Forward difference, dx for j, dy for i
% dxf = @(x) [diff(x,1,2), zeros(m,1)];
% dyf = @(x) [diff(x,1,1); zeros(1,n)];
% Backward difference
% dxb = @(x) [x(:,1), diff(x,1,2)];
% dyb = @(x) [x(1,:); diff(x,1,1);];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = f;
y = f;

maxits = 1e5;

ek = zeros(1, maxits);

ToL = 1e-10;

t = 1;
first = 1;

its = 1;
while(its<maxits)
    
    x_old = x;
    y_old = y;
    
    x = FBS(y, gamma, tau);
    
    t_old = t;
    t = (p + sqrt(q+r*t_old^2)) /2;
    a = (t_old-1) /t;
    y = x + a*(x-x_old);
    
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
            
            first = 0;
        end
        
        r = r * theta;
        if t >= 4*p/(4 - r)
            t = 4*p/(4 - r)/1;
            y = x;
        end
    end
    
    %%%%%%% stop?
    normE = norm(x_old - x, 'fro');
    if mod(its,1e2)==0
        itsprint(sprintf('        step %08d: norm(ek) = %.3e', its,normE), its);
    end
    
    ek(its) = normE;
    if (normE<ToL)||(normE>1e10); break; end
    
    its = its + 1;
    
end
fprintf('\n');

its = its - 1;
ek = ek(1:its);

% EoF