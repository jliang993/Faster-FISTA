function objF = objective_func(J, A,f, L, mu)


%%%%%%
switch J
        %%%%%%%%%%%%%%%%%%%%
    case 'lasso'
        
        R = @(x) norm(x, 1);
        
        %%%%%%%%%%%%%%%%%%%%
    case 'glasso'
        B = 4;
        
        R = @(x) norm(sqrt(sum(reshape(x,[B length(x)/B]).^2, 1)), 1);
        
        %%%%%%%%%%%%%%%%%%%%
    case 'infty'
        
        R = @(x) norm(x, 'inf');
        
        %%%%%%%%%%%%%%%%%%%%
    case 'tv'
        
        R = @(x) norm(diff(x), 1);
        
        %%%%%%%%%%%%%%%%%%%%
    case 'mc'
        
        R = @(x) norm(svd(reshape(x, L)), 1);
        
        %%%%%%%%%%%%%%%%%%%%
end

objF = @(x) mu*R(x) + norm(A*x-f, 2)^2/2;

% EOF