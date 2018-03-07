function r = alpha_est(p, gamma, x, A, J)

% PT = projTmtx(x, J);
% AT = A*PT;
% 
% AtA = (AT')*AT;

warning off;

n = numel(x);

tol = 1e-12;
switch J
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'lasso'
        I = (abs(x) >= tol);
        AT = A(:, I);
        
    case 'lassoh'
        I = abs(x) >= tol;
        AT = A(:, I);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'glasso'
        I = abs(x) >= tol;
        AT = A(:, I);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'infty'
        I = find(abs(abs(x)-max(abs(x))) <= tol);
        Ic= setdiff(1:n, I);
        sI= zeros(n,1);
        sI(I) = sign(x(I));
        
        PT = eye(n,n);
        PT = double([PT(:,Ic) sI/norm(sI)]);
        
        AT = A *PT;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'tv'
        I = find(abs(grad(x)) >= tol);
        Ic= setdiff(1:n,I);
        
        D = -diag(ones(n,1))+diag(ones(n-1,1),1);
        D(end,:) = 0;
        
        PT = null(D(Ic,:));
        
        AT = A*PT;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

AtA = (AT')*AT;

e = eig(AtA);
e = e(e>1e-4);

alpha = min(e);

r = 4*(1 - sqrt(alpha*gamma))^2/(1-alpha*gamma);




% function r = alpha_est(p, gamma, x, A, J)
% 
% % PT = projTmtx(x, J);
% % AT = A*PT;
% 
% warning off;
% 
% n = numel(x);
% 
% tol = 1e-12;
% switch J
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     case 'lasso'
%         I = (abs(x) >= tol);
%         AT = A(:, I);
%         
%     case 'lassoh'
%         I = abs(x) >= tol;
%         AT = A(:, I);
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     case 'glasso'
%         I = abs(x) >= tol;
%         AT = A(:, I);
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     case 'infty'
%         I = find(abs(abs(x)-max(abs(x))) <= tol);
%         Ic= setdiff(1:n, I);
%         sI= zeros(n,1);
%         sI(I) = sign(x(I));
%         
%         PT = eye(n,n);
%         PT = double([PT(:,Ic) sI/norm(sI)]);
%         
%         AT = A *PT;
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     case 'tv'
%         I = find(abs(grad(x)) >= tol);
%         Ic= setdiff(1:n,I);
%         
%         D = -diag(ones(n,1))+diag(ones(n-1,1),1);
%         D(end,:) = 0;
%         
%         PT = null(D(Ic,:));
%         
%         AT = A*PT;
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end
% 
% AtA = (AT')*AT;
% 
% % AtA = (AtA + (AtA'))/2;
% 
% G = eye(size(AtA)) - gamma*AtA;
% % e = real(eigs(G));
% % eta = min(e(1), 1);
% 
% h = ones(size(G,1), 1);
% for i=1:5e3
%     hnew = G*h;
%     h = hnew/norm(h);
% end
% 
% eta = ((hnew')*h)/((h')*h); 
% eta = min(real(eta), 1);
% % if strcmp(J, 'infty'); eta = (eta+1)/2; end
% 
% asol = (1 - sqrt(1-eta))^2/(eta);
% r = 4*p*(asol - (p-1)/p);