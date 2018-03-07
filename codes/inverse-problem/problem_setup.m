function [para, GradF,ProxJ, ObjPhi] = problem_setup(J)

% k-sparsity
k = 8;

oversampling = 6;

%%%%%%
switch J
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'lasso'
        k = 128;
        n = 2048; % length of signal
        m = oversampling*k; % # of measurements
        
        % original x
        xob = 10* randn(n,1);
        xob = sign(xob) .* min( max(abs(xob), 4), 16);
        mask = proj_mask(xob,k/n, 'r');
        xob = xob .* mask;
        
        % GradF = @(x, f) A'*(A*x-f);
        ProxJ = @(x, t) wthresh(x, 's', t);
        
        para.mu = 1;
        
        FuncJ = @(x) para.mu* sum(abs(x));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'glasso'
        k = 128;
        oversampling = 4;
        n = 2048; % length of signal
        m = oversampling*k; % # of measurements
        
        % blk size and # of blks
        B = 8;
        kB= k/B;
        k = B;
        
        % blocksparse = @(B,IB) randn(B,length(IB));
        normB = @(x,B) sqrt(sum(reshape(x,[B length(x)/B]).^2,1));
        % e = @(x,B) x./reshape(repmat(normB(x,B),[B 1]),[L 1]);
        
        % original x
        IB = randperm(n/B);
        IB = IB(1:kB);
        xob = zeros(B,n/B);
        for ib=1:length(IB)
            xob(:,IB(ib)) = SparseVector(B, k, 'GAUSSIAN', true);
        end
        xob = 4* xob(:);
        
        para.B = B;
        
        % GradF = @(x, f) A'*(A*x-f);
        ProxJ = @(x, t) max(1-t./reshape(repmat(normB(x,B),[B 1]),[n 1]),0).*x;
        
        para.mu = 1/2;
        
        FuncJ = @(x) para.mu* sum(normB(x, B));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'infty'
        n = 1024; % length of signal
        m = 1020; % # of measurements
        
        % original x
        xob = randn(n,1);
        
        max_ = max(abs(xob)) + 1;
        
        s = 10;
        Rnd = randperm(n);
        
        xob(Rnd(1:s)) = sign(xob(Rnd(1:s))) .*max_;
        
        % GradF = @(x, f) A'*(A*x-b);
        ProxJ = @(x, t) x - perform_l1ball_projection(x, t);
        
        para.mu = 1;
        
        FuncJ = @(x) para.mu* max(abs(x));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'tv'
        k = 128;
        n = 2048; % length of signal
        m = oversampling*k; % # of measurements
        
        % original x
        xob = 5*SparseVector(n, k, 'GAUSSIAN', true);
        xob = cumsum(xob);
        xob = xob - mean(xob);
        
        % GradF = @(x, f) A'*(A*x-b);
        ProxJ = @(x, t) perform_prox_tv1D(x, t);
        
        para.mu = 1;
        
        FuncJ = @(x) para.mu* sum(abs(diff(x)));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'mc'
        n = [32, 32]; % size of the mtx
        r = 2; % rank of the mtx
        
        m = 3*r*(sum(n) - r);
        % original x
        xob = rand(n(1),r)*diag(r:-1:1)*rand(r,n(2)); % matrix to be tested
        xob = xob(:);
        
        % GradF = @(x, f) A'*(A*x-f);
        ProxJ = @vecsvt;
        
        para.mu = 1;
        
        FuncJ = @(x) para.mu* sum(svd(x));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

tmpA = randn(m, prod(n));
A = tmpA /sqrt(m);

% observation
f0 = A* xob;
sigma = 1e-3 *std(f0);
f = f0 + sigma*randn(m,1);

% gradient
GradF = @(x) A'*(A*x-f);

ObjPhi = @(x) FuncJ(x) + 1/2*norm(A*x(:)-f)^2;

para.A = A;

para.beta = 1 /norm(A)^2;
para.n = n;
para.f = f;
para.xob = xob;

% EOF