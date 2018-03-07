function g = grad_pcp(x, theta, beta)
% gradient of Moreau Envelope
% fx = M - x
% just for NON-negative Robust PCA

x_ = max( wthresh(x, 's', theta), 0 );
% x_ = max(x_, 0);

g = -(x - x_)/beta;
