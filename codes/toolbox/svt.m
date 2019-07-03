function x = svt(mtx, tau)
% economy svd, s is square

[u,s,v] = svd(mtx, 'econ');

s = diag( wthresh(diag(s), 's', tau) );

x = u*s*(v');

x = x(:);