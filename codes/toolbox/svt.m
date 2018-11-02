function output = svt(mtx, tau)
% economy svd, s is square

[u,s,v] = svd(mtx, 'econ');

s = diag( wthresh(diag(s), 's', tau) );

output = u*s*(v');
