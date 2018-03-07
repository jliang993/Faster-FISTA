function [Dx, Dy] = tvMtx(N)

% For SQUARE mtx first
% should be sparse

A = -diag(ones(N,1))+diag(ones(N-1,1),1);
% A(end,:) = 0;

E = eye(N);

Dx = kron(A, E);
Dy = kron(E, A);