function H = con2tpzmtx1D(h, N)
%CON2TPZMTX convolution to toeplitz matrix
%   h: convolution kernel, size(h) should be ODD number
%   dim_: dimension of the image
%

reversh = h(end:-1:1);

L = length(h);

P = (-L:-1) + 1;

H = zeros(N,N);
for i = 1:N
    
    idx = pmod(N, P+i);
    
    H(i, idx) = reversh;
    
end

end
%%%%%%%%%%%%%%%%%%%%%%
function idxout = pmod(N, idxin)

idxout = N* (1/2 - sign(idxin)/2) + idxin;

end