function H = con2tpzmtx(h, dim_)
%CON2TPZMTX convolution to toeplitz matrix
%   h: convolution kernel, size(h) should be ODD number
%   dim_: dimension of the image
%	


n1 = dim_(1);
n2 = dim_(2);
N = n1*n2;

[m1,m2] = size(h);

c1 = (m1+1)/2;
c2 = (m2+1)/2;
pos_i = repmat([1:m1]',1,m2) - c1;
pos_j = repmat([1:m2],m1,1) - c2;

% H = sparse(N, N, 0);
H = zeros(N,N);
cnt = 1;
for j = 0:n2-1
    for i = 0:n1-1
        i_mod = mod(pos_i+i, n1);
        j_mod = mod(pos_j+j, n2);
        
        pos = j_mod *n1 + i_mod + 1;
        
        H(cnt,pos(:)) = h(:);
        
        cnt = cnt + 1;
    end
end