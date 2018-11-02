function mask = proj_mask(u, r, type)
%PROJ_MASK generates random projection mask
%   input arguments:
%               u: input object
%               r: data KNOWN ratio
%       type: data missing type
%                   'r': random missing rows
%                   'c': random missing columns
%                   'p': random missing pixels

% For ROW COLUMN missing cases, the max gap between two known
% rows/columns is at most 2(1/r-1) with probability r^2

% code duplicate for cases 'r' & 'c'

[m,n, Ch] = size(u);

mask = zeros(m,n, Ch);

switch type
    
    case 'r'
        
        gap = 1/r;
        gap_ = floor(gap);
        K = floor(m*r);
        
        for ch = 1:Ch
            for i=1:K
                j = floor((i-1)*gap) + randperm(gap_,1);
                mask(j,:, ch) = 1;
            end
        end
        
    case 'c'
        
        gap = 1/r;
        gap_ = floor(gap);
        K = floor(n*r);
        
        for ch = 1:Ch
            for i=1:K
                j = floor((i-1)*gap) + randperm(gap_,1);
                mask(:,j, ch) = 1;
            end
        end
        
    case 'p'
        
        pix = randperm(m*n*Ch);
        r = fix(r*m*n*Ch);
        pix = pix(1:r);
        mask(pix) = 1;
        
    otherwise % pixel-wise missing
        
        pix = randperm(m*n*Ch);
        r = fix(r*m*n*Ch);
        pix = pix(1:r);
        mask(pix) = 1;
        
end
