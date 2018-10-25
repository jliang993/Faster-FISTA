% Total variation denoising of 1-D signals, a.k.a. Fused lasso
% signal approximator, by Laurent Condat.
%
% Version 2.0, Aug. 3, 2017.
%
% Given a real vector y of length N and a real lambda>=0, the
% goal is to compute the real vector x minimizing
%    ||x-y||_2^2/2 + lambda.TV(x),
% where ||x-y||_2^2 = sum_{n=1}^{N} (x[n]-y[n])^2 and
% TV(x) = sum_{n=1}^{N-1} |x[n+1]-x[n]|.
%
% I proposed a fast and exact algorithm (say, the version 1.0)
% for this problem in L. Condat, "A direct algorithm for 1D
% total variation denoising," IEEE Signal Proc. Letters, vol.
% 20, no. 11, pp. 1054-1057, Nov. 2013.
% It has worst case complexity O(N^2) but it is recognized as
% the fastest in practice (using the C code on my webpage).
%
% The present code is a MATLAB implementation of a NEW
% algorithm, which combines the advantages of the v1 algorithm
% with the optimal O(N) complexity of the taut string algorithm.
% That is, it is exact, numerically robust (averages of values
% computed by Welford-Knuth running mean algorithm, not by sum
% divided by length), roughly as fast as the v1 algorithm, and
% it has linear time complexity.
%
% In a nutshell, the algorithm is based on the classical Pool
% Adjacent Violators Algorithm for isotonic regression, to
% maintain two nonincreasing and nondecreasing (instead of
% constant in the v1) lower and upper approximations of the
% signal.
%
% C code available on my webpage.
%
% If lambda=0, the algorithm returns x=y, but not exactly, only
% up to machine precision (an operation like c+(y[n]-c) does not
% always return y[n] exactly).


function x = prox_tv(y, lambda)
N = length(y);
if N<=1, x=y; return; end;
x = zeros(size(y)); % y can be a row or column vector.
indstart_low=zeros(1,N); % starting indices of constant
% segments of the lower approximation x_low
indstart_up=zeros(1,N); % starting indices of constant
% segments of the upper approximation x_up
j_low = 1; % index to count the segments of x_low
j_up = 1;  % same for x_up
jseg = 1;  % segment number of the current part under
% construction
indjseg = 1; % starting index of the current part
% we have indjseg = indstart_low(jseg) = indstart_up(jseg)
indstart_low(1) = 1; % starting index of the j_low-th
% segment of x_low
indstart_up(1) = 1; % same for x_up
x_low_first = y(1)-lambda; % value of the first segment
% of the part of x_low under construction
x_up_first = y(1)+lambda; % same for x_up
x_low_curr = x_low_first; % value of the last segment
% of the part of x_low under construction
x_up_curr = x_up_first; % same for x_up
% the constant value of x_low over the j-th segment is stored
% in x(indstart_low(j_low)), except for j=jseg, where the
% value is x_low_first. Same for x_up. Indeed, the parts of
% x_low and x_up under construction have distinct jump
% locations, but same starting index jseg.
for i = 2:N-1
    if y(i)>=x_low_curr
        if y(i)<=x_up_curr
            % fusion of x_up to keep it nondecreasing
            x_up_curr=x_up_curr+(y(i)-x_up_curr)/(i-indstart_up(j_up)+1);
            x(indjseg)=x_up_first;
            while (j_up>jseg)&&(x_up_curr<=x(indstart_up(j_up-1)))
                j_up=j_up-1;
                x_up_curr=x(indstart_up(j_up))+(x_up_curr-x(indstart_up(j_up)))*...
                    ((i-indstart_up(j_up+1)+1)/(i-indstart_up(j_up)+1));
            end
            if j_up==jseg,  % a jump in x downwards is possible
                % the fusion of x_low has not been done yet, but this is OK.
                while (x_up_curr<=x_low_first)&&(jseg<j_low)
                    % the second test should always be true if the first one
                    % is true and lambda>0, but this is a numerical safeguard.
                    % And it is necessary if lambda=0.
                    % validation of segments of x_low in x
                    jseg=jseg+1;
                    x(indjseg:indstart_low(jseg)-1)=x_low_first;
                    x_up_curr=x_up_curr+(x_up_curr-x_low_first)*...
                        ((indstart_low(jseg)-indjseg)/(i-indstart_low(jseg)+1));
                    indjseg=indstart_low(jseg);
                    x_low_first=x(indjseg);
                end
                x_up_first=x_up_curr;
                j_up=jseg;
                indstart_up(jseg)=indjseg;
            else, x(indstart_up(j_up))=x_up_curr; end
        else % we start a new segment in x_up
            j_up=j_up+1;
            indstart_up(j_up)=i;
            x(i)=y(i);
            x_up_curr=x(i);
        end
        % fusion of x_low to keep it nonincreasing
        x_low_curr=x_low_curr+(y(i)-x_low_curr)/(i-indstart_low(j_low)+1);
        x(indjseg)=x_low_first;
        while (j_low>jseg)&&(x_low_curr>=x(indstart_low(j_low-1)))
            j_low=j_low-1;
            x_low_curr=x(indstart_low(j_low))+(x_low_curr-x(indstart_low(j_low)))*...
                ((i-indstart_low(j_low+1)+1)/(i-indstart_low(j_low)+1));
        end
        if j_low==jseg  % a jump in x upwards is possible
            while (x_low_curr>=x_up_first)&&(jseg<j_up)
                % validation of segments of x_up in x
                jseg=jseg+1;
                x(indjseg:indstart_up(jseg)-1)=x_up_first;
                x_low_curr=x_low_curr+(x_low_curr-x_up_first)*...
                    ((indstart_up(jseg)-indjseg)/(i-indstart_up(jseg)+1));
                indjseg=indstart_up(jseg);
                x_up_first=x(indjseg);
            end
            x_low_first=x_low_curr;
            j_low=jseg;
            indstart_low(jseg)=indjseg;
            if indjseg==i, % this part is not mandatory, it is a kind
                % of reset to increase numerical robustness.
                % If we are here, this is just after a jump upwards has
                % been validated. We have x_up_first=y(i).
                x_low_first=x_up_first-2*lambda;
            end;
        else, x(indstart_low(j_low))=x_low_curr; end
    else
        % we start a new segment in x_low
        j_low = j_low+1;
        indstart_low(j_low) = i;
        x(i)=y(i);
        x_low_curr=x(i);
        % fusion of x_up to keep it nondecreasing
        x_up_curr=x_up_curr+(y(i)-x_up_curr)/(i-indstart_up(j_up)+1);
        x(indjseg)=x_up_first;
        while (j_up>jseg)&&(x_up_curr<=x(indstart_up(j_up-1)))
            j_up=j_up-1;
            x_up_curr=x(indstart_up(j_up))+(x_up_curr-x(indstart_up(j_up)))*...
                ((i-indstart_up(j_up+1)+1)/(i-indstart_up(j_up)+1));
        end
        if j_up==jseg  % a jump in x downwards is possible
            while (x_up_curr<=x_low_first)&&(jseg<j_low)
                % validation of segments of x_low in x
                jseg=jseg+1;
                x(indjseg:indstart_low(jseg)-1)=x_low_first;
                x_up_curr=x_up_curr+(x_up_curr-x_low_first)*...
                    ((indstart_low(jseg)-indjseg)/(i-indstart_low(jseg)+1));
                indjseg=indstart_low(jseg);
                x_low_first=x(indjseg);
            end
            x_up_first=x_up_curr;
            j_up=jseg;
            indstart_up(jseg)=indjseg;
            if indjseg==i, % this part is not mandatory, it is a kind
                % of reset to increase numerical robustness.
                x_up_first=x_low_first+2*lambda;
            end;
        else, x(indstart_up(j_up))=x_up_curr; end
    end
end
i=N;
if y(i)+lambda<=x_low_curr
    % the segments of x_low are validated
    while jseg<j_low
        jseg=jseg+1;
        x(indjseg:indstart_low(jseg)-1) = x_low_first;
        indjseg=indstart_low(jseg);
        x_low_first=x(indjseg);
    end
    x(indjseg:i-1) = x_low_first;
    x(i)=y(i)+lambda;
elseif y(i)-lambda>=x_up_curr
    % the segments of x_up are validated
    while jseg<j_up
        jseg=jseg+1;
        x(indjseg:indstart_up(jseg)-1) = x_up_first;
        indjseg=indstart_up(jseg);
        x_up_first=x(indjseg);
    end
    x(indjseg:i-1) = x_up_first;
    x(i)=y(i)-lambda;
else
    % fusion of x_low to keep it nonincreasing
    x_low_curr=x_low_curr+(y(i)+lambda-x_low_curr)/(i-indstart_low(j_low)+1);
    x(indjseg)=x_low_first;
    while (j_low>jseg)&&(x_low_curr>=x(indstart_low(j_low-1)))
        j_low=j_low-1;
        x_low_curr=x(indstart_low(j_low))+(x_low_curr-x(indstart_low(j_low)))*...
            ((i-indstart_low(j_low+1)+1)/(i-indstart_low(j_low)+1));
    end
    if j_low==jseg % the segments of x_up must be validated
        if x_up_first>=x_low_curr % same unique segment of x_low and x_up
            x(indjseg:i)=x_low_curr;
        else
            % fusion of x_up to keep it nondereasing
            x_up_curr=x_up_curr+(y(i)-lambda-x_up_curr)/(i-indstart_up(j_up)+1);
            x(indjseg)=x_up_first;
            while (j_up>jseg)&&(x_up_curr<=x(indstart_up(j_up-1)))
                j_up=j_up-1;
                x_up_curr=x(indstart_up(j_up))+(x_up_curr-x(indstart_up(j_up)))*...
                    ((i-indstart_up(j_up+1)+1)/(i-indstart_up(j_up)+1));
            end
            x(indstart_up(j_up):i)=x_up_curr;
            while jseg<j_up  % the segments of x_up are validated
                jseg=jseg+1;
                x(indjseg:indstart_up(jseg)-1) = x_up_first;
                indjseg=indstart_up(jseg);
                x_up_first=x(indjseg);
            end
        end
    else 	% the segments of x_low must be validated
        x(indstart_low(j_low):i)=x_low_curr;
        while jseg<j_low
            jseg=jseg+1;
            x(indjseg:indstart_low(jseg)-1) = x_low_first;
            indjseg=indstart_low(jseg);
            x_low_first=x(indjseg);
        end
    end
end
end