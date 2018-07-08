function y = Proj2X(x, B1,B2,b1,b2)

ProjB1 = @(x) x - (B1')/ (B1*B1') *(B1*x - b1);

ProjB2 = @(x) x - (B2')/ (B2*B2') *(B2*x - b2);

C1 = [-B1(2)/B1(1), 1];
c1 = - B1(1) + C1(1) + b1;
C2 = [-B2(2)/B2(1), 1];
c2 = - B2(1) + C2(1) + b2;

if (B1*x-b1>0)&&(B2*x-b2>0)
    y = x;
elseif (C1*x-c1<0)&&(C2*x-c2<0)
    y = [1; -B1(1)*1+b1];
elseif x(1)<=1
    y = ProjB1(x);
else
    y = ProjB2(x);
end



% if (B1*x-b1>0)
%     y = x;
% else
%     y = ProjB1(x);
% end