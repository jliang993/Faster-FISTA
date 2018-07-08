function y = Proj2Y(x, D1,D2,d1,d2)

ProjD1 = @(x) x - (D1')/ (D1*D1') *(D1*x - d1);

ProjD2 = @(x) x - (D2')/ (D2*D2') *(D2*x - d2);

E1 = [-D1(2)/D1(1), 1];
e1 = - D1(1) + E1(1)*0.8 + d1;
E2 = [-D2(2)/D2(1), 1];
e2 = - D2(1) + E2(1)*0.8 + d2;

if (D1*x-d1<0)&&(D2*x-d2<0)
    y = x;
elseif (E1*x-e1>0)&&(E2*x-e2>0)
    y = [0.8; -D1(1)*0.8+d1];
elseif x(1)<=0.8 % (E1*x-e1<0)&&(E2*x-e2>0)
    y = ProjD1(x);
else
    y = ProjD2(x);
end
