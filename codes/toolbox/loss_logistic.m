function loss = loss_logistic(x, W, Y)


% Logistic Loss
n = length(Y);
loss = 0;

for i=1:n
    w = W(i,:);
    y = Y(i);
    
    v = exp(-y* (w*x));
    loss = loss + log(1+v);
end

loss = loss /n;