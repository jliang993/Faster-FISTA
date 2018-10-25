function grad = grad_logistic(x, W, Y)


% Logistic Loss
n = length(Y);
grad = zeros(size(W,2), 1);

for i=1:n
    w = W(i,:);
    y = Y(i);
    
    v = exp(-y* (w*x));
    % grad = grad + (-y*v) /(1 + v) *(w');
    grad = grad + (-y) /(1 + 1/v) *(w');
end