function grad = igrad_logistic(x, i, W, Y)

w = W(i,:);
y = Y(i);

v = exp(-y* (w*x));
% grad = (-y*v) /(1 + v) *(w');
grad = (-y) /(1 + 1/v) *(w');
