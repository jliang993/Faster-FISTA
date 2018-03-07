function [y, t, a] = func_FISTA_step(p,k, t_old, x,x_old)

t = (k + p -1) /p;

a = (t_old - 1) /t;

y = x + a*(x - x_old);
