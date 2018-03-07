# Faster-FISTA

Matlab code to reproduce the results of the paper

Faster FISTA

Jingwei Liang, Carola-Bibiane Schönlieb, Preprint 2018



## Linear inverse problems

Consider solving the problem below
$$
\min_{x\in \mathbb{R}^n} ~ \mu R(x) + \frac{1}{2} \|Ax - f\|^2  .
$$

#### $\ell_{1}$-norm
 Relative error $\|x_{k}-x_{k-1}\|$          |  Objective function value $\Phi(x_{k}) - \Phi(x^\star)$
:-------------------------:|:-------------------------:
![l1-norm](codes/inverse-problem/cmp_fista_ek_lasso.png)  |  ![l1-norm](codes/inverse-problem/cmp_fista_phik_lasso.png)


#### $\ell_{1,2}$-norm
 Relative error $\|x_{k}-x_{k-1}\|$          |  Objective function value $\Phi(x_{k}) - \Phi(x^\star)$
:-------------------------:|:-------------------------:
![l12-norm](codes/inverse-problem/cmp_fista_ek_glasso.png)  |  ![l1-norm](codes/inverse-problem/cmp_fista_phik_glasso.png)


#### $\ell_{\infty}$-norm
 Relative error $\|x_{k}-x_{k-1}\|$          |  Objective function value $\Phi(x_{k}) - \Phi(x^\star)$
:-------------------------:|:-------------------------:
![linfty-norm](codes/inverse-problem/cmp_fista_ek_infty.png)  |  ![l1-norm](codes/inverse-problem/cmp_fista_phik_infty.png)



Copyright (c) 2018 Jingwei Liang