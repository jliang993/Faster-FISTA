# Faster-FISTA

Matlab code to reproduce the results of the paper

Faster FISTA

Jingwei Liang, Carola-Bibiane Schönlieb, Preprint 2018


## Quadratic problem

Consider solving the problem below
$$
\min_{x\in \mathbb{R}^n} ~ \frac{1}{2} \|Ax - f\|^2  ,
$$
where $A$ is the Laplacian operator
$$
A = 
\begin{bmatrix}
2 & -1 & & & & & &  \\
-1 & 2 & -1 & & & & &  \\
 & -1 & 2 & -1 & & & &  \\
 & & & \dotsm \\
 &  & & & -1 & 2 & -1 &  \\
 & &  & & & -1 & 2 & -1  \\
 & & &  & & & -1 & 2  \\
\end{bmatrix} .
$$

We set $n = 201$.

 Relative error $\|x_{k}-x_{k-1}\|$          |  Objective function value $\Phi(x_{k}) - \Phi(x^\star)$
:-------------------------:|:-------------------------:
![ ](codes/lse/cmp_lse_ek.png)  |  ![ ](codes/lse/cmp_lse_fk.png)

## Linear inverse problems

Consider solving the problem below
$$
\min_{x\in \mathbb{R}^n} ~ \mu R(x) + \frac{1}{2} \|Ax - f\|^2  .
$$

#### $\ell_{1}$-norm
 Relative error $\|x_{k}-x_{k-1}\|$          |  Objective function value $\Phi(x_{k}) - \Phi(x^\star)$
:-------------------------:|:-------------------------:
![ ](codes/inverse-problem/cmp_fista_ek_lasso.png)  |  ![ ](codes/inverse-problem/cmp_fista_phik_lasso.png)


#### $\ell_{1,2}$-norm
 Relative error $\|x_{k}-x_{k-1}\|$          |  Objective function value $\Phi(x_{k}) - \Phi(x^\star)$
:-------------------------:|:-------------------------:
![ ](codes/inverse-problem/cmp_fista_ek_glasso.png)  |  ![ ](codes/inverse-problem/cmp_fista_phik_glasso.png)


#### $\ell_{\infty}$-norm
 Relative error $\|x_{k}-x_{k-1}\|$          |  Objective function value $\Phi(x_{k}) - \Phi(x^\star)$
:-------------------------:|:-------------------------:
![ ](codes/inverse-problem/cmp_fista_ek_infty.png)  |  ![ ](codes/inverse-problem/cmp_fista_phik_infty.png)


## Total variation based image deblur

##### The codes only run under MacOS

 Original image    |   Blurred image  |   Recovered image          |  Performance comparison
:-------------------------:|:-------------------------:|:-------------------------:|:-------------------------:
![ ](codes/tv-deblur/original-img.png)  |  ![ ](codes/tv-deblur/original-blur.png)  |  ![ ](codes/tv-deblur/original-deblur.png)  |  ![ ](codes/tv-deblur/cmp_fista_tvdeblur.png)


## Principle component pursuit

#### Matrix example

 Mixture matrix    |   Sparse component  |   Low-rank component          |  Performance comparison
:-------------------------:|:-------------------------:|:-------------------------:|:-------------------------:
![ ](codes/pcp/observation.png)  |  ![ ](codes/pcp/sparse-mtx.png)  |  ![ ](codes/pcp/lowrank-mtx.png)  |  ![ ](codes/pcp/cmp_fista_pcp_mtx.png)


#### Video example

 Original frame    |   Foreground  |   Background          |  Performance comparison
:-------------------------:|:-------------------------:|:-------------------------:|:-------------------------:
![ ](codes/pcp/original-frame.png)  |  ![ ](codes/pcp/sparse-component.png)  |  ![ ](codes/pcp/lowrank-component.png)  |  ![ ](codes/pcp/cmp_fista_pcp.png)

Copyright (c) 2018 Jingwei Liang