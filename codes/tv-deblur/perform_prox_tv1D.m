function x = perform_prox_tv1D(y,lambda)
% Exact fast TV-mimization: Chambolle-Darbon algorithm for TV-regularized minimization problem (ROF):
%		min 0.5|| y - x ||^2_2 + lambda || grad(x) ||_1 , this is the anisotropic TV.
%
% Usage
%	x = tvexact(y,lambda)
% Input
%	y           nxn image
%	lambda	    regularization parameter
% Outputs
%	 x	    solution of the problem
%
% Fast TV-minimization by min graph-cuts
%  TV4 = with nearest neighbours interaction
%  TV8 = with also next nearest
%  TV16 = with 16 neighbours
%  
%  Based on the code by 
%      A. Chambolle and J. Darbon: On total variation
%      minimization and surface evolution using parametric maximum flows,
%      preprint (2008).
%  Their code implements Dorit Hochbaum's algorithm:
%     D. S. Hochbaum: An efficient algorithm for image segmentation,
%     Markov random fields and related problems. J. ACM, 48(4):686--701,
%     2001.	
%  
%

x = tvexact_mex(y,lambda);

	

