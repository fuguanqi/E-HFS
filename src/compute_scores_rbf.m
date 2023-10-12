function xnew = compute_scores_rbf(Data,CandPoint, lambda, gamma, rbf_flag )

% computing the scores of the candidate points from which we select the
% next sample point
% This function is called by the sampling strategies that use a candidate
% point sampling approach (cp.m, rs.m, cptv.m, cptvl.m)
% We compute 2 scores: one derived from the objective function value
% predictions by the RBF, one derived from the distance of the candidates
% to the set of already sampled points
%--------------------------------------------------------------------------
%Author information
%Juliane Mueller
%juliane.mueller2901@gmail.com
%--------------------------------------------------------------------------
%
%Input:
%Data - structure array with all problem information 
%CandPoint - matrix with candidate sample points
%w_r weight for the surrogate model criterion (weight for distance
%criterion is 1-w_r)
%lambda, gamma - RBF parameters
%rbf_flag - indicator of the type of RBF model wanted
%
%Output:
%xnew - new sample point at which the next function evaluation will be done
%--------------------------------------------------------------------------

[mX,nX]=size(CandPoint); %dimensions of the points where function value should be predicted
[mS,nS]=size(Data.S); %dimensions of sample site matrix (already evaluated points)  
R = pdist2(CandPoint,Data.S); %compute pairwise dstances between points in X and S. pdist2 is MATLAB built-in function
                              %pdist2(a*b,c*d)·µ»ØÒ»¸öa*c¾ØÕó
%compute RBF matrix values
if strcmp(rbf_flag,'cub') %cobic RBF
    Phi=R.^3;
elseif strcmp(rbf_flag,'lin') %linear RBF
    Phi=R;
elseif strcmp(rbf_flag,'tps') %thin-plate spline RBF
    Phi=R.^2.*log(R+realmin);
    Phi(logical(eye(size(R)))) = 0;
end
   
p1 = Phi*lambda; %first part of response surface - weighted sum of distances
p2 = [CandPoint,ones(mX,1)]*gamma; % polynomial tail of response surface
yhat=p1+p2; %predicted function value


[~,minID] = min(yhat);
xnew = CandPoint(minID,:);  %new sample point


end %function