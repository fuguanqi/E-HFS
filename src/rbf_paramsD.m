function [lambda, gamma] = rbf_params(Data, rbf_flag)

% computing the radial basis function parameters
% This function is called by all sampling strategies whenever a new data
% point (x, f(x)) has bee obtained.
%--------------------------------------------------------------------------
%Author information
%Juliane Mueller
%juliane.mueller2901@gmail.com
%--------------------------------------------------------------------------
%
%Input:
%Data - structure array with all problem information 
%rbf_flag - indicator of the type of RBF model wanted
%
%Output:
%lambda - weight parameters in the sum of RBF values (first term)
%gamma - parameters of polynomial tail
%--------------------------------------------------------------------------


% distances=mydist(Data.S,Data);
dist1=pdist2(Data.S(:,Data.category),Data.S(:,Data.category),"hamming");
dist1=dist1.*numel(Data.category);
non_cate=1:Data.dim;
non_cate(Data.category)=[];
dist2=pdist2(Data.S(:,non_cate),Data.S(:,non_cate));
distances=sqrt(dist2.^2+dist1);
% distances=pdist2(Data.S,Data.S); %compute pairwise dstances between points in S, pdist2 is MATLAB built-in function  成对的距离
if strcmp(rbf_flag,'cub') %cubic RBF  立方
    PairwiseDistance=distances.^3; 
elseif strcmp(rbf_flag,'lin') %linear RBF  线性
    PairwiseDistance=distances;    
elseif strcmp(rbf_flag,'tps') %thin-plate spline RBF  薄板样条
    PairwiseDistance=distances.^2.*log(distances+realmin);
    PairwiseDistance(logical(eye(size(distances)))) = 0;   %???
end


PHI(1:Data.m,1:Data.m)=PairwiseDistance; %matrix with radial bsais function values
P = [Data.S,ones(Data.m,1)];% sample site matrix equipped with vector of ones
[m,n]=size(Data.S); %determine how many points are in S and what is the dimension


A=[PHI,P;P', zeros(n+1,n+1)]; %set up matrix for solving for parameters  按模型列
RHS=[Data.Y;zeros(n+1,1)]; %right hand side of linear system                 %Data.Y？
params=A\RHS; %compute parameters
lambda=params(1:m); %weights     有m个lambda
gamma=params(m+1:end); %parameters of polynomial tail


end





