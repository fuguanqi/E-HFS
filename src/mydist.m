function dist_matrix = mydist(Samples,Data)
%MYDIST Summary of this function goes here
%   Detailed explanation goes here

n=size(Samples,1);
dist_matrix=zeros(n,n);
for i=1:n
    for j=i:n
        if i==j
            dist_matrix(i,j)=0;
        end
        dist=cal_dist(Samples(i,:),Samples(j,:),Data);
        dist_matrix(i,j)=dist;
        dist_matrix(j,i)=dist;
    end
end
end

function dist = cal_dist(x,y,Data)
res=y-x;
res=res.^2;
log=res~=0;
res(Data.category(1):Data.category(end))=double(log(Data.category(1):Data.category(end)));
dist=sqrt(sum(res));

% n=numel(x);
% dist=0;
% for i=1:n
%     if i>=Data.category(1) && i<=Data.category(end)
%         if x(i)~=y(i)
%             dist=dist+1;
%         end
%     else
%         dist=dist+(y(i)-x(i))^2;
%     end
% end
% dist=sqrt(dist);
end