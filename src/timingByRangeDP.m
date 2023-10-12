function [C,stage_EP,stage_ETI]=timingByRate(job_seq,assign,pr_speeds,A,j,r,prob)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
stage_EP=0;
stage_ETI=0;
C=zeros(prob.n,1);
mac_seqs=cell(prob.n_M(j),1);

for i=1:prob.n
    mac_seqs{assign(job_seq(i))}=[mac_seqs{assign(job_seq(i))},job_seq(i)];
end

for m=1:prob.n_M(j)
    mn=numel(mac_seqs{m});
    if mn>=1
        [mC,m_EP,m_ETI]=cal_1m(mac_seqs{m},j,m,pr_speeds(m),A,r,prob);
        stage_EP=stage_EP+m_EP;
        stage_ETI=stage_ETI+m_ETI;
        C(mac_seqs{m})=mC;
    end
end

end

function [C,m_EP,m_ETI]=cal_1m(job_seq,j,m,v,A,r,prob)
A=A(job_seq);
r=r(job_seq);
earliest_S=A;
OPA=prob.OPA(job_seq,j);
D=A+OPA;
d_minus=max(0,D-prob.window_width(job_seq,j)/2);
d_plus=D+prob.window_width(job_seq,j)/2;
alpha_=prob.alpha_(job_seq,j);
beta_=prob.beta_(job_seq,j);
a=prob.a(j,m,v);
b=prob.b(j,m);
p=prob.p(job_seq,j,m,v);
m_EP=prob.e(j,m,v).*prob.p(job_seq,j,m,v);
m_EP=sum(m_EP,"all");
memos=cell(numel(job_seq),1);
m_ETI=0;
f=[];
earliest_S(1)=A(1);

for i=1:numel(job_seq)
    if i>1
        earliest_S(i)=max(A(i),earliest_S(i-1)+p(i-1));
    end
    f=dp(earliest_S(i),d_minus(i),d_plus(i),alpha_(i),beta_(i),a,b,f);
    
    memos{i}=f;
    f(1,:)=f(1,:)+p(i);
end
ub=f(1,end-1);
% m_ETI=m_ETI+min(f(2,:))+b;
C=zeros(numel(job_seq),1);

for i=numel(job_seq):-1:1
    ub=ub-p(i);
    f=memos{i};
    rt=f(1,1)+r(i)*(ub-f(1,1));
    idx=find(f(1,:)<=rt+0.001);
    slope=(f(2,idx(end)+1)-f(2,idx(end)))/(f(1,idx(end)+1)-f(1,idx(end)));
    if slope==0
        C(i)=rt+p(i);
        ub=rt;
        m_ETI=m_ETI+b+f(2,idx(end));
    elseif slope>0
        if rt-f(1,idx(end))>constants.lambda_*(f(1,idx(end)+1)-f(1,idx(end)))
            C(i)=rt+p(i);
            ub=rt;
            if i==numel(job_seq)
                m_ETI=m_ETI+b+f(2,idx(end))+slope*(rt-f(1,idx(end)));
            end
        else
            C(i)=f(1,idx(end))+p(i);
            ub=f(1,idx(end));
            if i==numel(job_seq)
                m_ETI=m_ETI+b+f(2,idx(end));
            end
        end
    elseif slope<0
        if f(1,idx(end)+1)-rt>constants.lambda_*(f(1,idx(end)+1)-f(1,idx(end)))
            C(i)=rt+p(i);
            ub=rt;
            if i==numel(job_seq)
                m_ETI=m_ETI+b+f(2,idx(end))+slope*(rt-f(1,idx(end)));
            end
        else
            C(i)=f(1,idx(end)+1)+p(i);
            ub=f(1,idx(end)+1);
            if i==numel(job_seq)
                m_ETI=m_ETI+b+f(2,idx(end)+1);
            end
        end
    end
    if m_ETI<0
        error("error");
    end
    % f=f(:,idx);
    % f=[f,[ub;f(2,end)+slope*(ub-f(1,end))]];
    % x=find(f(2,:)==min(f(2,:)));
    % ub=f(1,(x(1)));
    % C(i)=f(1,x(1))+p(i);
end
end



function f=dp(earliest_S,d_minus,d_plus,alpha_,beta_,a,b,pre_f)
M=constants.M;
if size(pre_f,2)==0
    f=[0,d_minus,d_plus,M;d_minus*alpha_,0,0,beta_*(M-d_plus)];
    f=adjustST(f,earliest_S);
    return;
end
i=1;
while i<size(pre_f,2)
    if pre_f(2,i+1)>pre_f(2,i)
        new_f=minOf(pre_f(:,i:end),a,b);
        pre_f=[pre_f(:,1:i),new_f(:,2:end)];
    end
    i=i+1;
    pre_f=clear_func(pre_f);
end

for i=1:size(pre_f,2)-1
    if pre_f(1,i)<d_minus && pre_f(1,i+1)>d_minus
        pre_f=[pre_f(:,1:i),[d_minus;pre_f(2,i)+(d_minus-pre_f(1,i))*((pre_f(2,i+1)-pre_f(2,i))/(pre_f(1,i+1)-pre_f(1,i)))],pre_f(:,i+1:end)];
    end
end
for i=1:size(pre_f,2)-1
    if pre_f(1,i)<d_plus && pre_f(1,i+1)>d_plus
        pre_f=[pre_f(:,1:i),[d_plus;pre_f(2,i)+(d_plus-pre_f(1,i))*((pre_f(2,i+1)-pre_f(2,i))/(pre_f(1,i+1)-pre_f(1,i)))],pre_f(:,i+1:end)];
    end
end

for i=1:size(pre_f,2)
    if pre_f(1,i)<d_minus
        pre_f(2,i)=pre_f(2,i)+(d_minus-pre_f(1,i))*alpha_;
    elseif pre_f(1,i)>d_plus
        pre_f(2,i)=pre_f(2,i)+(pre_f(1,i)-d_plus)*beta_;
    end
end
f=adjustST(pre_f,earliest_S);
pre_f=clear_func(pre_f);
if pre_f(1,end-1)>constants.M*0.5
        error("error");
    end
end


function f=adjustST(old_f,ES)
% old_f(1,:)
% ES
if ES==old_f(1,1)
    f=old_f;
    return;
end
ind=find(old_f(1,:)>ES);
f=old_f(:,ind);
f=[[ES;f(2,1)-(f(1,1)-ES)*((f(2,1)-old_f(2,ind(1)-1))/(f(1,1)-old_f(1,ind(1)-1)))],f];
% f=clear_func(f);
end

function min_f=minOf(f,a,b)
y0=f(2,1);
x0=f(1,1);
x1=f(1,1)+b/a;
y1=f(2,1)+b;
x3=constants.M;
y3=y1;
min_f=[];
for i=2:size(f,2)
    % if e_SI(f(1,i-1),y0,x0,a,b)==f(2,i-1)
    %     if f(1,i-1)<x1
    %         min_f=[min_f,[f(1:2,i-1);min(b/a,f(3,i-1))]];
    %     else
    %         min_f=[min_f,[f(1:2,i-1);min(0,f(3,i-1))]];
    %     end
    %     flag=0;
    % elseif e_SI(f(1,i-1),y0,x0,a,b)>f(2,i-1)
    %     min_f=[min_f,[f(1:3,i-1)]];
    %     flag=1;
    % else
    %     flag=-1;
    % end
    min_f=[min_f,[f(1,i-1);min(f(2,i-1),e_SI(f(1,i-1),y0,x0,a,b))]];
    if f(1,i-1)<x1
        [is_intersect1,x,y]=intersect_at(f(1,i-1),f(2,i-1),f(1,i),f(2,i),x0,y0,x1,y1);
        if is_intersect1==1
            min_f=[min_f,[x;y]];
        end

        [is_intersect2,x,y]=intersect_at(f(1,i-1),f(2,i-1),f(1,i),f(2,i),x1,y1,x3,y3);
        if is_intersect2==1
            min_f=[min_f,[x1;min(y1,f(2,i-1)+(x1-f(1,i-1))*(f(2,i)-f(2,i-1))/(f(1,i)-f(1,i-1)))]];
            min_f=[min_f,[x;y]];
        end
        if is_intersect1==0 && is_intersect2==0
            min_f=[min_f,[x1;y1]];
        end
    end
    if f(1,i-1)>=x1
        [is_intersect2,x,y]=intersect_at(f(1,i-1),f(2,i-1),f(1,i),f(2,i),x1,y1,x3,y3);
        if is_intersect2==1
            min_f=[min_f,[x;y]];
        end
    end
end
min_f=[min_f,[f(1,end);min(f(2,end),e_SI(f(1,end),y0,x0,a,b))]];
min_f=clear_func(min_f);
end

function y=e_SI(x,y0,x0,a,b)
y=min(y0+a*(x-x0)+b,y0+b);
end

function f=clear_func(f)
if size(f,2)<3
    return;
end
indToDel=[];
for i=1:size(f,2)-1
    if abs(f(1,i+1)-f(1,i))<0.0001
        indToDel=[indToDel,i+1];
    end
    if f(1,i)>0.5*constants.M
        indToDel=[indToDel,i];
    end
end
f(:,indToDel)=[];

indToDel=[];
for i=1:size(f,2)-2
    if  abs((f(2,i+1)-f(2,i))/(f(1,i+1)-f(1,i)) - (f(2,i+2)-f(2,i+1))/(f(1,i+2)-f(1,i+1)))<0.001
        indToDel=[indToDel,i+1];
    end
end
f(:,indToDel)=[];
end







