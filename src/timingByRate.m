function [C,stage_EP,stage_ETI]=timingByRate(job_seq,assign,pr_speeds,A,j,ST_rate,prob)
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
        [mC,m_EP,m_ETI]=cal_1m(mac_seqs{m},j,m,pr_speeds(m),A,ST_rate,prob);
        stage_EP=stage_EP+m_EP;
        stage_ETI=stage_ETI+m_ETI;
        C(mac_seqs{m})=mC;
    end
end

end

function [C,m_EP,m_ETI]=cal_1m(job_seq,j,m,v,A,ST_rate,prob)
A=A(job_seq);
ST_rate=ST_rate(job_seq);
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
m_ETI=0;
earliest_S(1)=A(1);
S=zeros(numel(job_seq),1);
C=zeros(numel(job_seq),1);
for i=1:numel(job_seq)
    if i>1
        earliest_S(i)=max(A(i),C(i-1)+p(i-1));
    end
    if ST_rate(i)>1.0
        S(i)=earliest_S(i)+OPA(i)*(1/(2.001-ST_rate(i)));
    else
        S(i)=earliest_S(i)+OPA(i)*ST_rate(i);
    end
    m_ETI=m_ETI+max(0,alpha_(i)*(d_minus(i)-S(i)))+max(0,beta_(i)*(S(i)-d_plus(i)));
    C(i)=S(i)+p(i);
    if i>1
        m_ETI=m_ETI+min(b,a*(S(i)-C(i-1)));
    end
end
m_ETI=m_ETI+b;
end

