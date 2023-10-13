function y=datainput_direct(xx,prob) %objective function
x=xx'; % make sure vector is row vector

%test
% x=[2 3 4 5 2 3 4 5 1 1 2 3 3 1 2 3 4 5 1 1 1 1 1 1 1 1 1 1 1 1 1];
% x=[2:20 1 1 2 3 3 1 1 2 3 3 1 1 2 3 3 1 1 2 3 3 2 3 3 ];
%test


n=prob.n;
n_S=prob.n_S;
n_M=prob.n_M;
n_V=3;
n_T=3;

W=binvar(n,n,n_S,5,'full'); %direct precedence
W3=binvar(n,n,n_S,'full'); %direct precedence
X=binvar(n,n_S,'full'); %turn off
U=binvar(n,n_S,5,n_V,'full'); % assign and speed
V=binvar(n_S,5,n_V,'full'); % processing speed

S=sdpvar(n,n_S,'full'); %starting times
C=sdpvar(n,n_S,'full'); %completion time
A=sdpvar(n,n_S,'full'); % arrival time
D=sdpvar(n,n_S,'full'); % OPA end time
T_I=sdpvar(n,n,n_S,'full'); % idle time after o_ij
T_minus=sdpvar(n,n_S,'full'); % earliness time
T_plus=sdpvar(n,n_S,'full'); % tardiness time

bU=sdpvar(n,n_S,5,'full');% temp
E_P=sdpvar(n,n_S,'full');%processing energy
E_I=sdpvar(n,n_S,'full');% idle cost
E_SI=sdpvar(n,n_S,'full');%turning or idle cost
E_ET=sdpvar(n,n_S,'full'); % the overall ET cost

Constraints=[0<=S];
Constraints=[Constraints,0<=C];
Constraints=[Constraints,0<=D];
Constraints=[Constraints,0<=A];
Constraints=[Constraints,0<=T_minus];
Constraints=[Constraints,0<=T_plus];
Constraints=[Constraints,0<=T_I];
Constraints=[Constraints,0<=bU];
Constraints=[Constraints,0<=E_P];
Constraints=[Constraints,0<=E_I];
Constraints=[Constraints,0<=E_SI];
Constraints=[Constraints,0<=E_ET];

if n_S>=1
    for i=1:n
        Constraints=[Constraints, ...
            sum(U(i,1,4,:))==0,...
            sum(U(i,1,5,:))==0];
    end
end

if n_S>=4
    for i=1:n
        Constraints=[Constraints, ...
            sum(U(i,4,5,:))==0];
    end
end

if n_S>=5
    for i=1:n
        Constraints=[Constraints, ...
            sum(U(i,5,4,:))==0, ...
            sum(U(i,5,5,:))==0];
    end
end

if n_S>=6
    for i=1:n
        Constraints=[Constraints, ...
            sum(U(i,6,4,:))==0, ...
            sum(U(i,6,5,:))==0];
    end
end

if n_S>=7
    for i=1:n
        Constraints=[Constraints, ...
            sum(U(i,7,3,:))==0, ...
            sum(U(i,7,4,:))==0, ...
            sum(U(i,7,5,:))==0];
    end
end
if n_S>1
    Q=binvar(n,n_S-1,5,5,n_T,'full');%transportation
    E_T=sdpvar(n,n_S-1,'full');% transportation cost
    elQ=sdpvar(n,n_S-1,5,5,n_T,'full');% temp
    Constraints=[Constraints,0<=E_T];
    Constraints=[Constraints,0<=elQ];
    for i=1:n
        for j=1:n_S-1
            Constraints=[Constraints, ...
                E_T(i,j)==sum(elQ(i,j,:,:,:),"all"), ...
                ];
        end
    end
    for i=1:n
        for j=1:n_S-1
            for m=1:5
                for m1=1:5
                    Constraints=[Constraints, ...
                        sum(U(i,j,m,:),"all")>=sum(Q(i,j,m,m1,:),"all"), ...
                        sum(U(i,j+1,m1,:),"all")>=sum(Q(i,j,m,m1,:),"all"), ...
                        sum(U(i,j,m,:),"all")+sum(U(i,j+1,m1,:),"all")<=sum(Q(i,j,m,m1,:),"all")+1, ...
                        ];
                    for t=1:n_T
                        Constraints=[Constraints, ...
                            implies(Q(i,j,m,m1,t),A(i,j+1)==C(i,j)+prob.l(j,m,m1,t)), ...
                            elQ(i,j,m,m1,t)==prob.e_t(t)*prob.l(j,m,m1,t)*Q(i,j,m,m1,t), ...
                            ];
                    end
                end
            end
        end
    end

end

for i=1:n
    for j=1:n_S
        for m=1:5
            for v=1:n_V
                Constraints=[Constraints, ...
                    implies(U(i,j,m,v),S(i,j)==C(i,j)-prob.p(i,j,m,v)), ...
                    implies(U(i,j,m,v),E_P(i,j)==prob.e(j,m,v)*prob.p(i,j,m,v)), ...
                    implies(U(i,j,m,v),E_I(i,j)==prob.a(j,m,v)*sum(T_I(i,:,j))), ...
                    U(i,j,m,v)<=V(j,m,v), ...
                    ];
            end
        end
    end
end

for i=1:n
    for i1=1:n
        if i==i1
            Constraints=[Constraints, ...
                W(i,i1,:,:)==0, ...
                ];
        end
        for j=1:n_S
            Constraints=[Constraints, ...
                W3(i,i1,j)==sum(W(i,i1,j,:),"all"), ...
                implies(W3(i,i1,j),S(i1,j)>=C(i,j)), ...
                implies(W3(i,i1,j),T_I(i,i1,j)==S(i1,j)-C(i,j)), ...
                implies(1-W3(i,i1,j),T_I(i,i1,j)==0), ...
                ];
        end
    end
end

for i=1:n
    for j=1:n_S
        for m=1:5
            Constraints=[Constraints, ...
                sum(W(i,:,j,m),"all")<=sum(U(i,j,m,:),"all"), ...
                sum(W(:,i,j,m),"all")<=sum(U(i,j,m,:),"all"), ...
                bU(i,j,m)==prob.b(j,m)*sum(U(i,j,m,:),"all"), ...
                ];
        end
    end
end

for j=1:n_S
    for m=1:5
        Constraints=[Constraints, ...
            sum(W(:,:,j,m),"all")>=sum(U(:,j,m,:),"all")-1, ...
            sum(V(j,m,:),"all")==1, ...
            ];
    end
end

for i=1:n
    for j=1:n_S
        Constraints=[Constraints, ...
            sum(U(i,j,:,:),"all")==1, ...
            D(i,j)==A(i,j)+prob.OPA(i,j), ...
            S(i,j)>=A(i,j), ...
            T_minus(i,j)>=D(i,j)-prob.window_width(i,j)/2-S(i,j), ...
            T_plus(i,j)>=S(i,j)-D(i,j)-prob.window_width(i,j)/2, ...
            E_ET(i,j)==T_minus(i,j)*prob.alpha_(i,j)+T_plus(i,j)*prob.beta_(i,j), ...
            1-sum(W(i,:,j,:),"all")<=X(i,j), ...
            implies(1-X(i,j),E_SI(i,j)==E_I(i,j)), ...
            implies(X(i,j),E_SI(i,j)==sum(bU(i,j,:),"all")), ...
            ];
    end
end





Constraints=[Constraints, ...
    sum(A(:,1))==0, ...
    ];





for j=1:n_S
    mac_seqs=cell(prob.n_M(j),1);
    job_seq=[0,1,0];
    for i=1:n-1
        job_seq=[job_seq(1:x((j-1)*(n-1)+i)) i+1 job_seq(x((j-1)*(n-1)+i)+1:end)];
    end
    job_seq=job_seq(2:n+1);
    assign=x(n_S*(n-1)+(j-1)*n+1:n_S*(n-1)+j*n);
    pr_speed=x((n-1)*n_S+n*n_S+sum(n_M(1:j-1))+1:(n-1)*n_S+n*n_S+sum(n_M(1:j)));
    for i=1:prob.n
        mac_seqs{assign(job_seq(i))}=[mac_seqs{assign(job_seq(i))},job_seq(i)];
    end

    for m=1:prob.n_M(j)
        mn=numel(mac_seqs{m});
        m_seq=mac_seqs{m};
        Constraints=[Constraints, ...
            V(j,m,pr_speed(m))==1, ...
            ];
        if mn>=1
            % [mC,m_EP,m_ETI]=cal_1m(mac_seqs{m},j,m,pr_speeds(m),A,prob);
            for i=1:mn-1
                Constraints=[Constraints, ...
                    W(m_seq(i),m_seq(i+1),j,m)==1, ...
                    ];
            end

        end
    end


    if j<n_S
        tr_speed=x((n-1)*n_S+n*n_S+sum(n_M)+(j-1)*n+1:(n-1)*n_S+n*n_S+sum(n_M)+j*n);
        assign_next=x(n_S*(n-1)+j*n+1:n_S*(n-1)+(j+1)*n);
        
        for i=1:n
            Constraints=[Constraints, ...
                Q(i,j,assign(i),assign_next(i),tr_speed(i))==1, ...
            ];
            
        end
    end

end

% Define the objective
% Objective = 1;
if n_S>1
    Objective = sum(E_P,'all')+sum(E_T,'all')+sum(E_SI,'all')+sum(E_ET,'all');
else
    Objective = sum(E_SI,'all')+sum(E_ET,'all')+sum(E_P,'all');
end

% Set some options for YALMIP and solver
options = sdpsettings('solver','cplex','verbose',1);

% Solve the problem
sol = optimize(Constraints,Objective,options);

% Analyze error flags
if sol.problem == 0

    y=value(Objective);
else
    disp('Hmm, something went wrong!')
    sol.info
    yalmiperror(sol.problem)
end




end %
