clc;clear;
n=3;
n_S=2;
n_V=3;
n_T=3;
M=999999;

load (strcat('problems\prob_',num2str(n),'_',num2str(n_S),'.mat'));


%Define variables
Q=binvar(n,n_S-1,5,5,n_T,'full');%transportation
U=binvar(n,n_S,5,n_V,'full'); % processing and speed
W=binvar(n,n,n_S,5,'full'); %direct precedence
W3=binvar(n,n,n_S,'full'); %direct precedence
X=binvar(n,n_S,'full'); %turn off


S=sdpvar(n,n_S,'full'); %starting times
C=sdpvar(n,n_S,'full'); %completion time
A=sdpvar(n,n_S,'full'); % arrival time
D=sdpvar(n,n_S,'full'); % OPA end time
T_I=sdpvar(n,n,n_S,'full'); % idle time after o_ij
T_minus=sdpvar(n,n_S,'full'); % earliness time
T_plus=sdpvar(n,n_S,'full'); % tardiness time



E_T=sdpvar(n,n_S-1,'full');% transportation cost
elQ=sdpvar(n,n_S-1,5,5,n_T,'full');% temp
bU=sdpvar(n,n_S,5,'full');% temp
E_P=sdpvar(n,n_S,'full');%processing energy
E_I=sdpvar(n,n_S,'full');% idle cost
E_SI=sdpvar(n,n_S,'full');%turning or idle cost
E_ET=sdpvar(n,n_S,'full'); % the overall ET cost


%Constaints
% Constraints=[0<=S<=M];
% Constraints=[Constraints,0<=C<=M];
% Constraints=[Constraints,0<=D<=M];
% Constraints=[Constraints,0<=A<=M];
% Constraints=[Constraints,0<=T_minus<=M];
% Constraints=[Constraints,0<=T_plus<=M];
% Constraints=[Constraints,0<=T_I<=M];
% Constraints=[Constraints,0<=E_T<=M*999];
% Constraints=[Constraints,0<=elQ<=M*999];
% Constraints=[Constraints,0<=bU<=M*999];
% Constraints=[Constraints,0<=E_P<=M*999];
% Constraints=[Constraints,0<=E_I<=M*999];
% Constraints=[Constraints,0<=E_SI<=M*999];
% Constraints=[Constraints,0<=E_ET<=M*999];

%Constaints
Constraints=[0<=S];
Constraints=[Constraints,0<=C];
Constraints=[Constraints,0<=D];
Constraints=[Constraints,0<=A];
Constraints=[Constraints,0<=T_minus];
Constraints=[Constraints,0<=T_plus];
Constraints=[Constraints,0<=T_I];
Constraints=[Constraints,0<=E_T];
Constraints=[Constraints,0<=elQ];
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


for i=1:n
    for j=1:n_S
        for m=1:5
            for v=1:n_V
                Constraints=[Constraints, ...
                    implies(U(i,j,m,v),S(i,j)==C(i,j)-p(i,j,m,v)), ...
                    implies(U(i,j,m,v),E_P(i,j)==e(j,m,v)*p(i,j,m,v)), ...
                    implies(U(i,j,m,v),E_I(i,j)==a(j,m,v)*sum(T_I(i,:,j))), ...
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
                bU(i,j,m)==b(j,m)*sum(U(i,j,m,:),"all"), ...
                ];
        end
    end
end

for j=1:n_S
    for m=1:5
        Constraints=[Constraints, ...
            sum(W(:,:,j,m),"all")>=sum(U(:,j,m,:),"all")-1, ...
            ];
    end
end

for i=1:n
    for j=1:n_S
        Constraints=[Constraints, ...
            sum(U(i,j,:,:),"all")==1, ...
            D(i,j)==A(i,j)+OPA(i,j), ...
            S(i,j)>=A(i,j), ...
            T_minus(i,j)>=D(i,j)-window_width(i,j)/2-S(i,j), ...
            T_plus(i,j)>=S(i,j)-D(i,j)-window_width(i,j)/2, ...
            E_ET(i,j)==T_minus(i,j)*alpha_(i,j)+T_plus(i,j)*beta_(i,j), ...
            1-sum(W(i,:,j,:),"all")<=X(i,j), ...
            implies(1-X(i,j),E_SI(i,j)==E_I(i,j)), ...
            implies(X(i,j),E_SI(i,j)==sum(bU(i,j,:),"all")), ...
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
                        implies(Q(i,j,m,m1,t),A(i,j+1)==C(i,j)+l(j,m,m1,t)), ...
                        elQ(i,j,m,m1,t)==e_t(t)*l(j,m,m1,t)*Q(i,j,m,m1,t), ...
                        ];
                end
            end
        end
    end
end


for i=1:n
    for j=1:n_S-1
        Constraints=[Constraints, ...
            E_T(i,j)==sum(elQ(i,j,:,:,:),"all"), ...
            ];
    end
end

Constraints=[Constraints, ...
    sum(A(:,1))==0, ...
    ];



% Define the objective
% Objective = 1;
Objective = sum(E_P,'all')+sum(E_T,'all')+sum(E_SI,'all')+sum(E_ET,'all');

% Set some options for YALMIP and solver
options = sdpsettings('solver','cplex','verbose',1);

% Solve the problem
tic
sol = optimize(Constraints,Objective,options);
solve_time=toc;


% Analyze error flags
if sol.problem == 0
    % Extract and display value
    S=value(S);
    C=value(C);
    D=value(D);
    A=value(A);
    T_minus=value(T_minus);
    T_plus=value(T_plus);
    T_I=value(T_I);
    U=value(U);
    W=value(W);
    X=value(X);
    Q=value(Q);
    E_T=value(E_T);
    E_P=value(E_P);
    E_I=value(E_I);
    E_SI=value(E_SI);
    E_ET=value(E_ET);

    obj=value(Objective)
else
    disp('Hmm, something went wrong!')
    sol.info
    yalmiperror(sol.problem)
end


% save (strcat('cplex_ans\prob_',num2str(n),'_',num2str(n_S),'.mat'));


% jobId=zeros(n,n_S);
% startT=zeros(n,n_S);
% durationT=zeros(n,n_S);
% 
% for i=1:n
%     for j=1:n_S
%         for k=1:5
%             if (1-m<=U(i,j,k))&&(U(i,j,k)<=1+m)
%                 jobId(i,j)=sum(machine(1:j-1))+k;
%             end
%         end
%         startT(i,j)=S(i,j);
%         durationT(i,j)=C(i,j)-S(i,j);
%     end
% end
% 
% jobId=reshape(jobId',1,[]);
% startT=reshape(startT',1,[]);
% durationT=reshape(durationT',1,[]);
% 
% 
% pName{length(jobId)}='';
% for i=1:length(jobId)
%     pName(i)={[num2str(floor((i-1)/n_S)+1),'-',num2str(mod(i-1,n_S)+1)]};
% end
% GTC=ganttChart(startT,durationT,jobId,'String',pName);
% ax=gca;
% ax.YTickLabel={'M1-1','M1-2','M1-3','M2-1','M2-2','M2-3','M2-4','M2-5', ...
%     'M3-1','M3-2','M3-3','M3-4','M3-5','M4-1','M4-2','M4-3','M4-4', ...
%     'M5-1','M5-2','M5-3','M6-1','M6-2','M6-3','M7-1','M7-2','M8-1','M8-2','M8-3','M8-4','M8-5'};
% 
% 
% 
% 









