clc;clear;
warning("off");
n_array=[20,40,60,80,100,120];
repeat=5;
n_S_array=[1,2,4,8];
n_V=3;
n_T=3;
for rr=1:repeat
    for n_S=n_S_array
        for n=n_array
            prob=load (strcat('problems\prob_',num2str(n),'_',num2str(n_S),'.mat'));

            prob.n_M=prob.machine(1:n_S);
            n_M=prob.n_M;

            int_dim=n*n_S+n*(n_S-1);%  assign,  tr_speed
            con_dim=n*n_S;                   % seq
            Dim=int_dim+con_dim;

            Data.prob=prob;

            Data.number_startpoints=2*(Dim+1);
            Data.dim=Dim;
            Data.xlow=[zeros(1,con_dim),ones(1,int_dim)];

            Data.xup=[repmat([1],1,con_dim)];
            for i= 1:n_S
                Data.xup=[Data.xup,repelem(n_M(i),n)];
            end
            Data.xup=[Data.xup,repelem(n_T,n*(n_S-1))];

            Data.integer=[con_dim+1:Dim]; %indices of integer variables
            Data.category=[(n)*n_S+1:(n)*n_S+n*n_S]; %indices of category integer variables
            % Data.category=[1:(n-1)*n_S+n*n_S]; %indices of category integer variables
            % Data.category=[];
            Data.continuous=[1:con_dim]; %indices of continuous variables

            InitialPoints = slhd(Data);
            xlatin=repmat(Data.xlow, Data.number_startpoints,1) + repmat(Data.xup-Data.xlow,Data.number_startpoints,1).*InitialPoints;
            Data.S=xlatin;
            fixrank = false;
            if rank([Data.S,ones(size(Data.S,1),1)]) < Data.dim + 1
                fixrank = true;
            end
            while fixrank %rank is too small to fit RBF, add random points to initial design to satisfy rank condition
                n_new = Data.dim+1-rank([Data.S,ones(size(Data.S,1),1)]); %minimum number of new points needed
                randpoint = repmat(Data.xlow,n_new,1) + repmat((Data.xup-Data.xlow), n_new,1).*rand(n_new,Data.dim);
                temp=rank([[Data.S;randpoint], ones(size(Data.S,1)+n_new,1)]);
                if temp == Data.dim + 1
                    Data.S = [Data.S; randpoint];
                    fixrank = false;
                end
            end
            Data.S(:,Data.integer)=round(Data.S(:,Data.integer));


            % save (strcat('temp_data\Data1.mat'));
            for i=1:ceil(n*n_S/5)
                greedy_x=greedy_sol_D(Data.prob);
                Data.S=[Data.S;greedy_x];
            end

            for i=1:ceil(n*n_S/5)
                greedy_x=greedy_sol_Dl(Data.prob);
                Data.S=[Data.S;greedy_x];
            end

            for i=1:ceil(n*n_S/5)
                greedy_x=greedy_sol_Du(Data.prob);
                Data.S=[Data.S;greedy_x];
            end


            Iteration=Dim*100;


            % miso('datainput_dp',Iteration, 'rbf_c', [], 'slhd', 'cp4',[],Data); %SODA-ADM
            [xbest, fbest]=miso('datainput_real_dp_speed',Iteration, 'rbf_c', [], 'slhd', 'soda_adm_fu',[],Data); %the new SODA-ADM
            % [xbest, fbest] = miso('datainput_dp',Iteration, 'rbf_c', [], 'slhd', 'cp6',[],Data); %SODA-ADM-DP
        
        writematrix(xbest,strcat('datainput_real_dp_speed/best_X_n',num2str(Data.prob.n),'_',num2str(Data.prob.n_S)), 'WriteMode','append');

end
    end
    
end

function sol=greedy_sol_D(prob)
n=prob.n;
n_S=prob.n_S;
n_M=prob.n_M;
A=zeros(n,n_S);
A(:,1)=prob.r;
y=0;
C=zeros(n,n_S);
D=[prob.r+prob.OPA(1) [1:prob.n]'];
sol_seq=[];
sol_assign=[];
sol_pv=[];
sol_tv=ones(1,prob.n*(prob.n_S-1));
assign_next=[randi(prob.n_M(1),1,prob.n)];

for s=1:prob.n_S
    D=sortrows(D);
    stage_seq=D(:,2);
    stage_seq=[stage_seq [1:n]'];
    stage_seq=sortrows(stage_seq);
    stage_seq=(stage_seq(:,2)-1).*(1/n);
    sol_seq=[sol_seq stage_seq'];
    assign=assign_next;
    sol_assign=[sol_assign,assign];
    pv=randi(3,1,n_M(s));
    sol_pv=[sol_pv pv];
    [stage_C,stage_EP,stage_ETI]=timingByDPprspeed(D(:,2),assign,A(:,s),s,prob);
    C(:,s)=stage_C;
    if s<n_S
        l=prob.l;
       
        assign_next=[randi(prob.n_M(s+1),1,prob.n)];
        for i=1:n
            A(i,s+1)=C(i,s)+l(s,assign(i),assign_next(i),1);
        end
        D=[A(:,s+1)+prob.OPA(s) [1:prob.n]'];
    end
end
sol=[sol_seq sol_assign sol_tv];

end


function sol=greedy_sol_Du(prob)
n=prob.n;
n_S=prob.n_S;
n_M=prob.n_M;
A=zeros(n,n_S);
A(:,1)=prob.r;
y=0;
C=zeros(n,n_S);
D=[prob.r+prob.OPA(1)+prob.window_width(1) [1:prob.n]'];
sol_seq=[];
sol_assign=[];
sol_tv=ones(1,prob.n*(prob.n_S-1));
assign_next=[randi(prob.n_M(1),1,prob.n)];
sol_pv=[];
for s=1:prob.n_S
    D=sortrows(D);
    stage_seq=D(:,2);
    stage_seq=[stage_seq [1:n]'];
    stage_seq=sortrows(stage_seq);
    stage_seq=(stage_seq(:,2)-1).*(1/n);
    sol_seq=[sol_seq stage_seq'];
    assign=assign_next;
    sol_assign=[sol_assign,assign];
    pv=randi(3,1,n_M(s));
    sol_pv=[sol_pv pv];
    [stage_C,stage_EP,stage_ETI]=timingByDPprspeed(D(:,2),assign,A(:,s),s,prob);
    C(:,s)=stage_C;
    if s<n_S
        l=prob.l;
       
        assign_next=[randi(prob.n_M(s+1),1,prob.n)];
        for i=1:n
            A(i,s+1)=C(i,s)+l(s,assign(i),assign_next(i),1);
        end
        D=[A(:,s+1)+prob.OPA(s)+prob.window_width(s) [1:prob.n]'];
    end
end
sol=[sol_seq sol_assign sol_tv];

end

function sol=greedy_sol_Dl(prob)
n=prob.n;
n_S=prob.n_S;
n_M=prob.n_M;
A=zeros(n,n_S);
A(:,1)=prob.r;
y=0;
C=zeros(n,n_S);
D=[prob.r+prob.OPA(1)-prob.window_width(1) [1:prob.n]'];
sol_seq=[];
sol_assign=[];
sol_tv=ones(1,prob.n*(prob.n_S-1));
assign_next=[randi(prob.n_M(1),1,prob.n)];
sol_pv=[];
for s=1:prob.n_S
    D=sortrows(D);
    stage_seq=D(:,2);
    stage_seq=[stage_seq [1:n]'];
    stage_seq=sortrows(stage_seq);
    stage_seq=(stage_seq(:,2)-1).*(1/n);
    sol_seq=[sol_seq stage_seq'];
    assign=assign_next;
    sol_assign=[sol_assign,assign];
    pv=randi(3,1,n_M(s));
    sol_pv=[sol_pv pv];
    [stage_C,stage_EP,stage_ETI]=timingByDPprspeed(D(:,2),assign,A(:,s),s,prob);
    C(:,s)=stage_C;
    if s<n_S
        l=prob.l;
       
        assign_next=[randi(prob.n_M(s+1),1,prob.n)];
        for i=1:n
            A(i,s+1)=C(i,s)+l(s,assign(i),assign_next(i),1);
        end
        D=[A(:,s+1)+prob.OPA(s)-prob.window_width(s) [1:prob.n]'];
    end
end
sol=[sol_seq sol_assign sol_tv];

end