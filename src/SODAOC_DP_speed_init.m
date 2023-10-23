clc;clear;
warning("off");
n_array=[20,40,60,80,100,120];
repeat=15;
n_S_array=[1,2,4,8];
n_V=3;
n_T=3;
for rr=1:repeat
    for n_S=4
        for n=20
            prob=load (strcat('problems\prob_',num2str(n),'_',num2str(n_S),'.mat'));

            prob.n_M=prob.machine(1:n_S);
            n_M=prob.n_M;

            int_dim=(n-1)*n_S+n*n_S+n*(n_S-1);% seq, assign, tr_speed
            con_dim=0;
            Dim=int_dim+con_dim;

            Data.prob=prob;

            Data.number_startpoints=2*(Dim+1);
            Data.dim=Dim;
            Data.xlow=[ones(1,Dim)];

            Data.xup=[repmat(2:n,1,n_S)];
            for i= 1:n_S
                Data.xup=[Data.xup,repelem(n_M(i),n)];
            end
            Data.xup=[Data.xup,repelem(n_T,n*(n_S-1))];

            Data.integer=[1:Dim]; %indices of integer variables
            Data.category=[(n-1)*n_S+1:(n-1)*n_S+n*n_S]; %indices of category integer variables
            % Data.category=[1:(n-1)*n_S+n*n_S]; %indices of category integer variables
            % Data.category=[];
            Data.continuous=[]; %indices of continuous variables

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

            for i=1:ceil(n*n_S/10)
                greedy_x=greedy_sol_D(Data.prob);
                Data.S=[Data.S;greedy_x];
            end

            for i=1:ceil(n*n_S/10)
                greedy_x=greedy_sol_Dl(Data.prob);
                Data.S=[Data.S;greedy_x];
            end

            for i=1:ceil(n*n_S/10)
                greedy_x=greedy_sol_Du(Data.prob);
                Data.S=[Data.S;greedy_x];
            end
            

            Iteration=Dim*150;


            % miso('datainput_dp',Iteration, 'rbf_c', [], 'slhd', 'cp4',[],Data); %SODA-ADM
            miso('datainput_dp_speed',Iteration, 'rbf_c', [], 'slhd', 'soda_adm_fu',[],Data); %the new SODA-ADM

            % [xbest, fbest] = miso('datainput_dp',Iteration, 'rbf_c', [], 'slhd', 'cp6',[],Data); %SODA-ADM-DP
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
sol_tv=ones(1,prob.n*(prob.n_S-1));
assign_next=[randi(prob.n_M(1),1,prob.n)];

for s=1:prob.n_S
    D=sortrows(D);
    stage_seq=D(:,2);
    temp_sol=[];
    for i=prob.n:-1:2
        ind=find(stage_seq==i);
        stage_seq(ind)=[];
        temp_sol=[ind,temp_sol];
    end
    sol_seq=[sol_seq temp_sol];
    assign=assign_next;
    sol_assign=[sol_assign,assign];
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

for s=1:prob.n_S
    D=sortrows(D);
    stage_seq=D(:,2);
    temp_sol=[];
    for i=prob.n:-1:2
        ind=find(stage_seq==i);
        stage_seq(ind)=[];
        temp_sol=[ind,temp_sol];
    end
    sol_seq=[sol_seq temp_sol];
    assign=assign_next;
    sol_assign=[sol_assign,assign];
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

for s=1:prob.n_S
    D=sortrows(D);
    stage_seq=D(:,2);
    temp_sol=[];
    for i=prob.n:-1:2
        ind=find(stage_seq==i);
        stage_seq(ind)=[];
        temp_sol=[ind,temp_sol];
    end
    sol_seq=[sol_seq temp_sol];
    assign=assign_next;
    sol_assign=[sol_assign,assign];
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