clc;clear;
n_array=[20,40,60,80,100,120];
repeat=15;
n_S_array=[2,4,8];
n_V=3;
n_T=3;
for rr=1:repeat
    for n_S=n_S_array
        for n=n_array
            prob=load (strcat('problems\prob_',num2str(n),'_',num2str(n_S),'.mat'));


            prob.n_M=prob.machine(1:n_S);
            n_M=prob.n_M;

            int_dim=(n-1)*n_S+n*n_S+sum(n_M)+n*(n_S-1);% seq, assign, pr_speed, tr_speed
            con_dim=n*(n_S-1);
            Dim=int_dim+con_dim;

            Data.number_startpoints=2*(Dim+1);
            Data.dim=Dim;
            Data.xlow=[ones(1,int_dim),zeros(1,con_dim)];

            Data.xup=[repmat(2:n,1,n_S)];
            for i= 1:n_S
                Data.xup=[Data.xup,repelem(n_M(i),n)];
            end
            Data.xup=[Data.xup,repelem(n_V,sum(n_M)),repelem(n_T,n*(n_S-1))];
            Data.xup=[Data.xup,ones(1,con_dim)];


            Data.integer=[1:int_dim]; %indices of integer variables
            Data.category=[(n-1)*n_S+1:(n-1)*n_S+n*n_S]; %indices of category integer variables
            Data.continuous=[int_dim+1:Dim]; %indices of continuous variables
            prob.int_dim=int_dim;
            Data.prob=prob;


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


            Iteration=min(Dim*100,3000);

            % miso('datainput_dp',Iteration, 'rbf_c', [], 'slhd', 'cp4',[],Data); %SODA-ADM
            miso('datainput_rangedp',Iteration, 'rbf_c', [], 'slhd', 'soda_adm_fu',[],Data); %the new SODA-ADM

            % [xbest, fbest] = miso('datainput_dp',Iteration, 'rbf_c', [], 'slhd', 'cp6',[],Data); %SODA-ADM-DP

        end
    end
end