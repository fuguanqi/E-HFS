function Data =  soda_adm_fu(Data)%扰动半径没有上限

%% parameters and initialization for coordinate perturbation
cp.ncand= min(10*Data.dim,5000); % number of candidate points
cp.xrange = Data.xup - Data.xlow; %variable range in every dimension
% cp.minxrange = min(cp.xrange); %smallest variable range
cp.minxrange = median(cp.xrange);
% cp.sigma_stdev_default = cp.minxrange;    %perturbation range for generating candidate points
% cp.sigma_stdev = cp.sigma_stdev_default;     % current perturbation range
cp.sigma_stdev=ones(4,Data.dim);
cp.sigma_min=ones(4,Data.dim);
cp.sigma_max=ones(4,Data.dim);
%%----缩放1----%%%%
cp.sigma_stdev(1,:) = 0.6*cp.xrange; %160
cp.sigma_stdev(2,:) = 0.1*cp.xrange; %20
cp.sigma_stdev(3,:) = 0.01*cp.xrange; %2
cp.sigma_stdev(4,:) = 0.0025*cp.xrange; %0.5
cp.sigma_min(1,:)=0.2*(1/2)^6*cp.xrange; %0.015625   %扰动半径的下界，即阈值
cp.sigma_min(2,:)=0.2*(1/2)^6*cp.xrange;
cp.sigma_min(3,:)=0.2*(1/2)^6*cp.xrange;
cp.sigma_min(4,:)=0.2*(1/2)^6*cp.xrange;
cp.sigma_max=cp.sigma_stdev;   %largest allowed perturbation range


cp.maxshrinkparam = 5; % maximal number of perturbation range reduction
cp.failtolerance = round(sqrt(Data.dim)); %threshold for consecutively failed improvement trials
cp.succtolerance =3; %threshold for consecutively successful improvement trials
cp.iterctr = 0; % number of iterations
cp.shrinkctr = 0; % number of times perturbation range was reduced
cp.failctr = 0; % number of consecutive unsuccessful iterations
cp.succctr=0; % number of consecutive successful iterations
cp.count1=0;
cp.weightpattern=[0.3,0.5,0.8,0.95]; %weight pattern for scoring candidate points
fevalt = tic;
%% optimization iterations
% while Data.m <  (Data.maxeval-3*Data.dim)  %do until budget of function evaluations exhausted 0.5* Data.maxeval
while Data.m < Data.maxeval  %do until budget of function evaluations exhausted
    %compute RBF parameters
    if strcmp(Data.surrogate, 'rbf_c') %cubic RBF
        rbf_flag = 'cub';
        [lambda, gamma] = rbf_paramsD(Data, rbf_flag);
    elseif strcmp(Data.surrogate, 'rbf_l') %linear RBF
        rbf_flag = 'lin';
        [lambda, gamma] = rbf_paramsD(Data, rbf_flag);
    elseif strcmp(Data.surrogate, 'rbf_t') %thin plate spline RBF
        rbf_flag = 'tps';
        [lambda, gamma] = rbf_paramsD(Data, rbf_flag);
    else
        error('rbf type not defined')
    end

    cp.iterctr = cp.iterctr + 1; %increment iteration counter
    fprintf('\n Iteration: %d \n',cp.iterctr); %print some user information
    fprintf('\n Number of function evaluations: %d \n', Data.m);
    fprintf('\n Best function value so far: %d \n', Data.fbest);

    %select weight for score computation
    mw=mod(cp.iterctr,length(cp.weightpattern));
    if mw==0
        w_r=cp.weightpattern(end);
    else
        w_r=cp.weightpattern(mw);
    end

    % Compute perturbation probability
    %%%%%原始概率%%%%%
    pert_p = min(20/Data.dim,1)*(1-(log(Data.m-2*(Data.dim+1)+1)/log(Data.maxeval-2*(Data.dim+1))));
    %create candidate points by perturbing best point found
    %         if  pert_p <0.6
    %            pert_p= pert_p+0.2;
    %         end
    CandPoint1 = repmat(Data.xbest,cp.ncand,1);
    for ii =1:cp.ncand
        r=rand(1,Data.dim);
        ar = r<pert_p; %indices of variables to be perturbed
        if ~(any(ar)) %if no variable is to be perturbed, randomly select one
            r = randperm(Data.dim);
            ar(r(1))=1;
        end
        for jj =1:Data.dim %go through all dimensions and perturb variable values as necessary
            if ar(jj)==1
                if jj>=Data.integer(1) && jj<=Data.integer(end)
                    %                 if ismember(jj,Data.integer) %integer perturbation has to be at least 1 unit
                    rr=randn(1);
                    s_std=sign(rr)*max(1,abs(cp.sigma_stdev(1,jj)*rr));
                else
                    s_std= cp.sigma_stdev(1,jj)*randn(1);
                end
                CandPoint1(ii,jj) = CandPoint1(ii,jj) +s_std;

                    if ismember(jj,Data.category)
                        if CandPoint1(ii,jj) < Data.xlow(jj)
                            CandPoint1(ii,jj) = Data.xup(jj)- (Data.xlow(jj)-CandPoint1(ii,jj));
                            if CandPoint1(ii,jj) < Data.xlow(jj)
                                CandPoint1(ii,jj) = Data.xup(jj);
                            end
                        elseif CandPoint1(ii,jj) > Data.xup(jj)
                            CandPoint1(ii,jj) = Data.xlow(jj)+ (CandPoint1(ii,jj)-Data.xup(jj));
                            if CandPoint1(ii,jj) > Data.xup(jj)
                                CandPoint1(ii,jj) = Data.xlow(jj);
                            end
                        end
                    else
                        if CandPoint1(ii,jj) < Data.xlow(jj)
                            CandPoint1(ii,jj) = Data.xlow(jj)+ (Data.xlow(jj)-CandPoint1(ii,jj));
                            if CandPoint1(ii,jj) >Data.xup(jj)
                                CandPoint1(ii,jj) = Data.xlow(jj);
                            end
                        elseif CandPoint1(ii,jj) > Data.xup(jj)
                            CandPoint1(ii,jj) = Data.xup(jj)- (CandPoint1(ii,jj)-Data.xup(jj));
                            if CandPoint1(ii,jj) <Data.xlow(jj)
                                CandPoint1(ii,jj) = Data.xup(jj);
                            end
                        end
                    end

            end
        end
    end
    CandPoint2 = repmat(Data.xbest,cp.ncand,1);
    for ii =1:cp.ncand
        r=rand(1,Data.dim);
        ar = r<pert_p; %indices of variables to be perturbed
        if ~(any(ar)) %if no variable is to be perturbed, randomly select one
            r = randperm(Data.dim);
            ar(r(1))=1;
        end
        for jj =1:Data.dim %go through all dimensions and perturb variable values as necessary
            if ar(jj)==1
                if jj>=Data.integer(1) && jj<=Data.integer(end)
                    %                 if ismember(jj,Data.integer) %integer perturbation has to be at least 1 unit
                    rr=randn(1);
                    s_std=sign(rr)*max(1,abs(cp.sigma_stdev(2,jj)*rr));
                else
                    s_std= cp.sigma_stdev(2,jj)*randn(1);
                end
                CandPoint2(ii,jj) = CandPoint2(ii,jj) +s_std;
                if ismember(jj,Data.category)
                    if CandPoint2(ii,jj) < Data.xlow(jj)
                        CandPoint2(ii,jj) = Data.xup(jj)- (Data.xlow(jj)-CandPoint2(ii,jj));
                        if CandPoint2(ii,jj) < Data.xlow(jj)
                            CandPoint2(ii,jj) = Data.xup(jj);
                        end
                    elseif CandPoint2(ii,jj) > Data.xup(jj)
                        CandPoint2(ii,jj) = Data.xlow(jj)+ (CandPoint2(ii,jj)-Data.xup(jj));
                        if CandPoint2(ii,jj) > Data.xup(jj)
                            CandPoint2(ii,jj) = Data.xlow(jj);
                        end
                    end
                else
                    if CandPoint2(ii,jj) < Data.xlow(jj)
                        CandPoint2(ii,jj) = Data.xlow(jj)+ (Data.xlow(jj)-CandPoint2(ii,jj));
                        if CandPoint2(ii,jj) >Data.xup(jj)
                            CandPoint2(ii,jj) = Data.xlow(jj);
                        end
                    elseif CandPoint2(ii,jj) > Data.xup(jj)
                        CandPoint2(ii,jj) = Data.xup(jj)- (CandPoint2(ii,jj)-Data.xup(jj));
                        if CandPoint2(ii,jj) <Data.xlow(jj)
                            CandPoint2(ii,jj) = Data.xup(jj);
                        end
                    end
                end
            end
        end
    end
    CandPoint3 = repmat(Data.xbest,cp.ncand,1);
    for ii =1:cp.ncand
        r=rand(1,Data.dim);
        ar = r<pert_p; %indices of variables to be perturbed
        if ~(any(ar)) %if no variable is to be perturbed, randomly select one
            r = randperm(Data.dim);
            ar(r(1))=1;
        end
        for jj =1:Data.dim %go through all dimensions and perturb variable values as necessary
            if ar(jj)==1
                if jj>=Data.integer(1) && jj<=Data.integer(end)
                    %                 if ismember(jj,Data.integer) %integer perturbation has to be at least 1 unit
                    rr=randn(1);
                    s_std=sign(rr)*max(1,abs(cp.sigma_stdev(3,jj)*rr));
                else
                    s_std= cp.sigma_stdev(3,jj)*randn(1);
                end
                CandPoint3(ii,jj) = CandPoint3(ii,jj) +s_std;
                if ismember(jj,Data.category)
                    if CandPoint3(ii,jj) < Data.xlow(jj)
                        CandPoint3(ii,jj) = Data.xup(jj)- (Data.xlow(jj)-CandPoint3(ii,jj));
                        if CandPoint3(ii,jj) < Data.xlow(jj)
                            CandPoint3(ii,jj) = Data.xup(jj);
                        end
                    elseif CandPoint3(ii,jj) > Data.xup(jj)
                        CandPoint3(ii,jj) = Data.xlow(jj)+ (CandPoint3(ii,jj)-Data.xup(jj));
                        if CandPoint3(ii,jj) > Data.xup(jj)
                            CandPoint3(ii,jj) = Data.xlow(jj);
                        end
                    end
                else
                    if CandPoint3(ii,jj) < Data.xlow(jj)
                        CandPoint3(ii,jj) = Data.xlow(jj)+ (Data.xlow(jj)-CandPoint3(ii,jj));
                        if CandPoint3(ii,jj) >Data.xup(jj)
                            CandPoint3(ii,jj) = Data.xlow(jj);
                        end
                    elseif CandPoint3(ii,jj) > Data.xup(jj)
                        CandPoint3(ii,jj) = Data.xup(jj)- (CandPoint3(ii,jj)-Data.xup(jj));
                        if CandPoint3(ii,jj) <Data.xlow(jj)
                            CandPoint3(ii,jj) = Data.xup(jj);
                        end
                    end
                end
            end
        end
    end
    CandPoint4 = repmat(Data.xbest,cp.ncand,1);
    for ii =1:cp.ncand
        r=rand(1,Data.dim);
        ar = r<pert_p; %indices of variables to be perturbed
        if ~(any(ar)) %if no variable is to be perturbed, randomly select one
            r = randperm(Data.dim);
            ar(r(1))=1;
        end
        for jj =1:Data.dim %go through all dimensions and perturb variable values as necessary
            if ar(jj)==1
                if jj>=Data.integer(1) && jj<=Data.integer(end)
                    %                 if ismember(jj,Data.integer) %integer perturbation has to be at least 1 unit
                    rr=randn(1);
                    s_std=sign(rr)*max(1,abs(cp.sigma_stdev(4,jj)*rr));
                else
                    s_std= cp.sigma_stdev(4,jj)*randn(1);
                end
                CandPoint4(ii,jj) = CandPoint4(ii,jj) +s_std;
                if ismember(jj,Data.category)
                    if CandPoint4(ii,jj) < Data.xlow(jj)
                        CandPoint4(ii,jj) = Data.xup(jj)- (Data.xlow(jj)-CandPoint4(ii,jj));
                        if CandPoint4(ii,jj) < Data.xlow(jj)
                            CandPoint4(ii,jj) = Data.xup(jj);
                        end
                    elseif CandPoint4(ii,jj) > Data.xup(jj)
                        CandPoint4(ii,jj) = Data.xlow(jj)+ (CandPoint4(ii,jj)-Data.xup(jj));
                        if CandPoint4(ii,jj) > Data.xup(jj)
                            CandPoint4(ii,jj) = Data.xlow(jj);
                        end
                    end
                else
                    if CandPoint4(ii,jj) < Data.xlow(jj)
                        CandPoint4(ii,jj) = Data.xlow(jj)+ (Data.xlow(jj)-CandPoint4(ii,jj));
                        if CandPoint4(ii,jj) >Data.xup(jj)
                            CandPoint4(ii,jj) = Data.xlow(jj);
                        end
                    elseif CandPoint4(ii,jj) > Data.xup(jj)
                        CandPoint4(ii,jj) = Data.xup(jj)- (CandPoint4(ii,jj)-Data.xup(jj));
                        if CandPoint4(ii,jj) <Data.xlow(jj)
                            CandPoint4(ii,jj) = Data.xup(jj);
                        end
                    end
                end
            end
        end
    end
    %create second group of candidates by uniformly selecting sample points
    %     CandPoint2=repmat(Data.xlow, cp.ncand,1) + repmat(Data.xup-Data.xlow,cp.ncand,1).*rand(cp.ncand,Data.dim);
    CandPoint=[CandPoint1;CandPoint2;CandPoint3;CandPoint4];
    CandPoint(:,Data.integer)=round(CandPoint(:,Data.integer)); %round integer variable values

    if Data.m < 0.8* Data.maxeval
        xnew = compute_scores(Data,CandPoint,w_r, lambda, gamma, rbf_flag); %select the best candidate
    else
        xnew = compute_scores_rbf(Data,CandPoint, lambda, gamma, rbf_flag );
    end

    clear CandPoint;

    fnew = feval(Data.objfunction,xnew,Data.prob); %new function value                    %要改

    Data.m=Data.m+1; %update the number of function evaluations
    m=Data.m;
    Data.S(Data.m,:)=xnew; %update sample site matrix with new point
    Data.Y(Data.m,1)=fnew; %update vector with function values
    Data.cur_best(Data.m,1)=min(Data.cur_best(Data.m-1,1),fnew);



    if fnew < Data.fbest %update best point found so far if necessary
        if (Data.fbest - fnew) > (1e-3)*abs(Data.fbest)
            % "significant" improvement
            cp.failctr = 0; %update fai and success counters
            cp.count1=0;
            cp.succctr=cp.succctr+1;
            improvement=true;
        else
            %no "significant" improvement
            cp.failctr = cp.failctr + 1; %update fai and success counters
            cp.succctr=0;
            improvement=false;
        end
        xbest_old=Data.xbest;
        Data.xbest = xnew; %best point found so far
        Data.fbest = fnew; %best objective function value found so far
    else
        cp.failctr = cp.failctr + 1; %update fai and success counters
        cp.succctr=0;
        improvement=false;
    end


    if improvement==false
        r_flag=0;
    else
        r_flag=1;
    end

    %     if r_flag==0 && cp.failctr >= cp.failtolerance  %当没找到更优解，且失败超出阈值，则缩小半径
    %         for gg=1:4
    %         if ~isequal(cp.sigma_min(gg,:) , cp.sigma_stdev(gg,:))
    %            cp.sigma_stdev(gg,:)=max(cp.sigma_min(gg,:),cp.sigma_stdev(gg,:)./2);
    %            cp.failctr=0;
    %            cp.count1=0;
    %         elseif isequal(cp.sigma_min(gg,:) , cp.sigma_stdev(gg,:))
    %            cp.sigma_stdev(gg,:)=min(cp.sigma_max(gg,:),cp.sigma_stdev(gg,:).*100);  %8这个参数可以调
    %            cp.count1=0;
    %            cp.failctr=0;
    %         end
    %         end %for
    %     end
    %     if r_flag==1 && cp.succctr >= cp.succtolerance  %当找到更优解，且成功超出阈值，则扩大半径
    %         for gg=1:4
    %         cp.sigma_stdev(gg,:)=max(cp.sigma_min(gg,:),cp.sigma_stdev(gg,:).*2);
    %         cp.succctr=0;
    %         cp.count1=0;
    %         end
    %     end
    if r_flag==0 && cp.failctr >= cp.failtolerance  %当没找到更优解，且失败超出阈值，则缩小半径
        if ~isequal(cp.sigma_min , cp.sigma_stdev)
            cp.sigma_stdev=max(cp.sigma_min,cp.sigma_stdev./2);
            cp.failctr=0;
            cp.count1=0;
        elseif isequal(cp.sigma_min , cp.sigma_stdev)
            cp.sigma_stdev=min(cp.sigma_max,cp.sigma_stdev.*100);  %8这个参数可以调
            cp.count1=0;
            cp.failctr=0;
        end
    end
    if r_flag==1 && cp.succctr >= cp.succtolerance  %当找到更优解，且成功超出阈值，则扩大半径
        cp.sigma_stdev=max(cp.sigma_min,cp.sigma_stdev.*2);
        cp.succctr=0;
        cp.count1=0;
    end
    Data.T(Data.m,1) = toc(fevalt); %update vector with evaluation times
end

%% -------逐维扰动---------%
%% parameters and initialization for coordinate perturbation
cp.ncand= min(10*Data.dim,5000); % number of candidate points
cp.xrange = Data.xup - Data.xlow; %variable range in every dimension
% cp.minxrange = min(cp.xrange); %smallest variable range
cp.minxrange = median(cp.xrange);
% cp.sigma_stdev_default = 0.5*cp.minxrange;    %perturbation range for generating candidate points

cp.sigma_stdev(1,:) = 0.6*cp.xrange; %160
cp.sigma_stdev(2,:) = 0.2*cp.xrange; %20
cp.sigma_stdev(3,:) = 0.025*cp.xrange; %2
cp.sigma_stdev(4,:) = 0.005*cp.xrange; %0.5
cp.sigma_min(1,:)=0.6*(1/2)^4*cp.xrange; %0.015625   %扰动半径的下界，即阈值
cp.sigma_min(2,:)=0.2*(1/2)^4*cp.xrange;
cp.sigma_min(3,:)=0.025*(1/2)^4*cp.xrange;
cp.sigma_min(4,:)=0.005*(1/2)^4*cp.xrange;
cp.sigma_max=cp.sigma_stdev;                    %扰动半径上界

cp.maxshrinkparam = 5; % maximal number of perturbation range reduction
cp.failtolerance = round(sqrt(Data.dim)); %threshold for consecutively failed improvement trials
cp.succtolerance =3; %threshold for consecutively successful improvement trials
cp.iterctr = 0; % number of iterations
cp.shrinkctr = 0; % number of times perturbation range was reduced
cp.failctr = 0; % number of consecutive unsuccessful iterations
cp.succctr=0; % number of consecutive successful iterations
cp.weightpattern=[0.3,0.5,0.8,0.95]; %weight pattern for scoring candidate points

Data.S(Data.m,:)=Data.xbest;
Data.Fbest=Data.fbest;
%% -------判断阶段----------%
while Data.m < Data.maxeval
    m=Data.m;

    mw=mod(cp.iterctr,length(cp.weightpattern));
    if mw==0
        w_r=cp.weightpattern(end);
    else
        w_r=cp.weightpattern(mw);
    end
    % %-------逐维扰动----------%
    for kk=1:Data.dim
        PP=zeros(1,Data.dim);
        PP(kk)=1;
        %        SOL=cp_ME(Data,PP,cp);
        SOL=cprbf2(Data,PP,cp);
        %        SOL=cp_dafen(Data,PP,w_r,cp);
        cp=SOL.cp;
        Data.xbest=SOL.xbest;
        Data.fbest=SOL.fbest;
        %          f_p=feval(Data.objfunction,SOL.x_p);
        %            if f_p < Data.fbest
        %              Data.xbest=SOL.x_p;
        %              Data.fbest=f_p;
        %            end
        %     Data.S = [Data.S;Data.xbest]; %采样点集
        %       if length(SOL.S)==Data.S+1
        %           Data.S = [Data.S;Data.xbest];
        %       elseif length(SOL.S) > Data.S+1
        %           Data.S = SOL.S;
        %             if f_p > Data.fbest
        %                Data.S = [Data.S;Data.xbest];
        %             end
        %       end
        Data.S=SOL.S;
        Data.m = SOL.m;
        Data.Y = SOL.Y;
    end


    %------------更新数据:cp-------------%

    Data.Fbest=[Data.Fbest;Data.fbest];
    xbest= Data.xbest;
    fbest=Data.fbest;
    fprintf('\n Number of function evaluations: %d \n', Data.m);

end %while





end%function