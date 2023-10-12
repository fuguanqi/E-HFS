function Data =  cprbf2(Data,PP,cp)   %扰动半径4组

%% optimization iterations
% while Data.m < Data.maxeval  %do until budget of function evaluations exhausted 
    %compute RBF parameters
    if strcmp(Data.surrogate, 'rbf_c') %cubic RBF
        rbf_flag = 'cub';
        [lambda, gamma] = rbf_params(Data, rbf_flag);
    elseif strcmp(Data.surrogate, 'rbf_l') %linear RBF
        rbf_flag = 'lin';
        [lambda, gamma] = rbf_params(Data, rbf_flag);
    elseif strcmp(Data.surrogate, 'rbf_t') %thin plate spline RBF
        rbf_flag = 'tps';
        [lambda, gamma] = rbf_params(Data, rbf_flag);
    else
        error('rbf type not defined')
    end
    
    cp.iterctr = cp.iterctr + 1; %increment iteration counter
    fprintf('\n Iteration: %d \n',cp.iterctr); %print some user information
    fprintf('\n Number of function evaluations: %d \n', Data.m);
    fprintf('\n Best function value so far: %d \n', Data.fbest);
    


    % Compute perturbation probability
%     pert_p = min(20/Data.dim,1)*(1-(log(Data.m-2*(Data.dim+1)+1)/log(Data.maxeval-2*(Data.dim+1))));
    %create candidate points by perturbing best point found
    CandPoint1 = repmat(Data.xbest,cp.ncand,1);     %大扰动
    for ii =1:cp.ncand
%         r=rand(1,Data.dim);
%         ar = r<PP; %indices of variables to be perturbed
%         if ~(any(ar)) %if no variable is to be perturbed, randomly select one
%             r = randperm(Data.dim);
%             ar(r(1))=1;
%         end
        for jj =1:Data.dim %go through all dimensions and perturb variable values as necessary
%             if ar(jj)==1
            if rand(1)<= PP(jj)     %产生一个随机数与扰动概率比较，决定是否进行扰动
                if ismember(jj,Data.integer) %integer perturbation has to be at least 1 unit 
                    rr=randn(1);
                    s_std=sign(rr)*max(1,abs(cp.sigma_stdev(1,jj)*rr));
                else
                    s_std= cp.sigma_stdev(1,jj)*randn(1);
                end
                CandPoint1(ii,jj) = CandPoint1(ii,jj) +s_std;
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
        CandPoint2 = repmat(Data.xbest,cp.ncand,1);      %中扰动
    for ii =1:cp.ncand
        for jj =1:Data.dim %go through all dimensions and perturb variable values as necessary
%             if ar(jj)==1
            if rand(1)<= PP(jj)     %产生一个随机数与扰动概率比较，决定是否进行扰动
                if ismember(jj,Data.integer) %integer perturbation has to be at least 1 unit 
                    rr=randn(1);
                    s_std=sign(rr)*max(1,abs(cp.sigma_stdev(2,jj)*rr));
                else
                    s_std= cp.sigma_stdev(2,jj)*randn(1);
                end
                CandPoint2(ii,jj) = CandPoint2(ii,jj) +s_std;
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
            CandPoint3 = repmat(Data.xbest,cp.ncand,1);      %小扰动
    for ii =1:cp.ncand
        for jj =1:Data.dim %go through all dimensions and perturb variable values as necessary
%             if ar(jj)==1
            if rand(1)<= PP(jj)     %产生一个随机数与扰动概率比较，决定是否进行扰动
                if ismember(jj,Data.integer) %integer perturbation has to be at least 1 unit 
                    rr=randn(1);
                    s_std=sign(rr)*max(1,abs(cp.sigma_stdev(3,jj)*rr));
                else
                    s_std= cp.sigma_stdev(3,jj)*randn(1);
                end
                CandPoint3(ii,jj) = CandPoint3(ii,jj) +s_std;
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
         CandPoint4 = repmat(Data.xbest,cp.ncand,1);      %中扰动
    for ii =1:cp.ncand
        for jj =1:Data.dim %go through all dimensions and perturb variable values as necessary
%             if ar(jj)==1
            if rand(1)<= PP(jj)     %产生一个随机数与扰动概率比较，决定是否进行扰动
                if ismember(jj,Data.integer) %integer perturbation has to be at least 1 unit 
                    rr=randn(1);
                    s_std=sign(rr)*max(1,abs(cp.sigma_stdev(4,jj)*rr));
                else
                    s_std= cp.sigma_stdev(4,jj)*randn(1);
                end
                CandPoint4(ii,jj) = CandPoint4(ii,jj) +s_std;
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
%       CandPoint5 = repmat(Data.xbest,cp.ncand,1);      %中扰动
%     for ii =1:cp.ncand
%         for jj =1:Data.dim %go through all dimensions and perturb variable values as necessary
% %             if ar(jj)==1
%             if rand(1)<= PP(jj)     %产生一个随机数与扰动概率比较，决定是否进行扰动
%                 if ismember(jj,Data.integer) %integer perturbation has to be at least 1 unit 
%                     rr=randn(1);
%                     s_std=sign(rr)*max(1,abs(cp.sigma_stdev(5)*rr));
%                 else
%                     s_std= cp.sigma_stdev(5)*randn(1);
%                 end
%                 CandPoint5(ii,jj) = CandPoint5(ii,jj) +s_std;
%                 if CandPoint5(ii,jj) < Data.xlow(jj)
%                     CandPoint5(ii,jj) = Data.xlow(jj)+ (Data.xlow(jj)-CandPoint5(ii,jj));
%                     if CandPoint5(ii,jj) >Data.xup(jj)
%                         CandPoint5(ii,jj) = Data.xlow(jj);
%                     end
%                 elseif CandPoint5(ii,jj) > Data.xup(jj)
%                     CandPoint5(ii,jj) = Data.xup(jj)- (CandPoint5(ii,jj)-Data.xup(jj));
%                     if CandPoint5(ii,jj) <Data.xlow(jj)
%                         CandPoint5(ii,jj) = Data.xup(jj);
%                     end
%                 end
%             end
%         end
%     end
%        CandPoint6 = repmat(Data.xbest,cp.ncand,1);      %中扰动
%     for ii =1:cp.ncand
%         for jj =1:Data.dim %go through all dimensions and perturb variable values as necessary
% %             if ar(jj)==1
%             if rand(1)<= PP(jj)     %产生一个随机数与扰动概率比较，决定是否进行扰动
%                 if ismember(jj,Data.integer) %integer perturbation has to be at least 1 unit 
%                     rr=randn(1);
%                     s_std=sign(rr)*max(1,abs(cp.sigma_stdev(6)*rr));
%                 else
%                     s_std= cp.sigma_stdev(6)*randn(1);
%                 end
%                 CandPoint6(ii,jj) = CandPoint6(ii,jj) +s_std;
%                 if CandPoint6(ii,jj) < Data.xlow(jj)
%                     CandPoint6(ii,jj) = Data.xlow(jj)+ (Data.xlow(jj)-CandPoint6(ii,jj));
%                     if CandPoint6(ii,jj) >Data.xup(jj)
%                         CandPoint6(ii,jj) = Data.xlow(jj);
%                     end
%                 elseif CandPoint6(ii,jj) > Data.xup(jj)
%                     CandPoint6(ii,jj) = Data.xup(jj)- (CandPoint6(ii,jj)-Data.xup(jj));
%                     if CandPoint6(ii,jj) <Data.xlow(jj)
%                         CandPoint6(ii,jj) = Data.xup(jj);
%                     end
%                 end
%             end
%         end
%     end
    %create second group of candidates by uniformly selecting sample points
%     CandPoint4=repmat(Data.xlow, cp.ncand,1) + repmat(Data.xup-Data.xlow,cp.ncand,1).*rand(cp.ncand,Data.dim); 
    CandPoint=[CandPoint1;CandPoint2;CandPoint3;CandPoint4];
    CandPoint(:,Data.integer)=round(CandPoint(:,Data.integer)); %round integer variable values
    
    xnew = compute_scores_rbf(Data,CandPoint, lambda, gamma, rbf_flag ); %select the best candidate
    clear CandPoint;
    fevalt = tic; %start timer for function evaluation
    fnew = feval(Data.objfunction,xnew); %new function value                    %要改
    timer = toc(fevalt); %stop timer for function evaluation
    Data.m=Data.m+1; %update the number of function evaluations
    m=Data.m;
    Data.S(Data.m,:)=xnew; %update sample site matrix with new point
    Data.Y(Data.m,1)=fnew; %update vector with function values
    Data.T(Data.m,1) = timer; %update vector with evaluation times

    if fnew < Data.fbest %update best point found so far if necessary
        if (Data.fbest - fnew) > (1e-3)*abs(Data.fbest)
            % "significant" improvement
            cp.failctr = 0; %update fai and success counters
            cp.succctr=cp.succctr+1;                   %提升成功的次数加一
            improvement=true;
        else
            %no "significant" improvement
            cp.failctr = cp.failctr + 1; %update fai and success counters       %提升失败的次数加一
            cp.succctr=0;
            improvement=false;
        end  
        xbest_old=Data.xbest;
        Data.xbest = xnew; %best point found so far                                      %更新了xbest和fbest
        Data.fbest = fnew; %best objective function value found so far
    else
        
        cp.failctr = cp.failctr + 1; %update fai and success counters
        cp.succctr=0;
        improvement=false;
    end
    
   
%     if improvement %if improvement found, do local candidate point search on continuous variables only if integer variables had changed
%         if ~all(Data.xbest(Data.integer) - xbest_old(Data.integer)==0) && length(Data.continuous)~=0 %integer variables changed
%             if strcmp(Data.surrogate, 'rbf_c') %cubic RBF
%                 rbf_flag = 'cub';
%                 [lambda, gamma] = rbf_params(Data, rbf_flag);
%             elseif strcmp(Data.surrogate, 'rbf_l') %linear RBF
%                 rbf_flag = 'lin';
%                 [lambda, gamma] = rbf_params(Data, rbf_flag);
%             elseif strcmp(Data.surrogate, 'rbf_t') %thin plate spline RBF
%                 rbf_flag = 'tps';
%                 [lambda, gamma] = rbf_params(Data, rbf_flag);
%             else
%                 error('rbf type not defined')
%             end      
%            
%             cont_search_imp=false;
%             kl=1;
%             while kl <= 5 && Data.m < Data.maxeval
%                 % Perturbation probability
%                 pert_prob = min(20/Data.dim,1)*(1-(log(Data.m-2*(Data.dim+1)+1)/log(Data.maxeval-2*(Data.dim+1))));
%                 %create candidate points
%                 CandPoint = repmat(Data.xbest,cp.ncand,1); 
%                 for ii =1:cp.ncand
%                     r=rand(1,length(Data.continuous));
%                     ar=r<pert_prob;
%                     if ~(any(ar))
%                         r = randperm(length(Data.continuous));
%                         ar(r(1))=1;
%                     end
%                     for jj =1:length(Data.continuous)
%                         if ar(jj)==1
%                             sig = cp.sigma_stdev(3,jj)*randn(1);
%                             CandPoint(ii,Data.continuous(jj)) = CandPoint(ii,Data.continuous(jj)) +sig;
%                             if CandPoint(ii,Data.continuous(jj)) < Data.xlow(Data.continuous(jj))
%                                 CandPoint(ii,Data.continuous(jj)) = Data.xlow(Data.continuous(jj))+ (Data.xlow(Data.continuous(jj))-CandPoint(ii,Data.continuous(jj)));
%                                 if CandPoint(ii,Data.continuous(jj)) >Data.xup(Data.continuous(jj))
%                                     CandPoint(ii,Data.continuous(jj)) = Data.xlow(Data.continuous(jj));
%                                 end
%                             elseif CandPoint(ii,Data.continuous(jj)) > Data.xup(Data.continuous(jj))
%                                 CandPoint(ii,Data.continuous(jj)) = Data.xup(Data.continuous(jj))- (CandPoint(ii,Data.continuous(jj))-Data.xup(Data.continuous(jj)));
%                                 if CandPoint(ii,Data.continuous(jj)) <Data.xlow(Data.continuous(jj))
%                                     CandPoint(ii,Data.continuous(jj)) = Data.xup(Data.continuous(jj));
%                                 end
%                             end
%                         end
%                     end
%                 end
%                 xnew = compute_scores_rbf(Data,CandPoint, lambda, gamma, rbf_flag ); %select the best candidate
%                 fevalt = tic; % start timer for function evaluation 
%                 Fselected = feval(Data.objfunction,xnew);  %new function value           %要改
%                 timer = toc(fevalt); %stop timer for function evaluation
%                 Data.m=Data.m+1; %update the number of function evaluations
%                 Data.T(Data.m,1) = timer; %update the vector with evaluation times
%                 Data.S(Data.m,:)=xnew; %update the sample site matrix
%                 Data.Y(Data.m,1)=Fselected; %update the vector with function values
%                 if (Data.fbest - Fselected) > (1e-3)*abs(Data.fbest) %update best solution found if necessary
%                     cp.failctr = 0;
%                     cont_search_imp=true;
%                     Data.xbest = xnew; %best point found so far
%                     Data.fbest = Fselected; %best value found so far
%                 end
%                 if kl < 5 %update RBF parameters
%                     if strcmp(Data.surrogate, 'rbf_c') %cubic RBF
%                         rbf_flag = 'cub';
%                         [lambda, gamma] = rbf_params(Data, rbf_flag);
%                     elseif strcmp(Data.surrogate, 'rbf_l') %linear RBF
%                         rbf_flag = 'lin';
%                         [lambda, gamma] = rbf_params(Data, rbf_flag);
%                     elseif strcmp(Data.surrogate, 'rbf_t') %thin plate spline RBF
%                         rbf_flag = 'tps';
%                         [lambda, gamma] = rbf_params(Data, rbf_flag);
%                     else
%                         error('rbf type not defined')
%                     end  
%                     kl=kl+1;
%                 elseif kl==5 && cont_search_imp
%                     if strcmp(Data.surrogate, 'rbf_c') %cubic RBF
%                         rbf_flag = 'cub';
%                         [lambda, gamma] = rbf_params(Data, rbf_flag);
%                     elseif strcmp(Data.surrogate, 'rbf_l') %linear RBF
%                         rbf_flag = 'lin';
%                         [lambda, gamma] = rbf_params(Data, rbf_flag);
%                     elseif strcmp(Data.surrogate, 'rbf_t') %thin plate spline RBF
%                         rbf_flag = 'tps';
%                         [lambda, gamma] = rbf_params(Data, rbf_flag);
%                     else
%                         error('rbf type not defined')
%                     end  
%                     
%                     kl=1;
%                     cont_search_imp=false;
%                 else
%                     kl=kl+1;
%                 end
%             end
%         end
%     end
    
    %check if perturbation range needs to be decreased
%     cp.shrinkflag = 1;      
%     if cp.failctr >= cp.failtolerance %check how many consecutive failed improvement trials
%         if cp.failctr >= cp.failtolerance && cp.shrinkctr >= cp.maxshrinkparam
%             cp.shrinkflag = 0;
%             disp('Stopped reducing perturbation range because the maximum number of reduction has been reached.');
%         end
%         cp.failctr = 0;
% 
%         if cp.shrinkflag == 1 % decrease perturbation range
%             cp.shrinkctr = cp.shrinkctr + 1;
%             cp.sigma_stdev =  max(cp.sigma_min,cp.sigma_stdev/1.5);
%             disp('Reducing the perturbation range to half!');
%         end
%     end
%     if cp.succctr>=cp.succtolerance %check if number of consecutive improvements is large enough
%         cp.sigma_stdev=2*cp.sigma_stdev;%increase search radius
%         cp.succctr=0;
%     end

    
        if improvement==false
          r_flag=0;
        else
          r_flag=1;
        end
%未改进则缩小          
%         if r_flag==0 && cp.failctr >= cp.failtolerance  %当没找到更优解，且失败超出阈值，则缩小半径
%         if ~isequal(cp.sigma_min , cp.sigma_stdev)
%            cp.sigma_stdev=max(cp.sigma_min,cp.sigma_stdev/2);
%            cp.failctr=0;
%            cp.count1=0;
%         elseif isequal(cp.sigma_min , cp.sigma_stdev) && cp.count1~=3
%            cp.count1=cp.count1+1;
%         elseif isequal(cp.sigma_min , cp.sigma_stdev) && cp.count1==3
%            cp.sigma_stdev=min(cp.sigma_max,10*cp.sigma_stdev);  %8这个参数可以调
%            cp.count1=0;
%            cp.failctr=0;
%         end
%     end
%     if r_flag==1 && cp.succctr >= cp.succtolerance  %当找到更优解，且成功超出阈值，则扩大半径
%         cp.sigma_stdev=max(cp.sigma_min,2*cp.sigma_stdev);
%         cp.succctr=0;
%         cp.count1=0;
%     end
    
%  %未改进则扩大     
%     if r_flag==0 && cp.failctr >= cp.failtolerance  %当没找到更优解，且失败的次数超出阈值，则扩大半径
%         if ~isequal(cp.sigma_max , cp.sigma_stdev)
%            cp.sigma_stdev=min(cp.sigma_max,cp.sigma_stdev*2);  %扩大半径
%            cp.failctr=0;
%            cp.count1=0;
%         elseif isequal(cp.sigma_max , cp.sigma_stdev) && cp.count1~=3
%            cp.count1=cp.count1+1;
%         elseif isequal(cp.sigma_max , cp.sigma_stdev) && cp.count1==3
%            cp.sigma_stdev=max(cp.sigma_min,cp.sigma_stdev/8);  %缩小半径   8这个参数可以调
%            cp.count1=0;
%            cp.failctr=0;
%         end
%     end
%     if r_flag==1 && cp.succctr >= cp.succtolerance  %当找到更优解，且成功次数超出阈值，则缩小半径
%         cp.sigma_stdev=max(cp.sigma_min,0.5*cp.sigma_stdev);
%         cp.succctr=0;
%         cp.count1=0;
%     end       
% end
% %未改进则扩大 （版本三）    
    if r_flag==0 && cp.failctr >= cp.failtolerance  %当没找到更优解，且失败的次数超出阈值，则扩大半径
        if ~isequal(cp.sigma_max , cp.sigma_stdev) && ~isequal(cp.sigma_min , cp.sigma_stdev)  %若半径不是最大也不是最小，则缩小半径
           cp.sigma_stdev=max(cp.sigma_min,cp.sigma_stdev./2);  %缩小半径,直至最小
           cp.failctr=0;
        elseif isequal(cp.sigma_min , cp.sigma_stdev) 
            cp.sigma_stdev=min(cp.sigma_max,cp.sigma_stdev.*8);
            cp.failctr=0;
        elseif isequal(cp.sigma_max , cp.sigma_stdev) 
           cp.sigma_stdev=max(cp.sigma_min,cp.sigma_stdev./2);  %缩小半径   
           cp.failctr=0;
        end
    end
    if r_flag==1 && cp.succctr >= cp.succtolerance  %当找到更优解，且成功次数超出阈值，则缩小半径
        cp.sigma_stdev=max(cp.sigma_min,0.5*cp.sigma_stdev);
        cp.succctr=0;
    end              

    
%     Data.S(Data.m,:)=Data.xbest;
%     Data.x_p=Data.xbest;
%     Data.f_p=Data.fbest;
    Data.cp=cp;
    
% end  %while

end%function