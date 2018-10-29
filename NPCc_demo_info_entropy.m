
% function NPC_demo_info_entropy(N,rp,MOD,knots)


N=1000;
rp=5;
MOD=2;
knots=100;

parallel.defaultClusterProfile('local')

    
    rmpath(genpath('/home/hs258/Codes_Folder/HOUMAN_GITLAB_IIT/Copula_Arno'));
    addpath(genpath('/home/hs258/Codes_Folder/HOUMAN_GITLAB_IIT/Selmaan_Switch/glmnet'));
    addpath(genpath('/home/hs258/Codes_Folder/HOUMAN_GITLAB_IIT/TOOLS/ENTROPY_METHODS'));
    addpath(genpath('/home/hs258/Codes_Folder/Houman_Git/Ker_Copula'));
    
   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for WW=1%N=2.^[5:13]
        for rep=rp%(rp-1)*9+1:rp*9
                    rng('shuffle','v5uniform')
            
            switch MOD
                case 1
                    %%% model 1
                    
                    repp=mod(rep-1,9)+1;
                    ro=0.1*repp+1e-6*rand;%+0.1*rand(1)-0.05;
                    
                    vineP.margins = cell(2,1);
                    vineP.margins{1}.dist = 'norm';
                    vineP.margins{1}.theta = [0;1];
                    vineP.margins{1}.iscont = true; % Continuous margin
                    vineP.margins{2}.dist = 'norm';
                    vineP.margins{2}.theta = [0;1];
                    vineP.margins{2}.iscont = true; % Continuous margin
                    vineP.families = cell(2);
                    vineP.theta = cell(2);
                    vineP.families{1,2} = 'gaussian';
                    vineP.theta{1,2} = ro;
                    
                    X{1}=mixedvinernd(vineP,N);
                    
                case 2
                    %%% model 2
                    
                    repp=mod(rep-1,9)+1;
                    ro=0.1*repp+1e-6*rand;%+0.1*rand(1)-0.05)*2;
                    
                    vineP.margins = cell(2,1);
                    vineP.margins{1}.dist = 'norm';
                    vineP.margins{1}.theta = [0;1];
                    vineP.margins{1}.iscont = true; % Continuous margin
                    vineP.margins{2}.dist = 'norm';
                    vineP.margins{2}.theta = [0;1];
                    vineP.margins{2}.iscont = true; % Continuous margin
                    vineP.families = cell(2);
                    vineP.theta = cell(2);
                    vineP.families{1,2} = 'student';
                    vineP.theta{1,2} = [0;ro];          %%%%% ?ro
                    
                    X{2}=mixedvinernd(vineP,N);
                    %     figure;plot(X{2}(:,1),X{2}(:,2),'.')
                    
                case 3
                    %%% model 3
                    %     rep=3; %
                    %     repp=mod(rep-1,9)+1;
                    %     ro=0.1*repp;
                    repp=mod(rep-1,9)+1;
                    ro=0.1*repp+1e-6*rand;%+0.1*rand(1)-0.05)*2;
                    
                    vineP.margins = cell(2,1);
                    vineP.margins{1}.dist = 'gam';
                    vineP.margins{1}.theta = [0.1;10];
                    vineP.margins{1}.iscont = true; % Continuous margin
                    vineP.margins{2}.dist = 'gam';
                    vineP.margins{2}.theta = [0.1;10];
                    vineP.margins{2}.iscont = true; % Continuous margin
                    vineP.families = cell(2);
                    vineP.theta = cell(2);
                    vineP.families{1,2} = 'gaussian';
                    vineP.theta{1,2} = ro;
                    
                    X{3}=mixedvinernd(vineP,N);
                     
                case 4
                    %%% model 4
                    repp=mod(rep-1,9)+1;
                    ro=0.1*repp+1e-6*rand;%+0.1*rand(1)-0.05)*2;
                    
                    vineP.margins = cell(2,1);
                    vineP.margins{1}.dist = 'gam';
                    vineP.margins{1}.theta = [0.1;10];
                    vineP.margins{1}.iscont = true; % Continuous margin
                    vineP.margins{2}.dist = 'gam';
                    vineP.margins{2}.theta = [0.1;10];
                    vineP.margins{2}.iscont = true; % Continuous margin
                    vineP.families = cell(2);
                    vineP.theta = cell(2);
                    vineP.families{1,2} = 'student';
                    vineP.theta{1,2} = [0;ro];          %%%%% ?ro
                    
                    X{4}=mixedvinernd(vineP,N);
                    
                case 5
                    
                    repp=mod(rep-1,9)+1;
                    ro=0.1*repp+1e-6*rand;%+0.1*rand(1)-0.05)*2;

                    vineP.margins = cell(2,1);
                    vineP.margins{1}.dist = 'poiss';
                    vineP.margins{1}.theta = 100*ro;
                    vineP.margins{1}.iscont = false; % Continuous margin
                    vineP.margins{2}.dist = 'poiss';
                    vineP.margins{2}.theta = 100*ro;
                    vineP.margins{2}.iscont = false; % Continuous margin
                    vineP.families = cell(2);
                    vineP.theta = cell(2);
                    vineP.families{1,2} = 'student';
                    vineP.theta{1,2} = [0;0.5];
                    
                    X{5}=mixedvinernd(vineP,N);
                    
                    
                case 6
                    
                    repp=mod(rep-1,9)+1;
                    ro=0.1*repp+1e-6*rand;%+0.1*rand(1)-0.05)*2;
                    
                    vineP.margins = cell(2,1);
                    vineP.margins{1}.dist = 'bino';
                    vineP.margins{1}.theta = [50;ro/2];
                    vineP.margins{1}.iscont = false; % Continuous margin
                    vineP.margins{2}.dist = 'bino';
                    vineP.margins{2}.theta = [50;1-ro/2];
                    vineP.margins{2}.iscont = false; % Continuous margin
                    vineP.families = cell(2);
                    vineP.theta = cell(2);
                    vineP.families{1,2} = 'student';
                    vineP.theta{1,2} = [0;0.5];

                    X{6}=mixedvinernd(vineP,N);
                    
                case 7
                    
                    repp=mod(rep-1,9)+1;
                    ro=0.1*repp+1e-6*rand;%+0.1*rand(1)-0.05)*2;
                    
                    vineP.margins = cell(2,1);
                    vineP.margins{1}.dist = 'poiss';
                    vineP.margins{1}.theta = 100*ro;
                    vineP.margins{1}.iscont = false; % Continuous margin
                    vineP.margins{2}.dist = 'poiss';
                    vineP.margins{2}.theta = 100*ro;
                    vineP.margins{2}.iscont = false; % Continuous margin
                    vineP.families = cell(2);
                    vineP.theta = cell(2);
                    vineP.families{1,2} = 'gaussian';
                    vineP.theta{1,2} = 0.5;
                    
                    X{7}=mixedvinernd(vineP,N);
                    
                case 8
                    
                    repp=mod(rep-1,9)+1;
                    ro=0.1*repp+1e-6*rand;%+0.1*rand(1)-0.05)*2;
                    
                    vineP.margins = cell(2,1);
                    vineP.margins{1}.dist = 'poiss';
                    vineP.margins{1}.theta = 100*ro;
                    vineP.margins{1}.iscont = false; % Continuous margin
                    vineP.margins{2}.dist = 'poiss';
                    vineP.margins{2}.theta = 100*ro;
                    vineP.margins{2}.iscont = false; % Continuous margin
                    vineP.families = cell(2);
                    vineP.theta = cell(2);
                    vineP.families{1,2} = 'gaussian';
                    vineP.theta{1,2} = 0.5;
                    
                    X{8}=mixedvinernd(vineP,N);
                    
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            for MODEL=MOD
                
                EX=0;
                if MODEL<7
%                     EX= (exist(['/n/data2/hms/neurobio/harvey/Houman/Copula_project/Info_Entropy/infos/nNew_info_MOD_',num2str(MODEL),'_N_',num2str(N),'_rep_',num2str(rep),'_knots_',num2str(knots),'.mat'],'file')==0);
                    EX= (exist(['/n/data2/hms/neurobio/harvey/Houman/Copula_project/Info_Entropy/infos/RevisionSep_nNew_info_MOD_',num2str(MODEL),'_N_',num2str(N),'_rep_',num2str(rep),'_knots_',num2str(knots),'.mat'],'file')==0);
                else
%                     EX= (exist(['/n/data2/hms/neurobio/harvey/Houman/Copula_project/Info_Entropy/infos/nNew_info_MOD_',num2str(7),'_N_',num2str(N),'_rep_',num2str(rep),'_knots_',num2str(knots),'.mat'],'file')==0 | exist(['/n/data2/hms/neurobio/harvey/Houman/Copula_project/Info_Entropy/infos//nNew_info_MOD_',num2str(8),'_N_',num2str(N),'_rep_',num2str(rep),'_knots_',num2str(knots),'.mat'],'file')==0); 
                      EX= (exist(['/n/data2/hms/neurobio/harvey/Houman/Copula_project/Info_Entropy/infos/RevisionSep_nNew_info_MOD_',num2str(MODEL),'_N_',num2str(N),'_rep_',num2str(rep),'_knots_',num2str(knots),'.mat'],'file')==0);
              end
                
                if  EX==1%exist(['/n/data2/hms/neurobio/harvey/Houman/Copula_project/Info_Entropy/infos//nNew_info_MOD_',num2str(MODEL),'_N_',num2str(N),'_rep_',num2str(rep),'_knots_',num2str(knots),'.mat'],'file')==0 
                    
                    [N,rp,MOD,knots]
                    
                    switch MODEL
                        case 1
                            I_theoretical=-0.5 * log2(1-ro^2);
                        case 2
                            I_theoretical=Omega_Student(ro,0,2);%0.5 * log2(sx^2 * sy^2/det(SIGMA));
                        case 3
                            I_theoretical=-0.5 * log2(1-ro^2);
                        case 4
                            I_theoretical=Omega_Student(ro,0,2);%0.5 * log2(sx^2 * sy^2/det(SIGMA));
                            
                        case 5
                            I_theoretical=-0.5 * log2(1-ro^2);
                        case 6
                            I_theoretical=Omega_Student(ro,0,2);%0.5 * log2(sx^2 * sy^2/det(SIGMA));
                        case 7
                            I_theoretical=-0.5 * log2(1-ro^2);
                        case 8
                            I_theoretical=Omega_Student(ro,0,2);%0.5 * log2(sx^2 * sy^2/det(SIGMA));
                    end
                    
                    if MODEL<5
                        range(1:2,1)=min(X{MODEL}(:,1:2))-1e-10;
                        range(1:2,2)=max(X{MODEL}(:,1:2))+1e-10;
                        [vine]=prep_copula(X{MODEL}(:,[1 2]),{'kernel','kernel'},{'kercop' 'kercop';'kercop' 'kercop'},[1 1],'c-vine','rand',range([1 2],:));
                    else
                        range(1:2,1)=min(X{MODEL}(:,1:2));
                        range(1:2,2)=max(X{MODEL}(:,1:2))+1;
                        [vine]=prep_copula(X{MODEL}(:,[1 2]),{'kernel','kernel'},{'discrete' 'discrete';'discrete' 'discrete'},[1 1],'c-vine','rand',range([1 2],:));
                    end
                    
                    for i = 1:2
                        iscont(i) = true;
                    end
                    vineest = mixedvinefit(X{MODEL},vine.type,iscont,[],[],'gaussian');
                    %             vineest = mixedvinefit(X{MODEL},vine.type,iscont,[],[]);
                    vine.condition=0;
                    
                    for i=1:2
                        for j=1:2
                            if MODEL<=50
                                vine.METH{i,j}=[1 1];
                            else
                                %                             vine.METH{i,j}=[0 0];
                            end
                        end
                    end
                    
                    
                    [~,~,copula_LL1,~,~] = Fit_vCopula(vine,[1 1],100,'LL1',1,0,'rand',[],knots);
                    
                    if knots==100 & N==1024 & MOD<5
                    [~,~,copula_LL2,~,~] = Fit_vCopula(vine,[1 1],100,'LL2',1,0,'rand',[],knots);
                    end
                    
                    %             knots=500;
                    [~,~,copula_LL1,~,~] = Fit_vCopula(vine,[1 1],100,'LL1',-1,copula_LL1,'rand',[],knots);
                    
                    if knots==100 & N==1024 & MOD<5
                    [~,~,copula_LL2,~,~] = Fit_vCopula(vine,[1 1],100,'LL2',-1,copula_LL2,'rand',[],knots);                    
                    end
                    
                    [info(1),stderr_1,in_1] = kernelvineinfo(vine,copula_LL1,20000,50)
                    
                    if knots==100 & N==1024 & MOD<5
                    [info(2),stderr_2,in_2] = kernelvineinfo(vine,copula_LL2,20000,50)
                    else
                    info(2)=-10;
                    stderr_2=NaN;
                    in_2=NaN;
                    copula_LL2=copula_LL1;
                    end
                    
                    [info(3),stderr_3,in_3] = kernelvineinfo(vineest,copula_LL1,20000,50,[],[],[],[],[],'par')

                    
                    %%%%%%%%%%%%%%%%% COMPUTE NAIVE information and also information from integral over copula and not sampling
                    clear dat
                    [pts,GRID_u]= mk_grid(knots,'');
                    for d = 1:2
                        par{d}.fit=-1;
                        par{d}.max=max(GRID_u(:));
                        par{d}.min=min(GRID_u(:));
                        [dat(:,d),Mar_T{d}]=kernelcdf(vine.margins{d}.ker,vine.margins{d}.ker,par{d});
                    end
                    data.u=dat;
                    data.S=norminv(data.u,0,1);    %%%%norminv
                    [COEFF,data.X,~,~,~,mu] = pca(data.S);
                    data.X=data.S * COEFF - repmat(mu,size(data.S,1),1) * COEFF;
                    data.u=double(data.u);
                    data.S=double(data.S);
                    data.X=double(data.X);
                    Grid=GRID_Bands(data.u,knots,[1 1]);
                    lfit_1 = loclik_fit(copula_LL1{1,2}.fit.bw,data,Grid);
                    lfit_2 = loclik_fit(copula_LL2{1,2}.fit.bw,data,Grid);
                    
                    
                    [x1,~,s1]=unique(sort(Grid.S(:,1)));
                    [x2,~,s2]=unique(sort(Grid.S(:,2)));
                    P1=normpdf(x1);
                    P2=normpdf(x2);
                    NORM=(P1*P2')';
                    x1=unique(Grid.u(:,1));
                    xd1=[diff(x1)];xd1=[xd1;xd1(1)];
                    x2=unique(Grid.u(:,2));
                    xd2=[diff(x2)];xd2=[xd2;xd2(1)];
                    
                    F_1=reshape(lfit_1.Kergrid,knots,knots)./NORM;
                    F_2=reshape(lfit_2.Kergrid,knots,knots)./NORM;
                    F_1_naive=reshape(lfit_1.KergridNaive{1},knots,knots)./NORM;
                    F_2_naive=reshape(lfit_2.KergridNaive{1},knots,knots)./NORM;
                    
                    
                    %%%%% normalization
                    
                    for u=1:4
                        switch u
                            case 1
                                tnorm=F_1;
                            case 2
                                tnorm=F_2;
                            case 3
                                tnorm=F_1_naive;
                            case 4
                                tnorm=F_2_naive;
                        end
                        
                        for n=1:1000
                            %                             tnorm(tnorm<1e-300)=1e-300;
                            I1=sum(xd2'.*tnorm,2);
                            I2=sum(xd1.*tnorm,1);
%                             I1=trapz(x2,tnorm,2);
%                             I2=trapz(x1,tnorm,1);
                            K=I1.*I2;
                            tnorm=tnorm./K;
                        end                        
                        II=sum(xd1.*sum(xd2'.*tnorm,2));
%                         II=trapz(x2,trapz(x1,tnorm,1));
                        tnorm=tnorm/II;
                        
                        switch u
                            case 1
                                F_1=tnorm;
                            case 2
                                F_2=tnorm;
                            case 3
                                F_1_naive=tnorm;
                            case 4
                                F_2_naive=tnorm;
                        end
                    end
                    
                    %%%%%%%%%% infos
                    
                    xd=xd1.*xd2';
                    info(8)=sum(sum(F_1(F_1~=0).*log2(F_1(F_1~=0)).*xd(F_1~=0)));
                    info(9)=sum(sum(F_2(F_2~=0).*log2(F_2(F_2~=0)).*xd(F_2~=0)));
                    info(10)=sum(sum(F_1_naive(F_1_naive~=0).*log2(F_1_naive(F_1_naive~=0)).*xd(F_1_naive~=0)));
                    info(11)=sum(sum(F_2_naive(F_2_naive~=0).*log2(F_2_naive(F_2_naive~=0)).*xd(F_2_naive~=0)));
                    
                    %%%%%%%%%%%%%%%%% kNN informations
                    
                    if MOD<5
                        load('/home/hs258/Codes_Folder/HOUMAN_GITLAB_IIT/TOOLS/ENTROPY_METHODS/NPEET_LNC/alpha_vec.mat')
                        clear ii jj
                        
                        mm=0;
                        for k=[3 5]%10
                            mm=mm+1;
                            cd /home/hs258/Codes_Folder/HOUMAN_GITLAB_IIT/TOOLS/ENTROPY_METHODS/NPEET_LNC
                            alpha=alpha_vec(alpha_vec(:,2)==k,3);
                            info_knn=knn_info(X{MODEL}(:,1),X{MODEL}(:,2),alpha,k);
                            cd /home/hs258/Codes_Folder/HOUMAN_GITLAB_IIT/TOOLS/ENTROPY_METHODS/NPEET_LNC
                            ii(mm)=info_knn.LNC;
                            jj(mm)=info_knn.kras;
                        end
                        
                        info(4)=ii(2);
                        info(5)=jj(2);
                        info(6)=ii(1);
                        info(7)=jj(1);
                    else
                        
                        r1=X{MODEL}(:,1);
                        r2=X{MODEL}(:,2);
                        [~,~,rj]=unique(X{MODEL},'rows');
                        
                        [mm1, icts1] = multiplicitiesFromSamples(r1);
                        [mm2, icts2] = multiplicitiesFromSamples(r2);
                        [mmj, ictsj] = multiplicitiesFromSamples(rj);
                        
                        [H1, Hvar1] = computeH_PYM(mm1, icts1);
                        [H2, Hvar2] = computeH_PYM(mm2, icts2);
                        [Hj, Hvarj] = computeH_PYM(mmj, ictsj);
                        
                        info(4)=(H1+H2-Hj)/log(2);
                        
                        
                        x1gv = 0:100;
                        x2gv = 0:100;
                        [x1,x2] = ndgrid(x1gv,x2gv);
                        p = mixedvinepdf(vineP,[x1(:),x2(:)]);
                        p=p(:);
                        p1=marginpdf(vineP.margins{1},x1gv);
                        p2=marginpdf(vineP.margins{2},x2gv);
                        H1=-sum(p1(p1~=0).*log2(p1(p1~=0)));
                        H2=-sum(p2(p2~=0).*log2(p2(p2~=0)));
                        H3=-sum(p(p~=0).*log2(p(p~=0)));
                        
                        I_theoretical=H1+H2-H3;
                        
%                         x1gv = 0:100;
%                         x2gv = 0:100;
%                         [x1,x2] = ndgrid(x1gv,x2gv);
%                         p0 = mixedvinepdf(vineP1,[x1(:),x2(:)]);
%                         p0=p0(:);
%                         p1=marginpdf(vineP1.margins{1},x1gv);
%                         p2=marginpdf(vineP1.margins{2},x2gv);
%                         p10 = mixedvinepdf(vineP2,[x1(:),x2(:)]);
%                         p10=p10(:);
%                         p11=marginpdf(vineP2.margins{1},x1gv);
%                         p12=marginpdf(vineP2.margins{2},x2gv);
%                         
%                         p1=0.5*p11+0.5*p1;
%                         p2=0.5*p12+0.5*p2;
%                         p=0.5*p0+0.5*p10;
%                         
%                         H1=-sum(p1(p1~=0).*log2(p1(p1~=0)));
%                         H2=-sum(p2(p2~=0).*log2(p2(p2~=0)));
%                         H3=-sum(p(p~=0).*log2(p(p~=0)));
%                         
%                         I_theoretical=H1+H2-H3;
                        
                    end
                    
                    err=(info-I_theoretical)
                    
                    clear lfit_1 lfit_2 copula_LL1 copula_LL2 Grid data GRID_u Mar_T F_1 F_2 F_1_naive F_2_naive K NORM pd_grid 
                    clear s1 s2 tnorm Xd vine X dat vineest dat vineest vineP pts par  %%% sep 2018
                    
                    if MODEL<5
%                         save(['/n/data2/hms/neurobio/harvey/Houman/Copula_project/Info_Entropy/infos/nNew_info_MOD_',num2str(MODEL),'_N_',num2str(N),'_rep_',num2str(rep),'_knots_',num2str(knots),'.mat'])
                        save(['/n/data2/hms/neurobio/harvey/Houman/Copula_project/Info_Entropy/infos/RevisionSep_nNew_info_MOD_',num2str(MODEL),'_N_',num2str(N),'_rep_',num2str(rep),'_knots_',num2str(knots),'.mat'])
%                         save(['/n/data2/hms/neurobio/harvey/Houman/Copula_project/Info_Entropy/infos/RevisionSeptember_info_MOD_',num2str(MODEL),'_N_',num2str(N),'_rep_',num2str(rp),'_knots_',num2str(knots),'_knotS_',num2str(knots_S),'.mat'])
                    else
%                         save(['/n/data2/hms/neurobio/harvey/Houman/Copula_project/Info_Entropy/infos/nNew_info_MOD_',num2str(MODEL),'_N_',num2str(N),'_rep_',num2str(rep),'_knots_',num2str(knots),'.mat'])
                        save(['/n/data2/hms/neurobio/harvey/Houman/Copula_project/Info_Entropy/infos/RevisionSep_nNew_info_MOD_',num2str(MODEL),'_N_',num2str(N),'_rep_',num2str(rep),'_knots_',num2str(knots),'.mat'])
                    end
                    
                end
                
            end
            
            
        end
    end
    
    delete(gcp('nocreate'))

% catch
%     
%     s = lasterror;
%     ss=s.message;
%     aa=round(rand*1e4);
%     save(['/n/data2/hms/neurobio/harvey/Houman/Copula_project/Info_Entropy/infos/err_',num2str(aa),'.mat'],'ss','MOD','N','rep')
%     
% end


return

mmm=0;
% for knots_S=[10 30 50 100 150 200 500]
for NN=5:13%13%11:12%5:13
    for rep=1:3:1000
        for MOD=[1 2 3 4 5 7]%1:6
            for knots= 100%[150 200]%[10 30 50 100 150 200]
                N=2^NN;
                EX=0;
                for repn=(rep-1)*9+1:rep*9
                    %                         if exist(['/n/data2/hms/neurobio/harvey/Houman/Copula_project/Info_Entropy/infos/New_info_MOD_',num2str(MOD),'_N_',num2str(N),'_rep_',num2str(repn),'_knots_',num2str(knots),'.mat'],'file')==0
                    if exist(['/n/data2/hms/neurobio/harvey/Houman/Copula_project/Info_Entropy/infos/RevisionSep_nNew_info_MOD_',num2str(MOD),'_N_',num2str(N),'_rep_',num2str(repn),'_knots_',num2str(knots),'.mat'],'file')==0
                        %                         if exist(['/n/data2/hms/neurobio/harvey/Houman/Copula_project/Info_Entropy/infos/RevisionSeptember_info_MOD_',num2str(MOD),'_N_',num2str(N),'_rep_',num2str(repn),'_knots_',num2str(knots),'_knotS_',num2str(knotS),'.mat'],'file')==0
                        EX=EX+1;
                    end
                end
                %                     repn=(rep-1)*9+1;
                out='/n/data2/hms/neurobio/harvey/Houman/TEMP/SLURM';
                if EX~=0%exist(['/n/data2/hms/neurobio/harvey/Houman/Copula_project/Info_Entropy/infos/info_MOD_',num2str(MOD),'_N_',num2str(N),'_rep_',num2str(repn),'_knots_',num2str(knots),'.mat'],'file')==0
                    if NN<9
                        system(['sbatch -p short -n 1 -t 0-2:00 --mem-per-cpu=500M -o /n/data2/hms/neurobio/harvey/Houman/TEMP/SLURM/%N_%j.out --wrap="matlab -nodisplay -r  \"run /home/hs258/Codes_Folder/Houman_Git/Ker_Copula/demo_info_entropy(',num2str(N),',',num2str(rep),',',num2str(MOD),',',num2str(knots),')\""'])
                    else
%                         mmm=mmm+1;
%                         c=parcluster;
%                         c.AdditionalProperties.WallTime = '2:00:00';
%                         c.AdditionalProperties.QueueName = 'mpi';
%                         c.AdditionalProperties.AdditionalSubmitArgs = '--mem-per-cpu=1G';
%                         JJ{mmm}=c.batch(@demo_info_entropy,4,{(N),(rep),(MOD),(knots)},'Pool',2);
                        system(['sbatch -p short -n 1 -t 0-2:00 --mem-per-cpu=2G -o /n/data2/hms/neurobio/harvey/Houman/TEMP/SLURM/%N_%j.out --wrap="matlab -nodisplay -r  \"run /home/hs258/Codes_Folder/Houman_Git/Ker_Copula/demo_info_entropy(',num2str(N),',',num2str(rep),',',num2str(MOD),',',num2str(knots),')\""'])

%                           demo_info_entropy(N,rep,MOD,knots)
                    end
                end
            end
        end
    end
end
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%  Schematic figures for the copula and sklar
%%%%%%%%%%%%%%%%%%%%%%%%%  theorem

                    ro=0.4;                    
                    vineP1.margins = cell(2,1);
                    vineP1.margins{1}.dist = 'gam';
                    vineP1.margins{1}.theta = [0.7;10];
                    vineP1.margins{1}.iscont = true; % Continuous margin
                    vineP1.margins{2}.dist = 'norm';
                    vineP1.margins{2}.theta = [0;1];
                    vineP1.margins{2}.iscont = true; % Continuous margin
                    vineP1.families = cell(2);
                    vineP1.theta = cell(2);
                    vineP1.families{1,2} = 'student';
                    vineP1.theta{1,2} = [0;ro];          %%%%% ?ro

                    C=mixedvinernd(vineP1,2000); 
                    x1gv = linspace(0,50,200);
                    x2gv = linspace(-5,5,200);

                    knots=30;
                    range(1:2,1)=min(C(:,1:2))-1e-10;
                    range(1:2,2)=max(C(:,1:2))+1e-10;
                    [vine]=prep_copula(C(:,[1 2]),{'kernel','kernel'},{'kercop' 'kercop';'kercop' 'kercop'},[1 1],'c-vine','rand',range);

                    figure;plot(vine.margins{1}.ker,vine.margins{2}.ker,'.');axis off

                    for i = 1:2
                        iscont(i) = true;
                    end
                    vine.condition=0;

                    for i=1:2
                        for j=1:2
                                vine.METH{i,j}=[1 1];
                        end
                    end
                    
                    mm=0;clear points
                    for i=1:numel(x1gv)
                        for j=1:numel(x2gv)
                            mm=mm+1;
                            points(mm,:)=[x1gv(i) x2gv(j)];
                        end
                    end
                    
                    [~,~,copula_LL1,~,~] = Fit_vCopula(vine,[1 1],100,'LL1',1,0,'rand',[],knots);
                    [f1,f2,copula_LL1,~,pc] = Fit_vCopula(vine,points,100,'LL1',-1,copula_LL1,'rand',[],knots);
        
                    clear dat
                    [pts,GRID_u]= mk_grid(knots,'');
                    for d = 1:2
                        par{d}.fit=-1;
                        par{d}.max=max(GRID_u(:));
                        par{d}.min=min(GRID_u(:));
                        [dat(:,d),Mar_T{d}]=kernelcdf(vine.margins{d}.ker,vine.margins{d}.ker,par{d});
                    end
                    data.u=dat;
                    data.S=norminv(data.u,0,1);    %%%%norminv
                    [COEFF,data.X,~,~,~,mu] = pca(data.S);
                    data.X=data.S * COEFF - repmat(mu,size(data.S,1),1) * COEFF;
                    data.u=double(data.u);
                    data.S=double(data.S);
                    data.X=double(data.X);
                    Grid=GRID_Bands(data.u,knots,[1 1]);
                   
                    
                    X1=unique(Grid.u(:,1));
                    figure
                    surf(X1(5:end-5),X1(5:end-5),copula_LL1{1,2}.f_grid(5:end-5,5:end-5))
                    

                    c1=linspace(0,1-1e-100,200);
                    c2=linspace(0,1-1e-100,200);
                    xgv1=margininv(vineP1.margins{1},c1);
                    xgv2=margininv(vineP1.margins{2},c2);
                    pd1=marginpdf(vineP1.margins{1},xgv1);
                    pd2=marginpdf(vineP1.margins{2},xgv2);
                    cd1=margincdf(vineP1.margins{1},xgv1);
                    cd2=margincdf(vineP1.margins{2},xgv2);
                    
                    [x1,x2] = ndgrid(xgv1,xgv2);
                    
                    p = mixedvinepdf(vineP1,[x1(:),x2(:)]);
                    p1=marginpdf(vineP1.margins{1},x1(:));
                    p2=marginpdf(vineP1.margins{2},x2(:));
                    u1=margincdf(vineP1.margins{1},x1(:));
                    u2=margincdf(vineP1.margins{2},x2(:));
                    
                    pp=reshape(p./(p1.*p2),numel(x1gv),numel(x2gv));
                    pd=reshape(p,numel(x1gv),numel(x2gv));
                    
                    pdata = mixedvinepdf(vineP1,C);
                    COL=log2(pdata)
                    

                    figure(11)
                    subplot(3,3,5)
                    pcolor(c1(10:end-10),c2(10:end-10),pp(10:end-10,10:end-10)/max(pp(:)))
                    shading interp
%                     caxis([0 0.01])
                    caxis([0 0.04])
                    title('Student-t copula')
                    xlabel('u_1','fontsize',10)
                    ylabel('u_2','fontsize',10)
                    colormap('jet')
                    axis square
                    
                    subplot(3,3,1)
                    %                     [hAx,hLine1,hLine2] = plotyy(xgv1(xgv1<40 & xgv1>1e-1),pd1(xgv1<40 & xgv1>1e-1),xgv1(xgv1<40 & xgv1>1e-1),cd1(xgv1<40 & xgv1>1e-1))
                    [hAx,hLine1,hLine2] = plotyy(xgv1(xgv1<40),pd1(xgv1<40),xgv1(xgv1<40),cd1(xgv1<40))
                    ylabel(hAx(1),'PDF')
                    ylabel(hAx(2),'CDF')
                    %                     axis tight
                    %                     axis square
                    ylim([0 0.4])
                    hLine1.LineStyle = '--';
                    hLine2.LineStyle = '-';
                    hLine1.LineWidth = 2;
                    hLine2.LineWidth = 2;
                    title('Gamma marginal distribution')
                    xlabel('x_1','fontsize',10)
                    
                    subplot(3,3,7)
                    [hAx,hLine1,hLine2] = plotyy(xgv2,pd2,xgv2,cd2)
                    ylabel(hAx(1),'PDF')
                    ylabel(hAx(2),'CDF')
                    %                     axis tight
                    %                     axis square
                    hLine1.LineStyle = '--';
                    hLine2.LineStyle = '-';
                    hLine1.LineWidth = 2;
                    hLine2.LineWidth = 2;
                    title('Gaussian marginal distribution')
                    xlabel('x_2','fontsize',10)
                    
                    subplot(3,3,6)
                    scatter(C(:,1),C(:,2),20,COL,'.')
                    colormap('jet')
                    title('2D scatter')
                    axis square
                    xlim([-0.1 40])
                    xlabel('x_1','fontsize',10)
                    ylabel('x_2','fontsize',10)
                    
                    
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    ro=0.3;
                    vineP1.margins = cell(2,1);
                    vineP1.margins{1}.dist = 'norm';
                    vineP1.margins{1}.theta = [0;1];
                    vineP1.margins{1}.iscont = true; % Continuous margin
                    vineP1.margins{2}.dist = 'norm';
                    vineP1.margins{2}.theta = [0;1];
                    vineP1.margins{2}.iscont = true; % Continuous margin
                    vineP1.families = cell(2);
                    vineP1.theta = cell(2);
                    vineP1.families{1,2} = 'clayton';
                    vineP1.theta{1,2} = [2];          %%%%% ?ro

                    C=mixedvinernd(vineP1,5000);
               
                    knots=50;
                    range(1:2,1)=min(C(:,1:2))-1e-10;
                    range(1:2,2)=max(C(:,1:2))+1e-10;
                    [vine]=prep_copula(C(:,[1 2]),{'kernel','kernel'},{'kercop' 'kercop';'kercop' 'kercop'},[1 1],'c-vine','rand',range);
          
                    for i = 1:2
                        iscont(i) = true;
                    end
                    vine.condition=0;

                    for i=1:2
                        for j=1:2
                                vine.METH{i,j}=[1 1];
                        end
                    end
                    [~,~,copula_LL1,~,~] = Fit_vCopula(vine,[1 1],100,'LL1',1,0,'rand',[],knots);
                    [f1,f2,copula_LL1,~,pc] = Fit_vCopula(vine,points,100,'LL1',-1,copula_LL1,'rand',[],knots);
        

                    clear dat
                    [pts,GRID_u]= mk_grid(knots,'');
                    for d = 1:2
                        par{d}.fit=-1;
                        par{d}.max=max(GRID_u(:));
                        par{d}.min=min(GRID_u(:));
                        [dat(:,d),Mar_T{d}]=kernelcdf(vine.margins{d}.ker,vine.margins{d}.ker,par{d});
                    end
                    data.u=dat;
                    data.S=norminv(data.u,0,1);    %%%%norminv
                    [COEFF,data.X,~,~,~,mu] = pca(data.S);
                    data.X=data.S * COEFF - repmat(mu,size(data.S,1),1) * COEFF;
                    data.u=double(data.u);
                    data.S=double(data.S);
                    data.X=double(data.X);
                    Grid=GRID_Bands(data.u,knots,[1 1]);
                   
                    cdata1 = margincdf(vineP1.margins{1},C(:,1));
                    cdata2 = margincdf(vineP1.margins{2},C(:,2));
                                                                        
                    lfit_1 = loclik_fit(copula_LL1{1,2}.fit.bw,data,Grid);
                    pc1=normcdf(-2.1);
                    pc2=normcdf(-1.8);
                    a=ellipse(0.4,0.7,-pi/4,-2.1,-1.8,'r');
                    ec1=normcdf(a.XData);
                    ec2=normcdf(a.YData);
                    
                    pc3=normcdf(1.8);
                    pc4=normcdf(0);
                    a=ellipse(0.4,0.7,-pi/4,1.8,0,'r');
                    ec3=normcdf(a.XData);
                    ec4=normcdf(a.YData);
                    
                    
                    figure(12)
                    subplot(1,2,1)
                    %                     xticks(unique(Grid.u(:,1)))
                    %                     yticks(unique(Grid.u(:,1)))
                    hold on
                    g_y=unique(Grid.u(:,1));
                    g_x=unique(Grid.u(:,1));
                    for i=1:length(g_x)
                        plot([g_x(i) g_x(i)],[g_y(1) g_y(end)],'-','color',[0.4 0.4 0.4],'linewidth',0.3) %y grid lines
                        hold on
                    end
                    for i=1:length(g_y)
                        plot([g_x(1) g_x(end)],[g_y(i) g_y(i)],'-','color',[0.4 0.4 0.4],'linewidth',0.3) %y grid lines
                        hold on
                    end
%                     set(gca,'xticklabel',{[]})
%                     set(gca,'yticklabel',{[]})
%                     grid off
                    xlim([0 1])
                    ylim([0 1])
                    axis square
                    box off
                    title('(u_1,u_2)')
                    scatter(data.u(1:1000,1),data.u(1:1000,2),40,'k','.')
                    hold on
                    plot(pc1,pc2,'ok','markersize',3)
                    hold on
                    plot(ec1,ec2,'r')
                    hold on
%                     fill(ec1,ec2,'r')
                    alpha(0.5)
                    hold on
                    plot(pc3,pc4,'ok','markersize',3)
                    hold on
                    plot(ec3,ec4,'b')
                    hold on
%                     fill(ec3,ec4,'b')
                    alpha(0.4)
                    
                    subplot(1,2,2)
                    %                     xticks(unique(Grid.S(:,1)))
                    %                     yticks(unique(Grid.S(:,1)))
                    hold on
                    g_y=unique(Grid.S(:,1));
                    g_x=unique(Grid.S(:,1));
                    for i=1:length(g_x)
                        plot([g_x(i) g_x(i)],[g_y(1) g_y(end)],'-','color',[0.4 0.4 0.4],'linewidth',0.3) %y grid lines
                        hold on
                    end
                    for i=1:length(g_y)
                        plot([g_x(1) g_x(end)],[g_y(i) g_y(i)],'-','color',[0.4 0.4 0.4],'linewidth',0.3) %y grid lines
                        hold on
                    end
%                     set(gca,'xticklabel',{[]})
%                     set(gca,'yticklabel',{[]})
%                     grid off
                    xlim([-3.4 3.4])
                    ylim([-3.4 3.4])
                    axis square
                    box off
                    scatter(data.S(1:1000,1),data.S(1:1000,2),40,'k','.')
                    title('(\Phi^{-1}(u_1),\Phi^{-1}(u_2))')
                    hold on
                    a=ellipse(0.4,0.7,-pi/4,-2.1,-1.8,'r')
                    hold on
                    plot(-2.1,-1.8,'ok','markersize',3)
                    hold on
%                     fill(a.XData,a.YData,'r')
                    alpha(0.5)
                    hold on
                    a=ellipse(0.4,0.7,-pi/4,1.8,0,'b')
                    hold on
                    plot(1.8,0,'ok','markersize',3)
                    hold on
%                     fill(a.XData,a.YData,'b')
                    alpha(0.4)
                    ylim([min(Grid.S(:,1)) max(Grid.S(:,1))])
                    xlim([min(Grid.S(:,1)) max(Grid.S(:,1))])

                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    c1=unique(Grid.u(:,1));
                    c2=unique(Grid.u(:,2));
                    pp=lfit_1.Kergrid;
                    p=reshape(pp,knots,knots);
                    figure(13)                
                    subplot(1,6,1)
                    pcolor(c1,c2,p/max(p(:)))
                    shading interp
                    caxis([0 1])
                    title('Student-t copula')
                    xlabel('u_1','fontsize',10)
                    ylabel('u_2','fontsize',10)
                    colormap('jet')
                    axis square
                    for i=1:5
                        subplot(1,6,i+1)
                        pp=lfit_1.KergridNaive{i};
                        p=reshape(pp,knots,knots);
                        pcolor(c1,c2,p/max(p(:)))
                        shading interp
                        caxis([0 1])
                        title('Student-t copula')
                        xlabel('u_1','fontsize',10)
                        ylabel('u_2','fontsize',10)
                        colormap('jet')
                        axis square
                    end

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%% discrete to continous
                    
                    ro=0.1;
                    vineP2.margins = cell(2,1);
                    vineP2.margins{1}.dist = 'poiss';
                    vineP2.margins{1}.theta = 100*ro;
                    vineP2.margins{1}.iscont = false; % Continuous margin
                    vineP2.margins{2}.dist = 'poiss';
                    vineP2.margins{2}.theta = 100*ro;
                    vineP2.margins{2}.iscont = false; % Continuous margin
                    vineP2.families = cell(2);
                    vineP2.theta = cell(2);
                    vineP2.families{1,2} = 'student';
                    vineP2.theta{1,2} = [0;0.2];
                    
                    D=mixedvinernd(vineP2,2000);
                    
                    x1gv = 0:100;
                    x2gv = 0:100;
                    [x1,x2] = ndgrid(x1gv,x2gv);
                    p = mixedvinepdf(vineP2,[x1(:),x2(:)]);


                    pdata = mixedvinepdf(vineP2,D);
                    COL=log2(pdata);
                    
                    range(1:2,1)=min(D(:,1:2));
                    range(1:2,2)=max(D(:,1:2))+1;
                    [vine]=prep_copula(D(:,[1 2]),{'kernel','kernel'},{'discrete' 'discrete';'discrete' 'discrete'},[1 1],'c-vine','rand',range([1 2],:));
                    
                        
                    figure
                    subplot(1,2,1)
%                     plot(D(:,1),D(:,2),'Ok','markersize',3)
                    scatter(D(:,1),D(:,2),10,COL,'o')
                    axis square
                    subplot(1,2,2)
                    scatter(vine.margins{1}.ker,vine.margins{2}.ker,20,COL,'.')
                    colormap('jet')
                    axis square
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    
                    
                    
                    
                    