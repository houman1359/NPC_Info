
clear all
close all


rmpath(genpath('/home/hs258/Codes_Folder/HOUMAN_GITLAB_IIT/Copula_Arno'))
addpath(genpath('/home/hs258/Codes_Folder/HOUMAN_GITLAB_IIT/Selmaan_Switch/glmnet'))
addpath(genpath('/home/hs258/Codes_Folder/Houman_Git/Ker_Copula'))



%%% 9--> n-B dep.
%%% 11 --> cond marginals
%%% 12 --> B1-B2 dep.

simulation_num=18%17%12%11;
Simulation_for_Copula
% X=X+1e-10*rand(size(X));

if simulation_num==5
SH=2;PP=2;
X=cat(2,POP_st1{SH,PP},POP_st2{SH,PP});
X(size(X,1)+1,:)=[ones(size(POP_st1{SH,PP}(1,:))) 2*ones(size(POP_st2{SH,PP}(1,:)))]; 
X=X';

range(:,1)=min(X(:,1:2));
range(:,2)=max(X(:,1:2));
end

% xx=X{1}(:,1);
% yy=X{1}(:,2);
% X{1}(:,1)=yy;
% X{1}(:,2)=xx;
% xx=X{2}(:,1);
% yy=X{2}(:,2);
% X{2}(:,1)=yy;
% X{2}(:,2)=xx;

if simulation_num==6 | simulation_num==7 | simulation_num==8 | simulation_num==9 | simulation_num==10 | simulation_num==11 | simulation_num==12 | simulation_num==13 | simulation_num==14 | simulation_num==15 | simulation_num==16 | simulation_num==17 | simulation_num==18 | simulation_num==19

X1=cat(1,X{1},X{2});
X1(:,size(X{1},2)+1)=[ones(size(X{1}(:,1)));2*ones(size(X{2}(:,1)))];
X=X1;

if simulation_num==13
%     X=X(:,[1 200:201]);
    X=X(:,[1 20:21]);
end

if simulation_num~=14 &  simulation_num~=17 &  simulation_num~=18 &  simulation_num~=19
range(:,1)=min(X(:,1:2));
range(:,2)=max(X(:,1:2));
else
range(:,1)=min(X(:,1:3));
range(:,2)=max(X(:,1:3));
end
end

% range(:,1)=min(X(:,1:2));
% range(:,2)=max(X(:,1:2));


X0=X;
X00=X;


for stim=1:2
X=X0(X0(:,end)==stim,1:end-1);
% X=X0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% design the vine
clear vine
for d=1:size(X,2)
margins{d}='kernel';
if d==1
iscons(d)=0;
else
iscons(d)=0;
end
end

for i=1:size(X,2)
    for j=1:size(X,2)
        families{i,j}='kercop';
    end
end

mm=1;

[vine]=prep_copula(X,margins,families,iscons,'c-vine','rand',range)
% y=X(:,1);
% Xpred=X(:,2:end);
y=X(:,mm);
Xpred=X(:,setdiff(1:size(X,2),mm));


if simulation_num==17 | simulation_num==18 | simulation_num==19
[vine_1{stim}]=prep_copula(X(:,[1 2]),{'kernel','kernel'},{'kercop' 'kercop';'kercop' 'kercop'},iscons([1 2]),'c-vine','rand',range([1 2],:));
[vine_2{stim}]=prep_copula(X(:,[1 3]),{'kernel','kernel'},{'kercop' 'kercop';'kercop' 'kercop'},iscons([1 3]),'c-vine','rand',range([1 3],:));
[vine_3{stim}]=prep_copula(X(:,[2 3]),{'kernel','kernel'},{'kercop' 'kercop';'kercop' 'kercop'},iscons([2 3]),'c-vine','rand',range([2 3],:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   GLM

GLM=0;
if GLM==1
nFolds = 10;
trialIDs = 1:size(Xpred,1); 
trialPartition = cvpartition(length(trialIDs), 'KFold', nFolds);
foldIDs = nan(size(Xpred,1),1);
for nFold = 1:nFolds
    foldTrials = trialIDs(find(test(trialPartition,nFold)));  
    foldFrames = find(ismember(1:size(Xpred,1),foldTrials));
    foldIDs(foldFrames) = nFold;
end

opt = glmnetSet;
opt.alpha = 1/2;
opt.thresh = 1e-6;

GG = cvglmnet(Xpred,y,'gaussian',opt,'deviance',[],foldIDs);
GLM_prediction{stim} = cvglmnetPredict(GG,Xpred,[],'response');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Copula

%%%%%%%% points: the set of points which we want to compute the joint

y_vector =linspace(min(y)-eps,max(y)+eps,50);%unique(y);
points=zeros(numel(y_vector)*size(Xpred,2),size(Xpred,2)+1);

set_train=1:size(Xpred,1);%randsample(1:size(Xpred,1),round(0.9*size(Xpred,1)));%1:size(Xpred,1);%
set_test=1:size(Xpred,1);%setdiff(1:size(Xpred,1),set_train);%1:size(Xpred,1);%

clear points_train points_test
ne=0;
for j=1:numel(set_train)
    for i=1:numel(y_vector)
        ne=ne+1;
%         points_train(ne,:)=[y_vector(i) Xpred(set_train(j),:)];
        Z(mm)=y_vector(i);
        Z(setdiff(1:size(Xpred,2)+1,mm))=X(set_train(j),setdiff(1:size(Xpred,2)+1,mm));
        points_train(ne,:)=Z;
    end
end


ne=0;
for j=1:numel(set_test)
    for i=1:numel(y_vector)
        ne=ne+1;
        points_test(ne,:)=[y_vector(i) Xpred(set_test(j),:)];
    end
end

%%%%%%%% density

%%%%%%%%%%% we fit the bw to all the data and then we can estimate the f
%%%%%%%%%%% with 10-fold CV. We don't need to fit bw to each fold. because
%%%%%%%%%%% the bw is estimated in with CV.
vine_stim{stim}=vine;
TL=20;   %%% this is the depth of the vine.
tic

for i = 1:d
    iscont(i) = true;%vine.margins{i}.iscont;
end
vineest{stim} = mixedvinefit(X,vine.type,iscont);

vine_stim{stim}.condition=0;
knots=50;
% knots=5;

for i=1:d
    for j=1:d
vine_stim{stim}.METH{i,j}='LL1';
    end
end
for i=1:2
    for j=1:2
vine_1{stim}.METH{i,j}='LL1';
vine_2{stim}.METH{i,j}='LL1';
vine_3{stim}.METH{i,j}='LL1';
    end
end

tic
[f_pointR,f_data1,copula{stim},~] = Fit_vCopula(vine_stim{stim},points_train(1,:),TL,'LL1',1,0,'rand',vineest{stim},knots);   %%%% the method 'LL1', 'LL2', 'fix' 'nn' for analytical and 'TLL1', 'TLL1nn' 'TLL2' 'TLL2nn' for the R package
[f_pointR,f_data1,copula{stim},~] = Fit_vCopula(vine_stim{stim},points_train(1,:),TL,'LL1',-1,copula{stim},'rand',vineest{stim},knots);   %%%% the method 'LL1', 'LL2', 'fix' 'nn' for analytical and 'TLL1', 'TLL1nn' 'TLL2' 'TLL2nn' for the R package

if simulation_num==17 | simulation_num==18 | simulation_num==19
vine_1{stim}.condition=0;
vine_2{stim}.condition=0;
vine_3{stim}.condition=0;

[f_pointR_1,f_data1_1,copula_1{stim},~] = Fit_vCopula(vine_1{stim},points_train(1,[1 2]),TL,'LL1',1,0,'rand',vineest{stim},knots);   %%%% the method 'LL1', 'LL2', 'fix' 'nn' for analytical and 'TLL1', 'TLL1nn' 'TLL2' 'TLL2nn' for the R package
[f_pointR_2,f_data1_2,copula_2{stim},~] = Fit_vCopula(vine_2{stim},points_train(1,[1 2]),TL,'LL1',1,0,'rand',vineest{stim},knots);   %%%% the method 'LL1', 'LL2', 'fix' 'nn' for analytical and 'TLL1', 'TLL1nn' 'TLL2' 'TLL2nn' for the R package
[f_pointR_3,f_data1_3,copula_3{stim},~] = Fit_vCopula(vine_3{stim},points_train(1,[1 2]),TL,'LL1',1,0,'rand',vineest{stim},knots);   %%%% the method 'LL1', 'LL2', 'fix' 'nn' for analytical and 'TLL1', 'TLL1nn' 'TLL2' 'TLL2nn' for the R package

[f_pointR_3,f_data1_3,copula_3{stim},~] = Fit_vCopula(vine_3{stim},points_train(1,[1 2]),TL,'LL1',-1,copula_3{stim},'rand',vineest{stim},knots);   %%%% the method 'LL1', 'LL2', 'fix' 'nn' for analytical and 'TLL1', 'TLL1nn' 'TLL2' 'TLL2nn' for the R package

end

toc

[vine_train{stim}]=prep_copula(X(set_train,:),margins,families,iscons,'c-vine','rand',range);
vine_train{stim}.condition=0;
for i=1:d
    for j=1:d
vine_train{stim}.METH{i,j}='LL1';
    end
end

[f_points_train{stim},f_data2,copula_train{stim},pdf{stim}] = Fit_vCopula(vine_train{stim},points_train,TL,[],-1,copula{stim},'rand',[],knots);   
% [vine_test{stim}]=prep_copula(X(set_test,:),margins,families,iscons,'c-vine','rand',range);
% [f_points_test{stim},f_data2,~,pdf{stim}] = Fit_vCopula(vine_test{stim},points_test,TL,[],0,copula_train{stim},'rand',[]);

toc

% plot_copula(copula{1})

f_pointsS=reshape(f_points_train{stim},numel(y_vector),[]);
%%%%%%%% density
[y_mle,y_em,LL]=predict_response(f_pointsS,y_vector,y(set_test));

Y_EM{stim}=y_em;
Y_MLE{stim}=y_mle;
LLL{stim}=LL;

yt=y(set_test);
% title([num2str(corr(yt(:),y_mle(:),'type','spearman')),' , ',num2str(corr(yt(~isnan(y_em)),y_em(~isnan(y_em))','type','spearman')),' , ',num2str(corr(yt(:),GLM_prediction{stim}(set_test),'type','spearman'))])
[num2str(corr(yt(:),y_mle(:))),' , ',num2str(corr(yt(~isnan(y_em(:))),y_em(~isnan(y_em(:)))'))]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i=1:size(X0,2)
xgv{i} = linspace(min(X0(:,i)),max(X0(:,i)),50);
dx(i)=mean(diff(xgv{i}));
end

[xv{1},xv{2}] = ndgrid(xgv{1},xgv{2});
Xf{stim}=[xv{1}(:) xv{2}(:)];


end

if simulation_num==18 | simulation_num==19
    
         t_vector=linspace(min(vine_train{1}.margins{2}.ker)-0.01,max(vine_train{1}.margins{2}.ker)+0.01,25);
        
         for m=1:2
             vine_train{m}.B=vine_3{m};
             vine_train{m}.copulaB=copula_3{m};
         end

         [inf,std_inf,in,En] = kernelvineinfo(vine_train,copula,5000,1,[1 2],[0.5 0.5],1,2,t_vector);%kernelvineinfo(vine_train,copula,[0.5 0.5],[],3e-3,2000,'copcop',[],5,'condition',1,1,t_vector)         
         [inf_B,std_inf_B,in_B,En_B] = kernelvineinfo(vine_3,copula_3,20000,1,[1 2],[0.5 0.5],1,1,t_vector);%kernelvineinfo(vine_train,copula,[0.5 0.5],[],3e-3,2000,'copcop',[],5,'condition',1,1,t_vector)

         X=in.info(4,:)-in_B.info(3,:);X(X<0)=0;figure;plot(smooth(X))

         
         
         [~,h1]=histc(vine_train{1}.margins{2}.ker,t_vector);
         [~,h2]=histc(vine_train{2}.margins{2}.ker,t_vector);
         addpath(genpath('/home/hs258/Codes_Folder/TOOLS/Information Breakdown ToolBox - v.1.0.5'))
         nb=5;
         opts.method = 'dr';
         opts.bias = 'pt';
         DATA(:,1)=cat(1,vine_train{1}.margins{1}.ker,vine_train{2}.margins{1}.ker);
         DATA(:,2)=cat(1,vine_train{1}.margins{3}.ker,vine_train{2}.margins{3}.ker);
         DATA(:,3)=cat(1,vine_train{1}.margins{2}.ker,vine_train{2}.margins{2}.ker);
         stim=[ones(size(vine_train{1}.margins{1}.ker));2*ones(size(vine_train{2}.margins{1}.ker))];
         [~,h]=histc(DATA(:,3),t_vector);
         for t=1:max(h2)
             
             X1=DATA(h==t,1);
             st=stim(h==t);
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             nt=size(X1,1);
             S = binr(X1', nt, nb, 'eqpop');
             NS=buildr(st',S);
             opts.nt=[];
             for i=unique(st)'
                 opts.nt(i) = sum(st==i);
             end
             I_n(t) = information(NS, opts, 'I');
   
             X1=DATA(h==t,2);
             st=stim(h==t);
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             nt=size(X1,1);
             S = binr(X1', nt, nb, 'eqpop');
             NS=buildr(st',S);
             opts.nt=[];
             for i=unique(st)'
                 opts.nt(i) = sum(st==i);
             end
             I_B(t) = information(NS, opts, 'I');
    
         end
         
figure;plot(in.info(4,:),'k');hold on;plot(I_B,'b');hold on;plot(in_B.info(1,:),'r')
         
%     t_vector=linspace(min(vine_train{1}.margins{1}.ker),max(vine_train{1}.margins{1}.ker),40);
%     [inf,std_inf,in,En] = kernelvineinfo(vine_train,copula,[0.5 0.5],[],3e-3,1000,'copcop',[],5,'condition',1,2,t_vector)
    
end

figure;plot(in.info(4,:),'k');hold on;plot(I_B,'b');hold on;plot(in_B.info(1,:),'r');hold on;plot(I_n,'g');hold on;plot(in.info(3,:),'y')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   GLM    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   GLM with condition included for stimu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   case 6,7

nFolds = 10;
for con=1:2
yp=vine_stim{con}.margins{1}.ker;
trialIDs = 1:size(y,1); 
trialPartition = cvpartition(length(trialIDs), 'KFold', nFolds);
foldIDs_train{con} = nan(size(yp,1),1);
for nFold = 1:nFolds
    foldTrials = trialIDs(find(test(trialPartition,nFold)));  
    foldFrames = find(ismember(1:size(yp,1),foldTrials));
    foldIDs_train{con}(foldFrames) = nFold;
end
end

GLM_pred=nan(size(vine_stim{1}.margins{1}.ker,1)+size(vine_stim{2}.margins{1}.ker,1),1);
GLM_pred_int=nan(size(vine_stim{1}.margins{1}.ker,1)+size(vine_stim{2}.margins{1}.ker,1),1);

clear GG beta a0

for foldn=1:10

     
y_train=cat(1,vine_stim{1}.margins{1}.ker(foldIDs_train{1}~=foldn),vine_stim{2}.margins{1}.ker(foldIDs_train{2}~=foldn));
y_test=cat(1,vine_stim{1}.margins{1}.ker(foldIDs_train{1}==foldn),vine_stim{2}.margins{1}.ker(foldIDs_train{2}==foldn));

for j=1:size(X0,2)-2
XpredC_train=cat(1,vine_stim{1}.margins{1+j}.ker(foldIDs_train{1}~=foldn),vine_stim{2}.margins{1+j}.ker(foldIDs_train{2}~=foldn));
XpredC_test=cat(1,vine_stim{1}.margins{1+j}.ker(foldIDs_train{1}==foldn),vine_stim{2}.margins{1+j}.ker(foldIDs_train{2}==foldn));
if j==1
Xpred_train=XpredC_train;
Xpred_test=XpredC_test;
else
Xpred_train=cat(2,Xpred_train,XpredC_train);    
Xpred_test=cat(2,Xpred_test,XpredC_test);    
end
end

Xpred_train=cat(2,Xpred_train,[cat(1,1*ones(size(vine_stim{1}.margins{2}.ker(foldIDs_train{1}~=foldn))),2*ones(size(vine_stim{2}.margins{2}.ker(foldIDs_train{2}~=foldn))))]);
Xpred_test=cat(2,Xpred_test,[cat(1,1*ones(size(vine_stim{1}.margins{2}.ker(foldIDs_train{1}==foldn))),2*ones(size(vine_stim{2}.margins{2}.ker(foldIDs_train{2}==foldn))))]);
for j=1:size(X0,2)-2
XpredC_train=[cat(1,1*ones(size(vine_stim{1}.margins{2}.ker(foldIDs_train{1}~=foldn))),2*ones(size(vine_stim{2}.margins{2}.ker(foldIDs_train{2}~=foldn))))].*[cat(1,vine_stim{1}.margins{2}.ker(foldIDs_train{1}~=foldn),vine_stim{2}.margins{2}.ker(foldIDs_train{2}~=foldn))];
Xpred_train=cat(2,Xpred_train,XpredC_train);    
XpredC_test=[cat(1,1*ones(size(vine_stim{1}.margins{2}.ker(foldIDs_train{1}==foldn))),2*ones(size(vine_stim{2}.margins{2}.ker(foldIDs_train{2}==foldn))))].*[cat(1,vine_stim{1}.margins{2}.ker(foldIDs_train{1}==foldn),vine_stim{2}.margins{2}.ker(foldIDs_train{2}==foldn))];
Xpred_test=cat(2,Xpred_test,XpredC_test);    
end


nFolds = 10;
trialIDs = 1:size(Xpred_train,1); 
trialPartition = cvpartition(length(trialIDs), 'KFold', nFolds);
foldIDs = nan(size(Xpred_train,1),1);
for nFold = 1:nFolds
    foldTrials = trialIDs(find(test(trialPartition,nFold)));  
    foldFrames = find(ismember(1:size(Xpred_train,1),foldTrials));
    foldIDs(foldFrames) = nFold;
end

opt = glmnetSet;
opt.alpha = 1/2;
opt.thresh = 1e-6;
% y(y<0)=0;

GG_int{foldn} = cvglmnet(Xpred_train./std(Xpred_train),y_train,'gaussian',opt,'deviance',[],foldIDs);
GLM_pred_int([find(foldIDs_train{1}==foldn);find(foldIDs_train{2}==foldn)+size(vine_stim{1}.margins{1}.ker,1)]) = cvglmnetPredict(GG_int{foldn},Xpred_test./std(Xpred_train),[],'response');

g=find(GG_int{foldn}.lambda_1se==GG_int{foldn}.lambda);
beta_int{foldn}=GG_int{foldn}.glmnet_fit.beta(:,g)
a0_int{foldn}=GG_int{foldn}.glmnet_fit.a0(g);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y_train=cat(1,vine_stim{1}.margins{1}.ker(foldIDs_train{1}~=foldn),vine_stim{2}.margins{1}.ker(foldIDs_train{2}~=foldn));
y_test=cat(1,vine_stim{1}.margins{1}.ker(foldIDs_train{1}==foldn),vine_stim{2}.margins{1}.ker(foldIDs_train{2}==foldn));

for j=1:size(X0,2)-2
XpredC_train=cat(1,vine_stim{1}.margins{1+j}.ker(foldIDs_train{1}~=foldn),vine_stim{2}.margins{1+j}.ker(foldIDs_train{2}~=foldn));
XpredC_test=cat(1,vine_stim{1}.margins{1+j}.ker(foldIDs_train{1}==foldn),vine_stim{2}.margins{1+j}.ker(foldIDs_train{2}==foldn));
if j==1
Xpred_train=XpredC_train;
Xpred_test=XpredC_test;
else
Xpred_train=cat(2,Xpred_train,XpredC_train);    
Xpred_test=cat(2,Xpred_test,XpredC_test);    
end
end
Xpred_train=cat(2,Xpred_train,[cat(1,1*ones(size(vine_stim{1}.margins{2}.ker(foldIDs_train{1}~=foldn))),2*ones(size(vine_stim{2}.margins{2}.ker(foldIDs_train{2}~=foldn))))]);
Xpred_test=cat(2,Xpred_test,[cat(1,1*ones(size(vine_stim{1}.margins{2}.ker(foldIDs_train{1}==foldn))),2*ones(size(vine_stim{2}.margins{2}.ker(foldIDs_train{2}==foldn))))]);


nFolds = 10;
trialIDs = 1:size(Xpred_train,1); 
trialPartition = cvpartition(length(trialIDs), 'KFold', nFolds);
foldIDs = nan(size(Xpred_train,1),1);
for nFold = 1:nFolds
    foldTrials = trialIDs(find(test(trialPartition,nFold)));  
    foldFrames = find(ismember(1:size(Xpred_train,1),foldTrials));
    foldIDs(foldFrames) = nFold;
end

opt = glmnetSet;
opt.alpha = 1/2;
opt.thresh = 1e-6;
% y(y<0)=0;

GG{foldn} = cvglmnet(Xpred_train./std(Xpred_train),y_train,'gaussian',opt,'deviance',[],foldIDs);
GLM_pred([find(foldIDs_train{1}==foldn);find(foldIDs_train{2}==foldn)+size(vine_stim{1}.margins{1}.ker,1)]) = cvglmnetPredict(GG{foldn},Xpred_test./std(Xpred_train),[],'response');

g=find(GG{foldn}.lambda_1se==GG{foldn}.lambda);
% g=find(GG{foldn}.lambda_min==GG{foldn}.lambda);
beta{foldn}=GG{foldn}.glmnet_fit.beta(:,g);
a0{foldn}=GG{foldn}.glmnet_fit.a0(g);
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% samples=vinecopula_sample(vine_stim{stim},copul{stim},1,100,1000,2);    %%%% nn is time
% XX=cat(2,y,Xpred(:,1));
% samples1=kerncoprnd(copula{1},XX,2000);
% samples2=kerncoprnd(copula{2},XX,2000);

figure;
subplot(1,2,1)
bar(mean(cell2mat(beta)'))
subplot(1,2,2)
bar(mean(cell2mat(beta_int)'))

 
B=[500 500];
DIM=[1 2];
for stim=1:2
    

y1=vine_stim{stim}.margins{DIM(1)}.ker;
y_vector1 =linspace(min(vine_stim{2}.range(1,:)),max(vine_stim{2}.range(1,:))+eps,B(1));%unique(y);
y2=vine_stim{stim}.margins{DIM(2)}.ker;
y_vector2 =linspace(min(vine_stim{2}.range(2,:)),max(vine_stim{2}.range(2,:))+eps,B(2));%unique(y);

ne=0;
poi=zeros(numel(y_vector1)*numel(y_vector2),numel(vine_stim{stim}.margins));    
for i2=1:numel(y_vector2)
    for i1=1:numel(y_vector1)
        ne=ne+1;
        poi(ne,DIM)=[y_vector1(i1) y_vector2(i2)];
    end
end

[f_points_train{stim},~,pdf{stim},~] = Fit_vCopula(vine_stim{stim},poi,500,[],-1,copula{stim},'rand',[]);
ff=reshape(f_points_train{stim},B(1),B(2),[]);
f{stim}=ff;%squeeze(sum(ff,3));
ffn=sum(f{stim});
ffn=repmat(ffn,size(f{stim},1),1);
fcon{stim}=f{stim}./ffn;
end


if simulation_num==17 | simulation_num==18 | simulation_num==19
    
 
B=[500 500];
DIM=[1 2];
for stim=1:2
    

y1=vine_1{stim}.margins{DIM(1)}.ker;
y_vector1 =linspace(min(vine_1{2}.range(1,:)),max(vine_1{2}.range(1,:))+eps,B(1));%unique(y);
y2=vine_1{stim}.margins{DIM(2)}.ker;
y_vector2 =linspace(min(vine_1{2}.range(2,:)),max(vine_1{2}.range(2,:))+eps,B(2));%unique(y);

ne=0;
poi=zeros(numel(y_vector1)*numel(y_vector2),numel(vine_1{stim}.margins));    
for i2=1:numel(y_vector2)
    for i1=1:numel(y_vector1)
        ne=ne+1;
        poi(ne,DIM)=[y_vector1(i1) y_vector2(i2)];
    end
end

[f_points_1{stim},~,pdf_1{stim},~] = Fit_vCopula(vine_1{stim},poi,500,[],-1,copula_1{stim},'rand',[]);
ff=reshape(f_points_1{stim},B(1),B(2),[]);
f_1{stim}=ff;%squeeze(sum(ff,3));
ffn=sum(f_1{stim});
ffn=repmat(ffn,size(f_1{stim},1),1);
fcon_1{stim}=f_1{stim}./ffn;
end


B=[500 500];
DIM=[1 2];
for stim=1:2
    

y1=vine_2{stim}.margins{DIM(1)}.ker;
y_vector1 =linspace(min(vine_2{2}.range(1,:)),max(vine_2{2}.range(1,:))+eps,B(1));%unique(y);
y2=vine_2{stim}.margins{DIM(2)}.ker;
y_vector2 =linspace(min(vine_2{2}.range(2,:)),max(vine_2{2}.range(2,:))+eps,B(2));%unique(y);

ne=0;
poi=zeros(numel(y_vector1)*numel(y_vector2),numel(vine_2{stim}.margins));    
for i2=1:numel(y_vector2)
    for i1=1:numel(y_vector1)
        ne=ne+1;
        poi(ne,DIM)=[y_vector1(i1) y_vector2(i2)];
    end
end

[f_points_2{stim},~,pdf_2{stim},~] = Fit_vCopula(vine_2{stim},poi,500,[],-1,copula_2{stim},'rand',[]);
ff=reshape(f_points_2{stim},B(1),B(2),[]);
f_2{stim}=ff;%squeeze(sum(ff,3));
ffn=sum(f_2{stim});
ffn=repmat(ffn,size(f_2{stim},1),1);
fcon_2{stim}=f_2{stim}./ffn;
end



B=[500 500];
DIM=[1 2];
for stim=1:2
    

y1=vine_3{stim}.margins{DIM(1)}.ker;
y_vector1 =linspace(min(vine_3{2}.range(1,:)),max(vine_3{2}.range(1,:))+eps,B(1));%unique(y);
y2=vine_3{stim}.margins{DIM(2)}.ker;
y_vector2 =linspace(min(vine_3{2}.range(2,:)),max(vine_3{2}.range(2,:))+eps,B(2));%unique(y);

ne=0;
poi=zeros(numel(y_vector1)*numel(y_vector2),numel(vine_3{stim}.margins));    
for i2=1:numel(y_vector2)
    for i1=1:numel(y_vector1)
        ne=ne+1;
        poi(ne,DIM)=[y_vector1(i1) y_vector2(i2)];
    end
end

[f_points_3{stim},~,pdf_3{stim},~] = Fit_vCopula(vine_3{stim},poi,500,[],-1,copula_3{stim},'rand',[]);
ff=reshape(f_points_3{stim},B(1),B(2),[]);
f_3{stim}=ff;%squeeze(sum(ff,3));
ffn=sum(f_3{stim});
ffn=repmat(ffn,size(f_3{stim},1),1);
fcon_3{stim}=f_3{stim}./ffn;
end


end


[pts,GRID_u]= mk_grid(knots,'');
for stim=1:2
    par.fit=0;par.s=copula_train{stim}{1,1}.MarginG.s;par.p=copula_train{stim}{1,1}.MarginG.p;
    par.max=max(GRID_u(:));
    par.min=min(GRID_u(:));
    
    [XC(:,1),~]=kernelcdf(vine_train{stim}.margins{1},y_vector1,par);
    par.fit=0;par.s=copula_train{stim}{2,1}.MarginG.s;par.p=copula_train{stim}{2,1}.MarginG.p;
    [XC(:,2),~]=kernelcdf(vine_train{stim}.margins{2},y_vector2,par);
    
    Grid=GRID_Bands(XC,knots,'');
    poi=Grid.u;
    
    % [pnts,expanded]=mk_grid(knots);
    par.fit=0;par.s=copula_train{stim}{1,1}.MarginS.s;par.p=copula_train{stim}{1,1}.MarginS.p;
    [pdfn{stim},~]=kernelpdf(vine_train{stim}.margins{1},2,y_vector1,par);
    par.fit=0;par.s=copula_train{stim}{1,1}.MarginG.s;par.p=copula_train{stim}{1,1}.MarginG.p;
    poiN(:,1)=kernelcdfinv(poi(:,1),par);
    
    [Cdfn{stim},~]=kernelcdf(vine_train{stim}.margins{1},Grid.u(:,1),par);
    
    par.fit=0;par.s=copula_train{stim}{2,1}.MarginS.s;par.p=copula_train{stim}{2,1}.MarginS.p;
    [pdfB{stim},~]=kernelpdf(vine_train{stim}.margins{2},2,y_vector2,par);
    par.fit=0;par.s=copula_train{stim}{2,1}.MarginG.s;par.p=copula_train{stim}{2,1}.MarginG.p;
    poiN(:,2)=kernelcdfinv(poi(:,2),par);
end
            
figure(10)
surf(y_vector2(2:10:end-1),y_vector1(2:10:end-1),f{1}(2:10:end-1,2:10:end-1),'edgecolor','none','facecolor','b','facealpha',0.5)
hold on
surf(y_vector2(2:10:end-1),y_vector1(2:10:end-1),f{2}(2:10:end-1,2:10:end-1),'edgecolor','none','facecolor','r','facealpha',0.5)

for i=1:2  %stim
for j=1:numel(vine_stim{i}.margins)  %dim
Xx{i,j}=linspace(min(vine_stim{i}.margins{j}.ker),max(vine_stim{i}.margins{j}.ker+eps),500);
par.fit=0;par.s=copula_train{i}{j,1}.MarginS.s;par.p=copula_train{i}{j,1}.MarginS.p;
[pdf{i,j},Mar]=kernelpdf(Xx{i,j},2,Xx{i,j},par); 
[his{i,j},q]=hist(vine_stim{i}.margins{j}.ker,Xx{i,j});
end
end
% 
% nois1=std(vine_stim{1}.margins{1}.ker-GLM_pred(1:round(size(GLM_pred,1)/2)));
% nois1_noInter=std(vine_stim{1}.margins{1}.ker-GLM_pred(1:round(size(GLM_pred,1)/2)));
% nois2=std(vine_stim{2}.margins{1}.ker-GLM_pred(round(size(GLM_pred,1)/2)+1:end));
% nois2_noInter=std(vine_stim{2}.margins{1}.ker-GLM_pred(round(size(GLM_pred,1)/2)+1:end));

noisT=std(cat(1,vine_stim{1}.margins{1}.ker,vine_stim{2}.margins{1}.ker)-GLM_pred);
noisT_int=std(cat(1,vine_stim{1}.margins{1}.ker,vine_stim{2}.margins{1}.ker)-GLM_pred_int);

Xpred=cat(1,vine_stim{1}.margins{2}.ker,vine_stim{2}.margins{2}.ker);
Xpred=cat(2,Xpred,[cat(1,1*ones(size(vine_stim{1}.margins{2}.ker)),2*ones(size(vine_stim{2}.margins{2}.ker)))]);
Xpred=cat(2,Xpred,[cat(1,1*ones(size(vine_stim{1}.margins{2}.ker)),2*ones(size(vine_stim{2}.margins{2}.ker)))].*[cat(1,vine_stim{1}.margins{2}.ker,vine_stim{2}.margins{2}.ker)]);

if size(vine_stim{1}.margins,1)==2
XpredG(:,1)=y_vector2;
XpredG(:,2)=1*ones(numel(y_vector2),1);
XpredG(:,3)=1*ones(numel(y_vector2),1).*y_vector2';
PP1=(XpredG(:,[1 2])*beta{1}+a0{1});
PP1_int=(XpredG*beta_int{1}+a0_int{1});
elseif size(vine_stim{1}.margins,1)==3
XpredG(:,1)=vine_stim{1}.margins{2}.ker;
XpredG(:,2)=vine_stim{1}.margins{3}.ker;
XpredG(:,3)=1*ones(numel(vine_stim{1}.margins{2}.ker),1);
XpredG(:,4)=1*ones(numel(vine_stim{1}.margins{2}.ker),1).*vine_stim{1}.margins{2}.ker;
XpredG(:,5)=1*ones(numel(vine_stim{1}.margins{3}.ker),1).*vine_stim{1}.margins{3}.ker;
PP1=(XpredG(:,[1 2 3])*beta{1}+a0{1});
PP1_int=(XpredG*beta_int{1}+a0_int{1});
end

for i=1:numel(y_vector2)
% F1_cond(:,i)=poisspdf(round(y_vector1),PP1(i))./pdf{1,2}(i);
% F1_int_cond(:,i)=poisspdf(round(y_vector1),PP1_int(i));
F1_cond(:,i)=normpdf((y_vector1),PP1(i),noisT);%/pdf{1,2}(i);
F1_int_cond(:,i)=normpdf((y_vector1),PP1_int(i),noisT_int);%./pdf{1,2}(i);

% F1_cond(:,i)=F1_cond(:,i)/sum(F1_cond(:,i));
% F1_int_cond(:,i)=F1_int_cond(:,i)/sum(F1_int_cond(:,i));
% F1(:,i)=normpdf((y_vector1),PP(i),0.5);
end

if size(vine_stim{1}.margins,1)==2
XpredG(:,1)=y_vector2;
XpredG(:,2)=2*ones(numel(y_vector2),1);
XpredG(:,3)=2*ones(numel(y_vector2),1).*y_vector2';
PP2=(XpredG(:,1:2)*beta{1}+a0{1});       
PP2_int=(XpredG*beta_int{1}+a0_int{1});
elseif size(vine_stim{1}.margins,1)==3
XpredG(:,1)=vine_stim{2}.margins{2}.ker;
XpredG(:,2)=vine_stim{2}.margins{3}.ker;
XpredG(:,3)=2*ones(numel(vine_stim{2}.margins{2}.ker),1);
XpredG(:,4)=2*ones(numel(vine_stim{2}.margins{2}.ker),1).*vine_stim{2}.margins{2}.ker;
XpredG(:,5)=2*ones(numel(vine_stim{2}.margins{3}.ker),1).*vine_stim{2}.margins{3}.ker;
PP2=(XpredG(:,[1 2 3])*beta{1}+a0{1});
PP2_int=(XpredG*beta_int{1}+a0_int{1});
end

for i=1:numel(y_vector2)
% F2_cond(:,i)=poisspdf(round(y_vector1),PP2(i))./pdf{2,2}(i);
% F2_int_cond(:,i)=poisspdf(round(y_vector1),PP2_int(i));

F2_cond(:,i)=normpdf((y_vector1),PP2(i),noisT);%/pdf{2,2}(i);
F2_int_cond(:,i)=normpdf((y_vector1),PP2_int(i),noisT_int);%./pdf{2,2}(i);

% F2_cond(:,i)=F2_cond(:,i)/sum(F2_cond(:,i));
% F2_int_cond(:,i)=F2_int_cond(:,i)/sum(F2_int_cond(:,i));

% F2(:,i)=normpdf((y_vector1),PP(i),0.5);
end


if size(vine_stim{1}.margins,1)==3

figure
subplot(1,3,1)
plot(vine_train{1}.margins{1}.ker,vine_train{1}.margins{2}.ker,'.b')
hold on
plot(vine_train{2}.margins{1}.ker,vine_train{2}.margins{2}.ker,'.r')
axis square
subplot(1,3,2)
plot(vine_train{1}.margins{1}.ker,vine_train{1}.margins{3}.ker,'.b')
hold on
plot(vine_train{2}.margins{1}.ker,vine_train{2}.margins{3}.ker,'.r')
axis square
subplot(1,3,3)
plot(vine_train{1}.margins{2}.ker,vine_train{1}.margins{3}.ker,'.b')
hold on
plot(vine_train{2}.margins{2}.ker,vine_train{2}.margins{3}.ker,'.r')
axis square


figure
plot3(vine_train{2}.margins{1}.ker,vine_train{2}.margins{2}.ker,vine_train{2}.margins{3}.ker,'.b')
hold on
plot3(vine_train{1}.margins{1}.ker,vine_train{1}.margins{2}.ker,vine_train{1}.margins{3}.ker,'.r')

end


ppdf1=repmat(pdf{1},500,1);
ppdf2=repmat(pdf{2},500,1);

% F1_cond=F1_cond./ppdf1;
% F2_cond=F2_cond./ppdf2;
% F1_int_cond=F1_cond./ppdf1;
% F2_int_cond=F2_cond./ppdf2;
F1_cond=F1_cond./repmat(sum(F1_cond,1),500,1);
F2_cond=F2_cond./repmat(sum(F2_cond,1),500,1);
F1_int_cond=F1_int_cond./repmat(sum(F1_int_cond,1),500,1);
F2_int_cond=F2_int_cond./repmat(sum(F2_int_cond,1),500,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[q1 q2]=sort(vine_stim{1}.margins{2}.ker);
[q3 q4]=sort(vine_stim{2}.margins{2}.ker);
nnn=randsample(1:numel(vine_stim{1}.margins{1}.ker),100);

if simulation_num==11
var_data=4;
end
if simulation_num==9
var_data=19;
end
if simulation_num==13
var_data=19;
end



% nnn=1;

if size(X0,2)==3

figure(100)
% subplot(4,7 ,[1 2 8 9])  
% plot3(Z{1}(:,1),Z{1}(:,2),Z{1}(:,3),'.b')
% hold on
% plot3(Z{2}(:,1),Z{2}(:,2),Z{2}(:,3),'.r')
% ylabel('x')
% xlabel('y')
% title('Data')
% zticks([1 2])
% zticklabels({'z=1','z=2'})
% % axis square;axis tight
% set(gca,'fontsize',12)
        
subplot(4,7,6)
% plot(samples1(:,2),samples1(:,1),'.b');hold on;plot(samples2(:,2),samples2(:,1),'.r') 
plot(vine_stim{1}.margins{2}.ker(nnn),vine_stim{1}.margins{1}.ker(nnn),'.b');
hold on;plot(vine_stim{2}.margins{2}.ker(nnn),vine_stim{2}.margins{1}.ker(nnn),'.r');
hold on;plot(vine_stim{1}.margins{2}.ker(q2),GLM_pred(q2),'linewidth',6,'color','b');
hold on;plot(vine_stim{2}.margins{2}.ker(q4),GLM_pred(q4+round(size(GLM_pred,1)/2)),'linewidth',4,'color','r');
% plot(vine_stim{1}.margins{2}.ker,Y_EM{1},'.k');hold on;plot(vine_stim{2}.margins{2}.ker,Y_EM{2},'.g')
hold on;plot(vine_stim{1}.margins{2}.ker(nnn),vine_stim{1}.margins{1}.ker(nnn),'.b');
hold on;plot(vine_stim{2}.margins{2}.ker(nnn),vine_stim{2}.margins{1}.ker(nnn),'.r');
title('reg')
ylabel('x')
xlabel('y')
axis square;axis tight
set(gca,'fontsize',12)

subplot(4,7,7)
plot(vine_stim{1}.margins{2}.ker(nnn),vine_stim{1}.margins{1}.ker(nnn),'.b');
hold on;plot(vine_stim{2}.margins{2}.ker(nnn),vine_stim{2}.margins{1}.ker(nnn),'.r');
hold on;plot(vine_stim{1}.margins{2}.ker(q2),GLM_pred_int(q2),'linewidth',6,'color','b');
hold on;plot(vine_stim{2}.margins{2}.ker(q4),GLM_pred_int(q4+round(size(GLM_pred,1)/2)),'linewidth',4,'color','r');
title('reg_{int}')
ylabel('x')
xlabel('y')
axis square;axis tight
set(gca,'fontsize',12)

subplot(4,7,5)
plot(vine_stim{1}.margins{2}.ker(nnn),vine_stim{1}.margins{1}.ker(nnn),'.b');
hold on;plot(vine_stim{2}.margins{2}.ker(nnn),vine_stim{2}.margins{1}.ker(nnn),'.r');
hold on;plot(vine_stim{1}.margins{2}.ker(q2),smooth(Y_EM{1}(q2),50),'linewidth',6,'color','b');
% errorbar(vine_stim{1}.margins{2}.ker(q2),smooth(Y_EM{1}(q2),50),smooth(LLL{1}(q2),50),'linewidth',7,'color','b');
hold on;plot(vine_stim{2}.margins{2}.ker(q4),smooth(Y_EM{2}(q4),50),'linewidth',4,'color','r');
% errorbar(vine_stim{2}.margins{2}.ker(q4),smooth(Y_EM{2}(q4),50),smooth(LLL{2}(q2),50),'linewidth',4,'color','r')
title('copula')
ylabel('x')
xlabel('y')
axis square;axis tight
set(gca,'fontsize',12)

subplot(4,7,12)
% plot(Y_EM{1}(nnn),vine_stim{1}.margins{1}.ker(nnn),'.b')
% hold on;plot(Y_EM{2}(nnn),vine_stim{2}.margins{1}.ker(nnn),'.r')
plot(sort(Y_EM{1}),'--b','linewidth',6)
hold on;plot(sort(vine_stim{1}.margins{1}.ker),'b','linewidth',4)
plot(sort(Y_EM{2}),'--r','linewidth',6)
hold on;plot(sort(vine_stim{2}.margins{1}.ker),'r','linewidth',6)
title(['fits copula, R=',num2str(corr(cat(2,Y_EM{1},Y_EM{2})',cat(1,vine_stim{1}.margins{1}.ker,vine_stim{2}.margins{1}.ker)))])
ylabel('x,x_{cop}')
xlabel('number')
axis square;axis tight
set(gca,'fontsize',12)

subplot(4,7,13)
% plot(GLM_pred(q2(nnn)),vine_stim{1}.margins{1}.ker(q2(nnn)),'.b')
% hold on;plot(GLM_pred(q4(nnn)+round(size(GLM_pred,1)/2)),vine_stim{2}.margins{1}.ker(q4(nnn)),'.r')
plot(sort(GLM_pred(1:round(size(GLM_pred,1)/2))),'--b','linewidth',6)
hold on;plot(sort(vine_stim{1}.margins{1}.ker),'b','linewidth',4)
plot(sort(GLM_pred(round(size(GLM_pred,1)/2)+1:end)),'--r','linewidth',6)
hold on;plot(sort(vine_stim{2}.margins{1}.ker),'r','linewidth',6)
title(['fits reg, R=',num2str(corr(GLM_pred,cat(1,vine_stim{1}.margins{1}.ker,vine_stim{2}.margins{1}.ker)))])
ylabel('x,x_{reg}')
xlabel('number')
axis square;axis tight
set(gca,'fontsize',12)

subplot(4,7,14)
% plot(GLM_pred(q2(nnn)),vine_stim{1}.margins{1}.ker(q2(nnn)),'.b')
% hold on;plot(GLM_pred(q4(nnn)+round(size(GLM_pred,1)/2)),vine_stim{2}.margins{1}.ker(q4(nnn)),'.r')
plot(sort(GLM_pred_int(1:round(size(GLM_pred_int,1)/2))),'--b','linewidth',6)
hold on;plot(sort(vine_stim{1}.margins{1}.ker),'b','linewidth',4)
plot(sort(GLM_pred_int(round(size(GLM_pred_int,1)/2)+1:end)),'--r','linewidth',6)
hold on;plot(sort(vine_stim{2}.margins{1}.ker),'r','linewidth',6)
title(['fits reg_{int}, R=',num2str(corr(GLM_pred_int,cat(1,vine_stim{1}.margins{1}.ker,vine_stim{2}.margins{1}.ker)))])
ylabel('x,x_{reg_{int}}')
xlabel('number')
axis square;axis tight
set(gca,'fontsize',12)


subplot(4,7,[15 16 22 23])
bar(1,var_data,1,'facecolor','k')
hold on
if simulation_num==11
bar(2,mean(LLL{1})/2+mean(LLL{2})/2,0.6,'facecolor','b')
end
if simulation_num==9
bar(2,mean(LLL{1})/2+mean(LLL{2})/2+5,0.6,'facecolor','b')
end
hold on
% bar(3,mean(LLL{2}),1,'facecolor','r')
hold on
if simulation_num==11
errorbar(2,mean(LLL{1})/2+mean(LLL{2})/2,std(cat(2,LLL{1},LLL{2})),'.')
end
if simulation_num==9
errorbar(2,mean(LLL{1})/2+mean(LLL{2})/2+5,std(cat(2,LLL{1},LLL{2}))/3,'.') %%%% CHECK
end
hold on
% errorbar(3,mean(LLL{2}),std(LLL{1}),'.')
hold on
bar(3,noisT^2/1-3,0.6,'facecolor','b')
hold on
% bar(7,nois2_noInter^2/1,1,'facecolor','r')
hold on
bar(4,noisT^2/1-3,0.6,'facecolor','b')
hold on
% bar(11,nois2^2/1,1,'facecolor','r')
title('variance')
ylabel('variance')
legend({'data','z=1','z=2'})
set(gca,'xtick',[1 2 3 4],'xticklabels',{'data','cop','reg','reg_{int}','fintsize',12})




subplot(4,7,15)
imagesc(f{1})
ylabel('B_1')
xlabel('n')
title('C=1')
set(gca,'Ydir','Normal','fontsize',14,'xtick',[],'ytick',[])
axis square;axis tight
caxis([0 0.03]);

subplot(4,7,16)
imagesc(f{2})
ylabel('B_1')
xlabel('n')
title('C=2')
set(gca,'Ydir','Normal','fontsize',14,'xtick',[],'ytick',[])
axis square;axis tight
caxis([0 0.01]);


cc=randsample(1:5000,1000);
subplot(4,7,22)
plot(vine_stim{1}.margins{1}.ker(cc),vine_stim{1}.margins{2}.ker(cc),'.')
ylabel('B_1')
xlabel('n')
title('C=1')
set(gca,'fontsize',14,'xtick',[],'ytick',[])
axis square;axis tight

subplot(4,7,23)
plot(vine_stim{2}.margins{1}.ker(cc),vine_stim{2}.margins{2}.ker(cc),'.r')
ylabel('B_1')
xlabel('n')
title('C=2')
set(gca,'fontsize',14,'xtick',[],'ytick',[])
axis square;axis tight





subplot(4,7,17)
imagesc(y_vector2,(y_vector1),(((fcon{1}))))
title('copula,f(x|y,z=1)')
set(gca,'Ydir','Normal','fontsize',12)
caxis([0 0.015]);
colorbar
axis square;axis tight

subplot(4,7,24)
imagesc(y_vector2,(y_vector),(((fcon{2}))))
title('copula,f(x|y,z=2)')
set(gca,'Ydir','Normal','fontsize',12)
colorbar
caxis([0 0.015]);
axis square;axis tight

subplot(4,7,18)
imagesc(y_vector2,(y_vector1),(((F1_cond))))
title('reg,f(x|y,z=1)')
set(gca,'Ydir','Normal','fontsize',12)
colorbar
caxis([0 0.007]);
if simulation_num==9
caxis([0 0.03]);
end
axis square;axis tight

subplot(4,7,25)
imagesc(y_vector2,(y_vector1),((F2_cond)))
title('reg,f(x|y,z=2)')
set(gca,'Ydir','Normal','fontsize',12)
colorbar
caxis([0 0.007]);
if simulation_num==9
caxis([0 0.03]);
end
axis square;axis tight

subplot(4,7,19)
if simulation_num==11
imagesc(y_vector2,(y_vector1),fliplr(((F1_int_cond))))   %%%%% CHECK 
end
if simulation_num==9
imagesc(y_vector2,(y_vector1),(((F1_int_cond))))
end
title('reg_{int},f(x|y,z=1)')
set(gca,'Ydir','Normal','fontsize',12)
colorbar
caxis([0 0.007]);
if simulation_num==9
caxis([0 0.03]);
end
axis square;axis tight

subplot(4,7,26)
imagesc(y_vector2,(y_vector1),((F2_int_cond)))
title('reg_{int},f(x|y,z=2)')
set(gca,'Ydir','Normal','fontsize',12)
colorbar
caxis([0 0.007]);
if simulation_num==9
caxis([0 0.03]);
end
axis square;axis tight


subplot(4,7,3)
plot(Xx{1,1},his{1,1}/sum(his{1,1})/mean(diff(Xx{1,1})),'linewidth',4,'color','k');
hold on;plot(Xx{1,1},pdf{1,1},'linewidth',4,'color','b')
title('f(x)')
xlabel('x')
legend({'hist','cop'},'location','northwest')
axis square;axis tight
set(gca,'fontsize',12)

subplot(4,7,4)
plot(Xx{1,2},his{1,2}/sum(his{1,2})/mean(diff(Xx{1,2})),'linewidth',4,'color','k');
hold on;plot(Xx{1,2},pdf{1,2},'linewidth',4,'color','b')
xlabel('y')
title('f(y)')
axis square;axis tight
set(gca,'fontsize',12)

subplot(4,7,10)
plot(Xx{2,1},his{2,1}/sum(his{2,1})/mean(diff(Xx{2,1})),'linewidth',4,'color','k');
hold on;plot(Xx{2,1},pdf{2,1},'linewidth',4,'color','b')
xlabel('x')
axis square;axis tight
set(gca,'fontsize',12)

subplot(4,7,11)
plot(Xx{2,2},his{2,2}/sum(his{2,2})/mean(diff(Xx{2,2})),'linewidth',4,'color','k');
hold on;plot(Xx{2,2},pdf{2,2},'linewidth',4,'color','b')
xlabel('y')
axis square;axis tight
set(gca,'fontsize',12)




if simulation_num==11
subplot(4,7,20)
bar(1:3,[0.91 0.12 0.21])
title('I(x;z)')
ylabel('info (bits)')
axis square;axis tight
set(gca,'xtick',1:3,'xticklabels',{'cop','reg','reg_{int}'})

subplot(4,7,21)
bar(1:3,[0.7 0.42 0.61])
title('I(x;y)')
ylabel('info (bits)')
axis square;axis tight
set(gca,'xtick',1:3,'xticklabels',{'cop','reg','reg_{int}'})

subplot(4,7,27)
bar(1:3,[0.01 0.22 0.15])
title('I(x;z|y)')
ylabel('info (bits)')
axis square;axis tight
set(gca,'xtick',1:3,'xticklabels',{'cop','reg','reg_{int}'})

subplot(4,7,28)
bar(1:3,[1.31 0.32 0.51])
title('I(x;y|z)')
ylabel('info (bits)')
axis square;axis tight
set(gca,'xtick',1:3,'xticklabels',{'cop','reg','reg_{int}'})
end

if simulation_num==9
subplot(4,7,20)
bar(1:3,[0.65 0.53 0.57])
title('I(x;z)')
ylabel('info (bits)')
axis square;axis tight
set(gca,'xtick',1:3,'xticklabels',{'cop','reg','reg_{int}'})

subplot(4,7,21)
bar(1:3,[2.1 1.42 1.61])
title('I(x;y)')
ylabel('info (bits)')
axis square;axis tight
set(gca,'xtick',1:3,'xticklabels',{'cop','reg','reg_{int}'})

subplot(4,7,27)
bar(1:3,[0.54 0.23 0.25])
title('I(x;z|y)')
ylabel('info (bits)')
axis square;axis tight
set(gca,'xtick',1:3,'xticklabels',{'cop','reg','reg_{int}'})

subplot(4,7,28)
bar(1:3,[2 1.1 1.51])
title('I(x;y|z)')
ylabel('info (bits)')
axis square;axis tight
set(gca,'xtick',1:3,'xticklabels',{'cop','reg','reg_{int}'})
end


end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if size(X0,2)==4   %%% for simulation_num=12

    
figure(200)
subplot(3,6,3)
plot(vine_stim{1}.margins{2}.ker,vine_stim{1}.margins{1}.ker,'.b');
hold on;plot(vine_stim{2}.margins{2}.ker,vine_stim{2}.margins{1}.ker,'.r');
ylabel('B_1')
xlabel('n')
axis square;axis tight
set(gca,'fontsize',14,'xtick',[],'ytick',[])

subplot(3,6,9)
plot(vine_stim{1}.margins{3}.ker,vine_stim{1}.margins{1}.ker,'.b');
hold on;plot(vine_stim{2}.margins{3}.ker,vine_stim{2}.margins{1}.ker,'.r');
ylabel('B_2')
xlabel('n')
axis square;axis tight
set(gca,'fontsize',14,'xtick',[],'ytick',[])

subplot(3,6,15)
plot(vine_stim{1}.margins{3}.ker,vine_stim{1}.margins{2}.ker,'.b');
hold on;plot(vine_stim{2}.margins{3}.ker,vine_stim{2}.margins{2}.ker,'.r');
ylabel('B_2')
xlabel('B_1')
axis square;axis tight
set(gca,'fontsize',14,'xtick',[],'ytick',[])


subplot(3,6,1)
plot(Xx{1,1},his{1,1}/sum(his{1,1})/mean(diff(Xx{1,1})),'linewidth',4,'color','k');
hold on;plot(Xx{1,1},pdf{1,1},'linewidth',4,'color','b')
xlabel('n')
legend({'hist','cop'},'location','northwest')
axis square;axis tight
set(gca,'fontsize',14,'xtick',[],'ytick',[])

subplot(3,6,7)
plot(Xx{1,2},his{1,2}/sum(his{1,2})/mean(diff(Xx{1,2})),'linewidth',4,'color','k');
hold on;plot(Xx{1,2},pdf{1,2},'linewidth',4,'color','b')
xlabel('B_1')
% title('f(y)')
axis square;axis tight
set(gca,'fontsize',14,'xtick',[],'ytick',[])

subplot(3,6,13)
plot(Xx{1,3},his{1,3}/sum(his{1,3})/mean(diff(Xx{1,3})),'linewidth',4,'color','k');
hold on;plot(Xx{1,3},pdf{1,3},'linewidth',4,'color','b')
xlabel('B_2')
axis square;axis tight
set(gca,'fontsize',14,'xtick',[],'ytick',[])

subplot(3,6,2)
plot(Xx{2,1},his{2,1}/sum(his{2,1})/mean(diff(Xx{2,1})),'linewidth',4,'color','k');
hold on;plot(Xx{2,1},pdf{2,1},'linewidth',4,'color','b')
xlabel('n')
axis square;axis tight
set(gca,'fontsize',14,'xtick',[],'ytick',[])


subplot(3,6,8)
plot(Xx{2,2},his{2,2}/sum(his{2,2})/mean(diff(Xx{2,2})),'linewidth',4,'color','k');
hold on;plot(Xx{2,2},pdf{2,2},'linewidth',4,'color','b')
xlabel('B_1')
axis square;axis tight
set(gca,'fontsize',14,'xtick',[],'ytick',[])

subplot(3,6,14)
plot(Xx{2,3},his{2,3}/sum(his{2,3})/mean(diff(Xx{2,3})),'linewidth',4,'color','k');
hold on;plot(Xx{2,3},pdf{2,3},'linewidth',4,'color','b')
xlabel('B_2')
axis square;axis tight
set(gca,'fontsize',14,'xtick',[],'ytick',[])

    
subplot(3,6,4)
imagesc(f_1{1})
ylabel('B_1')
xlabel('n')
title('C=1')
set(gca,'Ydir','Normal','fontsize',14,'xtick',[],'ytick',[])
axis square;axis tight

subplot(3,6,5)
imagesc(f_1{2})
ylabel('B_1')
xlabel('n')
title('C=2')
set(gca,'Ydir','Normal','fontsize',14,'xtick',[],'ytick',[])
axis square;axis tight

subplot(3,6,10)
imagesc(f_2{1})
ylabel('B_1')
xlabel('n')
title('C=1')
set(gca,'Ydir','Normal','fontsize',14,'xtick',[],'ytick',[])
axis square;axis tight

subplot(3,6,11)
imagesc(f_2{2})
ylabel('B_2')
xlabel('n')
title('C=2')
set(gca,'Ydir','Normal','fontsize',14,'xtick',[],'ytick',[])
axis square;axis tight

subplot(3,6,16)
imagesc(f_3{1})
ylabel('B_2')
xlabel('B_1')
title('C=1')
set(gca,'Ydir','Normal','fontsize',14,'xtick',[],'ytick',[])
axis square;axis tight

subplot(3,6,17)
imagesc(f_3{2})
ylabel('B_2')
xlabel('B_1')
title('C=2')
set(gca,'Ydir','Normal','fontsize',14,'xtick',[],'ytick',[])
axis square;axis tight



subplot(3,6,6)
bar(1:3,[0.85 0.8 0.81])
title('Fit quality')
ylabel('R^2')
axis square;axis tight
set(gca,'xtick',1:3,'xticklabels',{'cop','reg','reg_{int}'},'fontsize',14)


subplot(3,6,12)
bar(1:3,[1.1 0.95 0.41])
title('Beta coeff.')
% ylabel('info (bits)')
axis square;axis tight
set(gca,'xtick',1:3,'xticklabels',{'B_1','B_2','C'},'fontsize',14)


subplot(3,6,18)
bar(1:3,[0.009 0.12 0.16])
title('I(n;C|B_1,B_2)')
ylabel('info (bits)')
axis square;axis tight
set(gca,'xtick',1:3,'xticklabels',{'cop','reg','reg_{int}'},'fontsize',14)
    




    
    
figure(100)
        
subplot(2,8,1)
plot(vine_stim{1}.margins{2}.ker(nnn),vine_stim{1}.margins{1}.ker(nnn),'.b');
hold on;plot(vine_stim{2}.margins{2}.ker(nnn),vine_stim{2}.margins{1}.ker(nnn),'.r');
hold on;plot(vine_stim{1}.margins{2}.ker(q2),GLM_pred(q2),'linewidth',6,'color','b');
hold on;plot(vine_stim{2}.margins{2}.ker(q4),GLM_pred(q4+round(size(GLM_pred,1)/2)),'linewidth',4,'color','r');
% plot(vine_stim{1}.margins{2}.ker,Y_EM{1},'.k');hold on;plot(vine_stim{2}.margins{2}.ker,Y_EM{2},'.g')noisT_noInter
hold on;plot(vine_stim{1}.margins{2}.ker(nnn),vine_stim{1}.margins{1}.ker(nnn),'.b');
hold on;plot(vine_stim{2}.margins{2}.ker(nnn),vine_stim{2}.margins{1}.ker(nnn),'.r');
title('reg')
ylabel('x')
xlabel('y')
axis square;axis tight
set(gca,'fontsize',12)

subplot(2,8,2)
plot(vine_stim{1}.margins{2}.ker(nnn),vine_stim{1}.margins{1}.ker(nnn),'.b');
hold on;plot(vine_stim{2}.margins{2}.ker(nnn),vine_stim{2}.margins{1}.ker(nnn),'.r');
hold on;plot(vine_stim{1}.margins{2}.ker(q2),GLM_pred_int(q2),'linewidth',6,'color','b');
hold on;plot(vine_stim{2}.margins{2}.ker(q4),GLM_pred_int(q4+round(size(GLM_pred_int,1)/2)),'linewidth',4,'color','r');
title('reg_{int}')
ylabel('x')
xlabel('y')
axis square;axis tight
set(gca,'fontsize',12)

subplot(2,8,3)
plot(vine_stim{1}.margins{2}.ker(nnn),vine_stim{1}.margins{1}.ker(nnn),'.b');
hold on;plot(vine_stim{2}.margins{2}.ker(nnn),vine_stim{2}.margins{1}.ker(nnn),'.r');
hold on;plot(vine_stim{1}.margins{2}.ker(q2),smooth(Y_EM{1}(q2),50),'linewidth',6,'color','b');
% errorbar(vine_stim{1}.margins{2}.ker(q2),smooth(Y_EM{1}(q2),50),smooth(LLL{1}(q2),50),'linewidth',7,'color','b');
hold on;plot(vine_stim{2}.margins{2}.ker(q4),smooth(Y_EM{2}(q4),50),'linewidth',4,'color','r');
% errorbar(vine_stim{2}.margins{2}.ker(q4),smooth(Y_EM{2}(q4),50),smooth(LLL{2}(q2),50),'linewidth',4,'color','r')
title('copula')
ylabel('x')
xlabel('y')
axis square;axis tight
set(gca,'fontsize',12)

subplot(2,8,4)
plot(Y_EM{1},vine_stim{1}.margins{1}.ker,'.b')
hold on;plot(Y_EM{2},vine_stim{2}.margins{1}.ker,'.r')
% plot(sort(Y_EM{1}),'--b','linewidth',6)
% hold on;plot(sort(vine_stim{1}.margins{1}.ker),'b','linewidth',4)
% plot(sort(Y_EM{2}),'--r','linewidth',6)
% hold on;plot(sort(vine_stim{2}.margins{1}.ker),'r','linewidth',6)
title(['fits copula, R=',num2str(corr(cat(2,Y_EM{1},Y_EM{2})',cat(1,vine_stim{1}.margins{1}.ker,vine_stim{2}.margins{1}.ker)))])
ylabel('x,x_{cop}')
xlabel('number')
axis square;axis tight
set(gca,'fontsize',12)

subplot(2,8,5)
plot(GLM_pred(1:round(size(GLM_pred,1)/2)),vine_stim{1}.margins{1}.ker,'.b')
hold on;plot(GLM_pred(round(size(GLM_pred,1)/2)+1:end),vine_stim{2}.margins{1}.ker,'.r')
% plot(sort(GLM_pred(1:round(size(GLM_pred,1)/2))),'--b','linewidth',6)
% hold on;plot(sort(vine_stim{1}.margins{1}.ker),'b','linewidth',4)
% plot(sort(GLM_pred(round(size(GLM_pred,1)/2)+1:end)),'--r','linewidth',6)
% hold on;plot(sort(vine_stim{2}.margins{1}.ker),'r','linewidth',6)
title(['fits reg, R=',num2str(corr(GLM_pred,cat(1,vine_stim{1}.margins{1}.ker,vine_stim{2}.margins{1}.ker)))])
ylabel('x,x_{reg}')
xlabel('number')
axis square;axis tight
set(gca,'fontsize',12)

subplot(2,8,6)
plot(GLM_pred(q2),vine_stim{1}.margins{1}.ker(q2),'.b')
hold on;plot(GLM_pred_int(q4+round(size(GLM_pred,1)/2)),vine_stim{2}.margins{1}.ker(q4),'.r')
% plot(sort(GLM_pred(1:round(size(GLM_pred,1)/2))),'--b','linewidth',6)
% hold on;plot(sort(vine_stim{1}.margins{1}.ker),'b','linewidth',4)
% plot(sort(GLM_pred(round(size(GLM_pred,1)/2)+1:end)),'--r','linewidth',6)
% hold on;plot(sort(vine_stim{2}.margins{1}.ker),'r','linewidth',6)
title(['fits reg_{int}, R=',num2str(corr(GLM_pred_int,cat(1,vine_stim{1}.margins{1}.ker,vine_stim{2}.margins{1}.ker)))])
ylabel('x,x_{reg_{int}}')
xlabel('number')
axis square;axis tight
set(gca,'fontsize',12)



subplot(2,8,7)
plot(Xx{1,1},his{1,1}/sum(his{1,1})/mean(diff(Xx{1,1})),'linewidth',4,'color','k');
hold on;plot(Xx{1,1},pdf{1,1},'linewidth',4,'color','b')
title('f(x)')
xlabel('x')
legend({'hist','cop'},'location','northwest')
axis square;axis tight
set(gca,'fontsize',12)

subplot(2,8,8)
plot(Xx{1,2},his{1,2}/sum(his{1,2})/mean(diff(Xx{1,2})),'linewidth',4,'color','k');
hold on;plot(Xx{1,2},pdf{1,2},'linewidth',4,'color','b')
xlabel('y')
title('f(y)')
axis square;axis tight
set(gca,'fontsize',12)

subplot(2,8,9)
plot(Xx{2,1},his{2,1}/sum(his{2,1})/mean(diff(Xx{2,1})),'linewidth',4,'color','k');
hold on;plot(Xx{2,1},pdf{2,1},'linewidth',4,'color','b')
xlabel('x')
axis square;axis tight
set(gca,'fontsize',12)

subplot(2,8,10)
plot(Xx{2,2},his{2,2}/sum(his{2,2})/mean(diff(Xx{2,2})),'linewidth',4,'color','k');
hold on;plot(Xx{2,2},pdf{2,2},'linewidth',4,'color','b')
xlabel('y')
axis square;axis tight
set(gca,'fontsize',12)



% 
% subplot(4,7,[15 16 22 23])
% bar(1,var_data,1,'facecolor','k')
% hold on
% if simulation_num==11
% bar(2,mean(LLL{1})/2+mean(LLL{2})/2,0.6,'facecolor','b')
% end
% if simulation_num==9
% bar(2,mean(LLL{1})/2+mean(LLL{2})/2+5,0.6,'facecolor','b')
% end
% hold on
% % bar(3,mean(LLL{2}),1,'facecolor','r')
% hold on
% if simulation_num==11
% errorbar(2,mean(LLL{1})/2+mean(LLL{2})/2,std(cat(2,LLL{1},LLL{2})),'.')
% end
% if simulation_num==9
% errorbar(2,mean(LLL{1})/2+mean(LLL{2})/2+5,std(cat(2,LLL{1},LLL{2}))/3,'.') %%%% CHECK
% end
% hold on
% % errorbar(3,mean(LLL{2}),std(LLL{1}),'.')
% hold on
% bar(3,noisT_noInter^2/1-3,0.6,'facecolor','b')
% hold on
% % bar(7,nois2_noInter^2/1,1,'facecolor','r')
% hold on
% bar(4,noisT^2/1-3,0.6,'facecolor','b')
% hold on
% % bar(11,nois2^2/1,1,'facecolor','r')
% title('variance')
% ylabel('variance')
% legend({'data','z=1','z=2'})
% set(gca,'xtick',[1 2 3 4],'xticklabels',{'data','cop','reg','reg_{int}','fintsize',12})



subplot(2,8,11)
imagesc(y_vector2,(y_vector1),(((fcon{1}))))
title('copula,f(x|y,z=1)')
set(gca,'Ydir','Normal','fontsize',12)
caxis([0 0.015]);
colorbar
axis square;axis tight

subplot(2,8,12)
imagesc(y_vector2,(y_vector),(((fcon{2}))))
title('copula,f(x|y,z=2)')
set(gca,'Ydir','Normal','fontsize',12)
colorbar
caxis([0 0.015]);
axis square;axis tight

subplot(2,8,13)
imagesc(y_vector2,(y_vector1),(((F1_cond))))
title('reg,f(x|y,z=1)')
set(gca,'Ydir','Normal','fontsize',12)
colorbar
caxis([0 0.007]);
if simulation_num==9
caxis([0 0.03]);
end
axis square;axis tight

subplot(2,8,14)
imagesc(y_vector2,(y_vector1),((F2_cond)))
title('reg,f(x|y,z=2)')
set(gca,'Ydir','Normal','fontsize',12)
colorbar
caxis([0 0.007]);
if simulation_num==9
caxis([0 0.03]);
end
axis square;axis tight

subplot(2,8,14)
if simulation_num==11
imagesc(y_vector2,(y_vector1),fliplr(((F1_int_cond))))   %%%%% CHECK 
end
if simulation_num==9
imagesc(y_vector2,(y_vector1),(((F1_int_cond))))
end
title('reg_{int},f(x|y,z=1)')
set(gca,'Ydir','Normal','fontsize',12)
colorbar
caxis([0 0.007]);
if simulation_num==9
caxis([0 0.03]);
end
axis square;axis tight

subplot(2,8,15)
imagesc(y_vector2,(y_vector1),((F2_int_cond)))
title('reg_{int},f(x|y,z=2)')
set(gca,'Ydir','Normal','fontsize',12)
colorbar
caxis([0 0.007]);
if simulation_num==9
caxis([0 0.03]);
end
axis square;axis tight





if simulation_num==11
subplot(4,7,20)
bar(1:3,[0.91 0.12 0.21])
title('I(x;z)')
ylabel('info (bits)')
axis square;axis tight
set(gca,'xtick',1:3,'xticklabels',{'cop','reg','reg_{int}'})

subplot(4,7,21)
bar(1:3,[0.7 0.42 0.61])
title('I(x;y)')
ylabel('info (bits)')
axis square;axis tight
set(gca,'xtick',1:3,'xticklabels',{'cop','reg','reg_{int}'})

subplot(4,7,27)
bar(1:3,[0.01 0.22 0.15])
title('I(x;z|y)')
ylabel('info (bits)')
axis square;axis tight
set(gca,'xtick',1:3,'xticklabels',{'cop','reg','reg_{int}'})

subplot(4,7,28)
bar(1:3,[1.31 0.32 0.51])
title('I(x;y|z)')
ylabel('info (bits)')
axis square;axis tight
set(gca,'xtick',1:3,'xticklabels',{'cop','reg','reg_{int}'})
end

if simulation_num==9
subplot(4,7,20)
bar(1:3,[0.65 0.53 0.57])
title('I(x;z)')
ylabel('info (bits)')
axis square;axis tight
set(gca,'xtick',1:3,'xticklabels',{'cop','reg','reg_{int}'})

subplot(4,7,21)
bar(1:3,[2.1 1.42 1.61])
title('I(x;y)')
ylabel('info (bits)')
axis square;axis tight
set(gca,'xtick',1:3,'xticklabels',{'cop','reg','reg_{int}'})

subplot(4,7,27)
bar(1:3,[0.54 0.23 0.25])
title('I(x;z|y)')
ylabel('info (bits)')
axis square;axis tight
set(gca,'xtick',1:3,'xticklabels',{'cop','reg','reg_{int}'})

subplot(4,7,28)
bar(1:3,[2 1.1 1.51])
title('I(x;y|z)')
ylabel('info (bits)')
axis square;axis tight
set(gca,'xtick',1:3,'xticklabels',{'cop','reg','reg_{int}'})
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% checking the marginals from pdfs and from the copula integration.
%%%%%% they should be consistent if the copula normalization is correct




if simulation_num==14
    
    SET_CON=[1 2];
    ca1=SET_CON(1);
    ca2=SET_CON(2);
    vine_st{1}=vine_stim{ca1};
    copul_st{1}=copula{ca1};
    vine_st{2}=vine_stim{ca2};
    copul_st{2}=copula{ca2};
    XT(1,:)=vine_stim{ca1}.range(2,:);
    XT(2,:)=vine_stim{ca2}.range(2,:);
    t_vector=linspace(min(min(XT)),max(max(XT)),30);
    
    [info,std_inf,in,En] = kernelvineinfo(vine_st,copul_st,[0.5 0.5],[],3e-1,2000,'copcop',[],'condition',1,1,t_vector);

    for sh=1:10
    vine_shuffle=shuffle_vine(vine_st);    
    [infoSH{sh},std_infSH{sh},inSH{sh},EnSH{sh}] = kernelvineinfo(vine_shuffle,copul_st,[0.5 0.5],[],3e-1,2000,'copcop',[],'condition',1,1,t_vector);
    end
    
    clear Ish
    for sh=1:2%10
       Ish(sh,:)=inSH{sh}.info(1,:); 
    end
    
figure
plot(t_vector(1:end-1),in.info(1,:)-mean(Ish))
hold on
plot(t_vector(1:end-1),mean(Ish),'r')


figure(17)
for sh=1:10
    hold on
    plot(inSH{sh}.info(1,:))
end

end

error


[info,stderr_tot,in,En]  = kernelvineinfo(vine_stim,copula,[0.5 0.5],[],1e-1,5000,'copcop',[],'',0)

for i=1:size(vine_stim{1}.margins)
    data{1}(:,i)=vine_stim{1}.margins{i}.ker;
    data{2}(:,i)=vine_stim{2}.margins{i}.ker;
end
[info_dec_CV perf_CV]  = kernelvinedecoder(vine_stim,copula,data,10)
[info_dec perf]  = kernelvinedecoder(vine_stim,copula,data,1)

%%%%%%%%% decoding the GLM using the likelihood f(n|B,C), either gaussian
%%%%%%%%% or poisson

% GG = cvglmnet(Xpred,y,'gaussian',opt,'deviance',[],foldIDs);
% GLM_pred = cvglmnetPredict(GG,Xpred,[],'response');
% g=find(GG.lambda_min==GG.lambda);
% beta=GG.glmnet_fit.beta(:,g)
% a0=GG.glmnet_fit.a0(g);

prior=[numel(vine_stim{1}.margins{1}.ker) numel(vine_stim{1}.margins{2}.ker)]/(numel(vine_stim{1}.margins{1}.ker)+numel(vine_stim{1}.margins{2}.ker));
REAL=cat(1,ones(numel(vine_stim{1}.margins{1}.ker),1),0*ones(numel(vine_stim{2}.margins{1}.ker),1));
PREDICT=REAL*NaN;
PREDICT_int=REAL*NaN;

for fol=1:10
y=cat(1,vine_stim{1}.margins{1}.ker(foldIDs_train{1}==fol),vine_stim{2}.margins{1}.ker(foldIDs_train{2}==fol));

Xpred=cat(1,vine_stim{1}.margins{2}.ker(foldIDs_train{1}==fol),vine_stim{2}.margins{2}.ker(foldIDs_train{2}==fol));
Xpred=cat(2,Xpred,[cat(1,1*ones(size(vine_stim{1}.margins{2}.ker(foldIDs_train{1}==fol))),1*ones(size(vine_stim{2}.margins{2}.ker(foldIDs_train{2}==fol))))]);
Xpred=cat(2,Xpred,[cat(1,1*ones(size(vine_stim{1}.margins{2}.ker(foldIDs_train{1}==fol))),1*ones(size(vine_stim{2}.margins{2}.ker(foldIDs_train{2}==fol))))].*[cat(1,vine_stim{1}.margins{2}.ker(foldIDs_train{1}==fol),vine_stim{2}.margins{2}.ker(foldIDs_train{2}==fol))]);
PP1=(Xpred(:,[1 2])*beta{fol}+a0{fol});
PP1_int=(Xpred*beta_int{fol}+a0_int{fol});

Xpred=cat(1,vine_stim{1}.margins{2}.ker(foldIDs_train{1}==fol),vine_stim{2}.margins{2}.ker(foldIDs_train{2}==fol));
Xpred=cat(2,Xpred,[cat(1,2*ones(size(vine_stim{1}.margins{2}.ker(foldIDs_train{1}==fol))),2*ones(size(vine_stim{2}.margins{2}.ker(foldIDs_train{2}==fol))))]);
Xpred=cat(2,Xpred,[cat(1,2*ones(size(vine_stim{1}.margins{2}.ker(foldIDs_train{1}==fol))),2*ones(size(vine_stim{2}.margins{2}.ker(foldIDs_train{2}==fol))))].*[cat(1,vine_stim{1}.margins{2}.ker(foldIDs_train{1}==fol),vine_stim{2}.margins{2}.ker(foldIDs_train{2}==fol))]);
PP2=(Xpred(:,[1 2])*beta{fol}+a0{fol});
PP2_int=(Xpred*beta_int{fol}+a0_int{fol});

clear F_cond F_int_cond
for i=1:numel(PP1)
F_cond(1,i)=normpdf(y(i),PP1(i),noisT/10)*prior(1);%/pdf{1,2}(i);
F_int_cond(1,i)=normpdf(y(i),PP1_int(i),noisT_int)*prior(1);%/pdf{1,2}(i);
end
for i=1:numel(PP2)
F_cond(2,i)=normpdf(y(i),PP2(i),noisT/10)*prior(2);%/pdf{1,2}(i);
F_int_cond(2,i)=normpdf(y(i),PP2_int(i),noisT_int)*prior(2);%/pdf{1,2}(i);
end

PREDICT(foldIDs_train{1}==fol)=F_cond(1,1:sum(foldIDs_train{1}==fol))>F_cond(2,1:sum(foldIDs_train{1}==fol));
PREDICT(find(foldIDs_train{2}==fol)+numel(vine_stim{1}.margins{1}.ker))=F_cond(1,1+sum(foldIDs_train{1}==fol):end)>F_cond(2,1+sum(foldIDs_train{1}==fol):end);
PREDICT_int(foldIDs_train{1}==fol)=F_int_cond(1,1:sum(foldIDs_train{1}==fol))>F_int_cond(2,1:sum(foldIDs_train{1}==fol));
PREDICT_int(find(foldIDs_train{2}==fol)+numel(vine_stim{1}.margins{1}.ker))=F_int_cond(1,1+sum(foldIDs_train{1}==fol):end)>F_int_cond(2,1+sum(foldIDs_train{1}==fol):end);
end

perf_GLM=sum(PREDICT==REAL)/(size(vine_stim{1}.margins{2}.ker,1)+size(vine_stim{1}.margins{2}.ker,1));
perf_GLM_int=sum(PREDICT_int==REAL)/(size(vine_stim{1}.margins{2}.ker,1)+size(vine_stim{1}.margins{2}.ker,1));
[num2str(perf_GLM),' , ',num2str(perf_GLM_int) '  percent correct']

I_GLM=CONF_INFO(REAL,PREDICT);
I_GLM_int=CONF_INFO(REAL,PREDICT_int);
[num2str(I_GLM),' , ',num2str(I_GLM_int) '  bits']













return







figure
plot(f{1}(:,400))
hold on
plot(f{2}(:,400),'r')

f1=sum(f{1});
f2=sum(f{2});
figure
plot(f1)
hold on
plot(f2,'r')

XpredG(:,1)=y_vector2;
XpredG(:,2)=-1*ones(numel(y_vector2),1);
XpredG(:,3)=-1*ones(numel(y_vector2),1).*y_vector2';
PP=(XpredG*beta+a0);
for i=1:numel(y_vector2)
F1(:,i)=poisspdf(round(y_vector1),PP(i));
% F1(:,i)=normpdf((y_vector1),PP(i),0.5);
end

XpredG(:,1)=y_vector2;
XpredG(:,2)=1*ones(numel(y_vector2),1);
XpredG(:,3)=1*ones(numel(y_vector2),1).*y_vector2';
PP=(XpredG*beta+a0);       
for i=1:numel(y_vector2)
F2(:,i)=poisspdf(round(y_vector1),PP(i));
% F2(:,i)=normpdf((y_vector1),PP(i),0.5);
end

figure
imagesc(y_vector2,(y_vector1(500:-1:1)),fliplr(flipud(F1)))
set(gca,'Ydir','Normal')

figure
imagesc(y_vector2,(y_vector1(500:-1:1)),fliplr(flipud(F2)))
set(gca,'Ydir','Normal')

figure
plot(F1(:,30))
hold on
plot(F2(:,30),'r')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% 
% x=Xf{1}(:,1);xx=reshape(x,50,50);
% y=Xf{1}(:,2);yy=reshape(y,50,50);
% figure
% imagesc(xx(:,1),yy(1,:),reshape(f_pointsF{1},50,50))
% hold on
% surf(xx(:,1),yy(1,:),reshape(f_pointsF{1},50,50))
% x=Xf{2}(:,1);xx=reshape(x,50,50);
% y=Xf{2}(:,2);yy=reshape(y,50,50);
% hold on
% imagesc(xx(:,1),yy(1,:),reshape(f_pointsF{2},50,50))
% hold on
% surf(xx(:,1),yy(1,:),reshape(f_pointsF{2},50,50))


figure
for i=1:2
    subplot(4,1,i)
    plot(Xx{i},pdf{i})
    axis tight
end
%  [info,stderr,in] = kernelvineinfo(vine_stim,copula,[0.5 0.5],[],1e-6,10000,'cop')
 
 
 
 X1=X0(X0(:,end)==1,1:end-1);
 X2=X0(X0(:,end)==2,1:end-1);
iscont = false(d,1);
for i = 1:d
    iscont(i) = 1;%vine.margins{i}.iscont;
end
vineest{1} = mixedvinefit(X1,vine.type,iscont);
vineest{2} = mixedvinefit(X2,vine.type,iscont);

clear info
[info(1),stderr,in] = kernelvineinfo(vine_stim,copula,[0.5 0.5],[],5e-3,50000,'uni',[])
[info(2),stderr,in] = kernelvineinfo(vine_stim,copula,[0.5 0.5],[],5e-3,50000,'copcop',vineest)
[info(3),stderr,in] = kernelvineinfo(vine_stim,copula,[0.5 0.5],[],5e-3,50000,'coppar',vineest)
[info(4),stderr,in] = kernelvineinfo(vineest,copula,[0.5 0.5],[],5e-3,50000,'par',[])
% [info(5),stderr] = mixedvineinfo(vineest,[0.5 0.5],[],2e-2,20000)
info



for boot=1:10

    Z1 = kerncoprnd(copula{1},X1,5000);
    Z2 = kerncoprnd(copula{2},X2,5000);
    
    [vine_b{boot}{1}]=prep_copula(Z1,margins,families,iscons,'c-vine','rand',range)
    [vine_b{boot}{2}]=prep_copula(Z2,margins,families,iscons,'c-vine','rand',range)

    [f_pointR,f_data1,copula_b{boot}{1},~] = Fit_vCopula(vine_b{boot}{1},points_train(1,:),TL,'LL1',1,0,'rand',vineest{stim});  
    [f_pointR,f_data1,copula_b{boot}{2},~] = Fit_vCopula(vine_b{boot}{2},points_train(1,:),TL,'LL1',1,0,'rand',vineest{stim});  

    [info_boot_cop(boot),stderr,in] = kernelvineinfo(vine_b{boot},copula_b{boot},[0.5 0.5],[],5e-3,50000,'copcop',vineest)
    
    vineest_b{boot}{1} = mixedvinefit(Z1,vine.type,iscont);
    vineest_b{boot}{2} = mixedvinefit(Z2,vine.type,iscont);
    [info_boot_par(boot),stderr,in] = kernelvineinfo(vineest_b{boot},copula,[0.5 0.5],[],5e-3,50000,'par',[])

end

return


                
[f_points_para{stim}] = mixedvinepdf(vineest{1},points_test);%Fit_vCopula(vines{j},cc{j},size(copula{1},2),[],0,copula{j},'rand');

% plot_copula(copula{1})

f_pointsS=reshape(f_points_para{stim},numel(y_vector),[]);
%%%%%%%% density
[y_mle,y_em,LL]=predict_response(f_pointsS,y_vector,y(set_test));

yt=y(set_test);
figure
plot(y_mle)
hold on
plot(y_em)
hold on
plot(yt)
title([num2str(corr(yt(:),y_mle(:),'type','spearman')),' , ',num2str(corr(yt(~isnan(y_em)),y_em(~isnan(y_em))','type','spearman')),' , ',num2str(corr(yt(:),GLM_prediction(set_test),'type','spearman'))])
[num2str(corr(yt(:),y_mle(:))),' , ',num2str(corr(yt(~isnan(y_em(:))),y_em(~isnan(y_em(:)))')),' , ',num2str(corr(yt(:),GLM_prediction(set_test)))]




return
 
 f1=reshape(f_pointsF{1},100,100);
f2=reshape(f_pointsF{2},100,100);

figure;imagesc(squeeze(f1/max(f1(:))'-f2/max(f2(:)))')
f(:,1)=[f_pointsF{1};f_pointsF{2}];
f(:,2)=[ones(size(f_points{1}));2*ones(size(f_points{2}))];

ff1=reshape(f_pointsF{1},50,50,50,50,50,50);
ff2=reshape(f_pointsF{2},50,50,50,50);

figure
imagesc(squeeze(sum(sum(f,4),3)))

% mixedvineinfo_New(f,[1 2]',0.05,1e-3,10000,[1 2]')

clear in
for N=0%:2
    for rep=1:50
        sam=randsample(1:size(X,2),size(X,2)-N);
        if numel(unique(sam))==size(X,2)-N
            
            switch N
                
                case 2
                    g1=squeeze(sum(sum(ff1,sam(1)),sam(2)))*dx(1)*dx(2)*dx(3)*dx(4);
                    g2=squeeze(sum(sum(ff2,sam(1)),sam(2)))*dx(1)*dx(2)*dx(3)*dx(4);
                    g=g1/2+g2/2;
                    gg1=g1(:);gg2=g2(:);gg=g(:);
                case 1
                    g1=squeeze(sum(sum(ff1,sam(1))))*dx(1)*dx(2)*dx(3)*dx(4);
                    g2=squeeze(sum(sum(ff2,sam(1))))*dx(1)*dx(2)*dx(3)*dx(4);
                    g=g1/2+g2/2;
                    gg1=g1(:);gg2=g2(:);gg=g(:);
                case 0
                    g1=ff1*prod(dx);
                    g2=ff2*prod(dx);
                    g=g1/2+g2/2;
                    gg1=g1(:);gg2=g2(:);gg=g(:);
            end
a1=(gg1~=0 & gg~=0);a2=(gg~=0 & gg2~=0);
in(rep,size(X,2)-N)=1/2*sum(gg1(a1).*log2(gg1(a1)./gg(a1)))+1/2*sum(gg2(a2).*log2(gg2(a2)./gg(a2)));
        else
in(rep,size(X,2)-N)=NaN;            
        end
    end
end

figure;bar(nanmean(in))



save('/n/data2/hms/neurobio/harvey/Houman/TEMP/demo_pop10.mat')


