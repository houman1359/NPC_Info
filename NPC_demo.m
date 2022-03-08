
clear all
close all


addpath(genpath('/home/hs258/Codes_Folder/Houman_Git/NPC_Info/NPC_Info'))


load('group.mat')

opts.bw='LL1'; %% LL1 or LL2 bandwidth methods
opts.knots_fit=50; %% number of bins used in estimating the copula
opts.knots_est=opts.knots_fit;
opts.type={'cont','cont'}; %% continous or discrete marginals
% opts.type={'discrete','discrete'};
opts.parallel=0; %% 0=non paralle, 1=parallel computing
opts.alpha=0.05; %% alpha for monte-carlo sampling varinace
opts.erreps=1e-3; %% variance threshhold for monte-carlo sampling for information estimation
opts.iter=5; %% maximum number of monte-carlo iterations
opts.cases=20000; %% number of samples in each iteration
opts.plot=1; %% 0=no iteration plot (default), 1=with iteration plot
% Computing I(X(:,1);X(:,2)):

tic
for gr=1:3
    D=0;
X0=group{gr}';
if size(X0,1)>10000
    ind=randsample(1:size(X0,1)-10000,1);
    X0=X0(ind:ind+10000,:);
end

for DEL=0%-10:10

    [gr DEL]
    D=D+1;

X=X0;
X(:,1)=circshift(X(:,1),DEL);

range(1:2,1)=min(X(:,1:2))-1e-30;
range(1:2,2)=max(X(:,1:2))+1e-30;
[vine]=NPC_prep_copula(X(:,[1 2]),opts.type,range([1 2],:));

[ density_X , ~ , copula , ~ , ~ ] = NPC_Fit_vCopula(vine,X(1,:),opts.bw,1,0,opts.knots_fit,opts.parallel);
[ ~ , ~ , copula , ~ , ~ ] = NPC_Fit_vCopula(vine,X(1,:),opts.bw,-1,copula,opts.knots_est,opts.parallel);

[ info(gr,D) , ~ , ~ , ~ ] = NPC_kernelvineinfo(vine,copula,opts)

end
end
toc
% figure
% plot(info')

figure
bar(info)

