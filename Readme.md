The sample bivariate data:
load('NPC_X.mat')

The parameters which can be selected in the copula fitting and information estimation are as follows:

opts.bw='LL1';                               %% LL1 or LL2 bandwidth methods
opts.knots_fit=100;                      %% number of bins used in estimating the copula 
opts.knots_est=opts.knots_fit;
opts.type={'cont','cont'};              %% continous or discrete marginals
% opts.type={'discrete','discrete'};
opts.parallel=0;                             %% 0=non paralle, 1=parallel computing 
opts.alpha=0.05;                           %% alpha for monte-carlo sampling varinace
opts.erreps=1e-3;                         %% variance threshhold for monte-carlo sampling for information estimation
opts.iter=50;                                  %% maximum number of monte-carlo iterations 
opts.cases=20000;                        %% number of samples in each iteration 
opts.plot=0;                                   %% 0=no iteration plot (default), 1=with iteration plot

Computing I(X(:,1);X(:,2)):
First we build the vine structure
range(1:2,1)=min(X(:,1:2))-1e-10;
range(1:2,2)=max(X(:,1:2))+1e-10;
[vine]=NPC_prep_copula(X(:,[1 2]),opts.type,range([1 2],:));

The the bandwidths of copula can be fitted on the data:
[density_X,~,copula,~,~] = NPC_Fit_vCopula(vine,X(1,:),opts.bw,1,0,opts.knots_fit,opts.parallel);
The copula can be estimated over the grid and over the data points:
[~,~,copula,~,~] = NPC_Fit_vCopula(vine,X(1,:),opts.bw,-1,copula,opts.knots_est,opts.parallel);
Finally the mutual information can be estimated by callingL
[info,~,~,~] = NPC_kernelvineinfo(vine,copula,opts)