% Huggett, Ventura & Yaron (2006) - Human Capital and Earnings Distribution Dynamics
% This model is simple, just a deterministic life-cycle model.
% The main complication is that there are a lot of permanent types.
% Note that the permanent types are part of the calibration, and moreover the initial distribution 
% is part of the calibration and is a joint-distribution over the ability/type and the human capital levels.

% Could model human capital, h, as an experience asset. But this model is so simple we just brute force it.
n_d=101; % share of time to invest in new human capital
n_h=501; % Human capital [I call this n_h, as if I call it n_a I got confused about it being human capital not ability]
n_z=0;
N_i=21; % number of different ability levels
% Note: HVY2006, pg 277 says initial distribution was based on a 
% rectangular grid with 20 points on each of human capital and learning
% ability. Here I use 21 points for ability as I just prefer to use odd
% numbers when discretizing. [I use way more points for human capital.]

vfoptions.experienceasset=1; % Using an experience asset
simoptions.experienceasset=1;

N_j=39; % 39 or 49, depending on exercise [starting age 20, starting age 10, respectively]
if N_j==39
    Params.agejshifter=19; % starting age 20, to correspond to 39 periods
elseif N_j==49
    Params.agejshifter=9; % starting age 10, to correspond to 49 periods
end

figure_c=0;
CreateFigures=1 % Set to 0 to skip creating figures (so can be run headless)
DoDataWork=1 % Set to 0 to skip doing the original data work (just loads results)

% If you are running this with DoDataWork=1, you need to download the PSID
% data set that it uses:
% Download the data from PSID, then import
% Go to PSID and look for previous "carts"
% https://simba.isr.umich.edu/DC/c.aspx
% Search for and download the cart called 
% "325247 - HVY2006" associated with the email address "robertdkirkby@gmail.com"
% This will download a zip file, unzip it and put the files in the same
% folder as this current m-file.
% Rename the file which is .sps to _sps.txt
% Rename the file which is _format.sps to _formate_sps.txt
% (if you don't rename these ImportPSIDdata will just give a warning telling you to do so)


%% Parameters
Params.J=N_j;
Params.agej=1:1:Params.J; % the model period

% Human capital
Params.delta=0.0114; % depreciation rate of human capital
Params.alpha=0.7; % returns to scale of human capital production fn
% Note: HVY2006 try five different values of alpha: 0.5,0.6,0.7,0.8,0.9

% Interest rate
Params.r=0.04; % real interest rate (used to determine the discount factor)
Params.beta=1/(1+Params.r); % the discount factor (HVY2006 just call it 1/(1+r) )

% Wage growth
Params.g=0.0014;
Params.w=(1+Params.g).^(Params.agej-1); % wage rate (per unit of human-capital-hour)

%% Parameters for ability and for initial dist (on ability and human capital)

Params.mean_h1=92.3;
Params.coeffofvariation_h1=0.481;
Params.mean_ability=0.209;
Params.coeffofvariation_ability=0.347;
Params.corr_h1ability=0.781;
% These parameters are taken from Table 5; they are for parametric case with initial age 20 and alpha=0.7

% From these we can get the standard deviations
Params.stddev_h1=Params.coeffofvariation_h1*Params.mean_h1; % coeff. of variation=std. dev./mean, is the definition of the coeff. of variation
Params.stddev_ability=Params.coeffofvariation_ability*Params.mean_ability;

% But everything is really done around log(ability) and log(human capital)
% So we actually need the mean, std dev and correlations for the logs.
% There is no way to do this for std dev and correlation, so just have to
% stick in some guesses instead (we are anyway going to calibrate all these
% parameters later so no big deal).

% We want to discretize the log-normal distribution, currently we have mean
% and std. dev. and corr for the variables, we need to figure out the same
% for the logs (as that is what we will discretize as normal).
% Cannot just do this from formulas (maybe you can, you could with a single
% log-normal, the bivariate complicates things, I am not aware of formula,
% maybe there is one but I am fairly confident there is not).
Params.mean_logh1=log(Params.mean_h1);
Params.stddev_logh1=0.1; % PLACEHOLDER
Params.mean_logability=log(Params.mean_ability);
Params.stddev_logability=0.1; % PLACEHOLDER
Params.corr_logh1logability=0.5; % PLACEHOLDER

%%
% maxh is going to be coming from delta*h=ability*(h*l)^alpha
% Rearrange gives h=(ability/delta)^(1/(1-alpha))
% max of ability is roughly 0.4 (take a look at ability_grid below)
hmax=(0.4/Params.delta)^(1/(1-(Params.alpha)));
% Doing this gives hmax of around 10000, but solving model a few times is
% clear noone ever goes anywhere near this so I just use hmax=1000
h_grid=linspace(0.1,1000,n_h)'; % grid on human capital

l_grid=linspace(0,1,n_d)'; % grid on human capital investment time

% rename in toolkit notation
n_a=n_h;
d_grid=l_grid;
a_grid=h_grid;
z_grid=[];
pi_z=[];

%%  Set up the permanent types
% Ability is correlated with initial human capital, but we just ignore this correlation while we set up the grid on ability
% (we will consider the correlation later when putting agents onto the grids)

[ability_grid,pi_ability]=discretizeAR1_FarmerToda(Params.mean_logability,0,Params.stddev_logability,N_i);
ability_grid=exp(ability_grid);
pi_ability=pi_ability(1,:)'; % iid
% Comment: originally I just created a grid on ability, rather than log of ability, but this gave some points a negative 
% value for ability, so switched to current setup of creating the grid on log of ability, then take exponent.
% Note that the initial distribution is anyway log-normal, so working with a grid on log-ability makes more sense for that anyway.

% Ability is a parameter that depends on permanent type
Params.ability=ability_grid;

PTypeDistParamNames={'abilitydist'};
Params.abilitydist=pi_ability; % Note: This is not really used, as is replaced with the masses from the initial distribution which accounts for the correlation between initial human capital and initial assets


%% Set up the experience asset, we need to aprimeFn which gives value of aprime given d2 and a2 (d2 is the decision variable relevant to experience asset, a2 is the experience asset)
% f_hl=ability*(h*l)^alpha; % new production of human capital
% hprime=f_hl+h*(1-delta);

% aprimeFn: hprime(l,h,parameters)
vfoptions.aprimeFn=@(l,h,ability,alpha,delta) ability*(h*l)^alpha+h*(1-delta);
simoptions.aprimeFn=vfoptions.aprimeFn;
simoptions.a_grid=a_grid; % Needed when using experience asset
simoptions.d_grid=d_grid; % Needed when using experience asset

%% Following shows how much h can increase before maxing out. Looking at
% this there is no max h (there is, but it is large), but because we have finite periods there is going
% to be a maximum that really just comes from the N_j
if CreateFigures==1
    figure_c=figure_c+1;
    figure(figure_c)
    subplot(2,1,1); plot(h_grid,Params.ability(1)*(h_grid.*1').^Params.alpha+h_grid.*(1-Params.delta), h_grid, h_grid)
    subplot(2,1,2); plot(h_grid,Params.ability(end)*(h_grid.*1').^Params.alpha+h_grid.*(1-Params.delta), h_grid, h_grid)
    legend('human capital prodn','45 degree')
end

%% 
DiscountFactorParamNames={'beta'};

ReturnFn=@(l,h,w)...
    HuggettVenturaYaron2006_ReturnFn(l,h,w);

%% Now solve the value function iteration problem, just to check that things are working before we go to General Equilbrium
tic;
[V, Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,N_i, d_grid, a_grid, ability_grid, pi_ability, ReturnFn, Params, DiscountFactorParamNames, vfoptions);
vftime=toc


%% HVY2006 consider both a parametric and a non-parametric initial distribution over initial human  capital, h1 and ability.
% Here we will just do the parametric, which uses a bi-variate log-normal distribution.

Params.logh1logability_Mean=[Params.mean_logh1,Params.mean_logability];

% Create covariance matrix for std devs and correlation
CorrMatrix=[1, Params.corr_logh1logability; Params.corr_logh1logability, 1];
Params.logh1logability_CoVarMatrix = corr2cov([Params.stddev_logh1,Params.stddev_logability],CorrMatrix);

% Now we can find the probabilities of a bivariate log-normal distribution over the existing grids
tic;
P=MVNormal_ProbabilitiesOnGrid([log(h_grid);log(ability_grid)],Params.logh1logability_Mean, Params.logh1logability_CoVarMatrix, [n_h,N_i]);
% Note: P is on defined on (h,ability).
mvntime=toc

tic;
P=gpuArray(MVNormal_ProbabilitiesOnGrid(gather([log(h_grid);log(ability_grid)]),gather(Params.logh1logability_Mean), gather(Params.logh1logability_CoVarMatrix), [n_h,N_i]));
% Note: P is on defined on (h,ability).
mvntime2=toc % This is actually way way faster for some reason

% Make sure all the agent ptypes have positive (non-zero) probabilities [as otherwise this would cause problems]
sum(P,1) % with N_i=21 the highest and lowest are roughly 10 to minus eight
% sum(P,1)==0 shows they are all non-zero (only just!)

% Initial distribution of agents at birth (j=1)
% Before we plot the life-cycle profiles we have to define how agents are
% at age j=1. We will give them all zero assets.
% jequaloneDist=zeros(n_a,N_i,'gpuArray'); % Put no households anywhere on grid (I just leave this here to show what size it is)
jequaloneDist=P; % joint log-normal distribution onto our existing grids


%% Now replace this with initial agent distribution as a function as we need that so that the parameters of this can be part of our calibration
% jequaloneDist as a function has first four inputs (a_grid,z_grid,n_a,n_z,...) followed by any parameters
jequaloneDistFn=@(a_grid,z_grid,n_a,n_z,mean_logh1,mean_logability,stddev_logh1,stddev_logability,FTcorr_logh1logability, ability, Number_i) HuggettVenturaYaron2006_jequaloneDistFn(a_grid,z_grid,n_a,n_z,mean_logh1,mean_logability,stddev_logh1,stddev_logability,FTcorr_logh1logability, ability, Number_i);
% To use jequaloneDist as a function you must include a_grid and z_grid in simoptions
simoptions.a_grid=h_grid;
simoptions.z_grid=[]; % Even though we don't use z_grid, it is still needed in simoptions to be able to call jequaloneDist as a function
% Params.ability=ability_grid; % so we can pass it as an input to jequaloneDistFn, We already set this above
Params.Number_i=N_i; % so we can pass it as an input to jequaloneDistFn

% Our parametrization of the covariance matrix follows AH2021:
% So we need to set up initial guess for
Params.FTcorr_logh1logability=0.5*log((1+Params.corr_logh1logability)/(1-Params.corr_logh1logability)); % This formula is exact when using AH2021 for 2x2 Covar matrix


%% We now compute the 'stationary distribution' of households
% Start with a mass of one at initial age, use the conditional survival
% probabilities sj to calculate the mass of those who survive to next
% period, repeat. Once done for all ages, normalize to one
Params.mewj=ones(1,Params.J); % Marginal distribution of households over age
Params.mewj=Params.mewj./sum(Params.mewj); % Normalize to one
AgeWeightParamNames={'mewj'}; % So VFI Toolkit knows which parameter is the mass of agents of each age
tic;
StationaryDist=StationaryDist_Case1_FHorz_PType(jequaloneDistFn,AgeWeightParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_z,N_j,N_i,[],Params,simoptions);
statdisttime=toc

%% Calculate the life-cycle profiles
FnsToEvaluate.humancapital=@(l,h) h; % h is human capital
FnsToEvaluate.timel=@(l,h) l; % l is fraction of time invested in human capital

tic;
AgeConditionalStats=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,[],simoptions);
acstime=toc

if CreateFigures==1
    figure_c=figure_c+1;
    figure(figure_c)
    subplot(2,1,1); plot(Params.agejshifter+Params.agej, AgeConditionalStats.humancapital.Mean)
    title('Human capital: mean life-cycle profile')
    subplot(2,1,2); plot(Params.agejshifter+Params.agej, AgeConditionalStats.timel.Mean)
    title('Fraction of time learning: mean life-cycle profile')
end
% Note, we can check our distribution of the initial distribuiton by
% looking at the means and correlation for h and ability

%% Let's do some data work to prepare for the model calibration
% Following code imports PSID data, then estimates the life-cycle profiles
% for earnings that we use for model calibration
if DoDataWork==1
    HVY2006_DataWork
    save ./SavedOutput/HVY2006_Data.mat Fig1_Mean Fig1_Median Fig1_Gini
else
    load ./SavedOutput/HVY2006_Data.mat 
end

%% Calibration exercise of HVY2006
% Choose the parameters that parametrize the initial distribution over
% human capital and ability to minimize the sum of squares of the log
% ratios for the life-cycle profile of earnings in mean, Gini, and
% ratio-of-mean-to-median.
% Here we will instead target mean, Gini and median. (Because of how VFI
% Toolkit works, median is an option, but ratio-of-mean-to-median is not.
% Given we target mean and median instead of mean and ratio-of-mean-to-median
% it seems likely that any difference should be very minimal.)

%% One more thing, we need to deal with the 'ability'
% ability is a parameter that depends on PType
% But it is also determined/parametrized by mean_logability and stddev_logability which
% are parameters to be calibrated.
% This idea of a PType parameter having a parametrized distribution across PTypes
% is common, and the way VFI Toolkit handles it is with PTypeParamFn

PTypeParamFn=@(Params) HVY2006_AbilityGrid(Params);
% PTypeParamFn will always have Params as it's first and only entry
% Because any parameter that depends on PType will have to be kept in
% Params, the output of the function will also just be Params.

%% Now, we calibrate the model parameters to target the dlife cycle profile of earnings in mean, median, and Gini

% To do this using the toolkit there are two things we need to setup
% First, just say all the model parameters we want to calibrate
% For exercise of HVY2006 these are all the parameters that define the
% initial distribution over human capital and ability.
CalibParamNames={'mean_logh1','mean_logability','stddev_logh1','stddev_logability','FTcorr_logh1logability'};
% CalibParamNames gives the names of all the model parameters that will be
% calibrated, these parameters must all appear in Params, and the values in
% Params will be used as the initial values.
% All other parameters in Params will remain fixed.

% Because the third and fourth parameters are standard deviations we want to constrain them to be positive
caliboptions.constrainpositive=[0,0,1,1,0];
% Because the fifth parameter is a correlation we want to constrain it to be 0 to 1
caliboptions.constrain0to1=[0,0,0,0,1];

% Second, we need to say which model statistics we want to target
% We can target any model stat generated by the AllStats, and LifeCycleProfiles commands
% We set up a structure containing all of our targets
TargetMoments.AgeConditionalStats.earnings.Mean=Fig1_Mean';
TargetMoments.AgeConditionalStats.earnings.RatioMeanToMedian=Fig1_Mean'./Fig1_Median'; % Note: Median, whereas HVY2006 target ratio-of-mean-to-median
TargetMoments.AgeConditionalStats.earnings.Gini=Fig1_Gini';
% Note: When setting up TargetMoments there are some rules you must follow
% There are two options TargetMomements.AgeConditionalStats and
% TargetMoments.AllStats (you can use both). Within these you must 
% follow the structure that you get when you run the commands
% AgeConditionalStats=LifeCycleProfiles_FHorz_Case1()
% and
% AllStats=EvalFnOnAgentDist_AggVars_FHorz_Case1()
% [If you are using PType, then these will be the PType equivalents]

% So we want a FnsToEvaluate which is just earnings (this will be faster as it only includes what we actually need)
clear FnsToEvaluate
FnsToEvaluate.earnings=@(l,h,w) w*h*(1-l);

% We want to set a few options, one relates to the objective function
caliboptions.metric='sum_logratiosquared';


% Prior to calibration, set the initial guess for the parameters
% Params.mean_logh1=4;
% Params.mean_logability=-2;
% Params.stddev_logability=0.5;
% Params.stddev_logh1=0.3;
% Params.FTcorr_logh1logability=0.5;
% I am going to overwrite with a better guess based on a few iterations of the calibration commmand
Params.mean_logh1=4.46;
Params.mean_logability=-1.5;
Params.stddev_logability=0.34;
Params.stddev_logh1=0.32;
Params.FTcorr_logh1logability=0.98;

save ./SavedOutput/HVY2006precalib.mat

% load ./SavedOutput/HVY2006precalib.mat


%% Done. Now just run the calibrate life-cycle model command
[CalibParams1,calibsummary1]=CalibrateLifeCycleModel_PType(CalibParamNames,TargetMoments,n_d,n_a,n_z,N_j,N_i,d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, jequaloneDistFn,AgeWeightParamNames, PTypeDistParamNames, PTypeParamFn, FnsToEvaluate, caliboptions, vfoptions,simoptions);
% CalibParams is the calibrated parameter values
% calibsummary is a structure containing various info on how the calibration went

save ./SavedOutput/HVY2006calib.mat CalibParams1 Params calibsummary1

% load ./SavedOutput/HVY2006calib.mat

% From now on, use our calibrated parameter values
for pp=1:length(CalibParamNames)
    Params.(CalibParamNames{pp})=CalibParams1.(CalibParamNames{pp});
end

%% Report some of the things in HVY2006
% We already replicated Figure 1 (in HVY2006_DataWork)
% We already replicated Figure 2 (in HVY2006_DataWork)
% We already replicated most of Figure 3 (in HVY2006_DataWork; we did time fixed effects and cohort fixed effects, we skipped restricted time fixed effects)

% Table 1 is just the parameters

% Not needed for HVY2006 plots, but I want to see these
FnsToEvaluate.humancapital=@(l,h) h; % h is human capital
FnsToEvaluate.timel=@(l,h) l; % l is fraction of time invested in human capital


% Figure 4 is same as 5, but for non-parametric calibration instead of parametrized calibration.
% Figure 5
Params=HVY2006_AbilityGrid(Params);
[V, Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,N_i,d_grid, a_grid, [], [], ReturnFn, Params, DiscountFactorParamNames, vfoptions);
jequaloneDist=jequaloneDistFn(a_grid,z_grid,n_a,n_z,Params.mean_logh1,Params.mean_logability,Params.stddev_logh1,Params.stddev_logability,Params.FTcorr_logh1logability, Params.ability, Params.Number_i);
StationaryDist=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_z,N_j,N_i,[],Params,simoptions);
AgeConditionalStats=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,[],simoptions);

% Not needed for anything, but just going to also create
AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1_PType(StationaryDist, Policy, FnsToEvaluate, Params,n_d,n_a,n_z,N_j,N_i,d_grid, a_grid, z_grid, simoptions);
AllStats=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist, Policy, FnsToEvaluate, Params,n_d,n_a,n_z,N_j,N_i,d_grid, a_grid, z_grid, simoptions);


if CreateFigures==1
    % Draw a Figure 5 that compares the model and target moments (add median alongside the other three)
    figure(5)
    subplot(4,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.Mean,Params.agejshifter+(1:1:N_j),TargetMoments.AgeConditionalStats.earnings.Mean)
    title('Mean Earnings')
    legend('Model','Data')
    subplot(4,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.Gini,Params.agejshifter+(1:1:N_j),TargetMoments.AgeConditionalStats.earnings.Gini)
    title('Gini Earnings')
    legend('Model','Data')
    subplot(4,1,3); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.RatioMeanToMedian,Params.agejshifter+(1:1:N_j),TargetMoments.AgeConditionalStats.earnings.RatioMeanToMedian)
    title('Skewness (Mean/Median) Earnings')
    legend('Model','Data')
    subplot(4,1,4); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.Median,Params.agejshifter+(1:1:N_j),Fig1_Median)
    title('Median Earnings')
    legend('Model','Data')


    % Do a plot of the initial dist to be able to see it
    figure(4)
    subplot(2,1,1); plot(Params.ability,cumsum(Params.abilitydist),Params.ability,cumsum(sum(jequaloneDist,1)))
    title('ability dist')
    legend('from Params', 'from jequaloneDist')
    subplot(2,1,2); plot(a_grid,cumsum(sum(jequaloneDist,2)))
    title('asset dist')

    % Plot mean earnings conditional on ptype
    figure(3)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype001.Mean)
    hold on
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype002.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype003.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype004.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype005.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype006.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype007.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype008.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype009.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype010.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype011.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype012.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype013.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype014.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype015.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype016.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype017.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype018.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype019.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype020.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype021.Mean)
    hold off
    title('Mean earnings conditional on ptype')
    % Plot median earnings conditional on ptype
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype001.Median)
    hold on
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype002.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype003.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype004.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype005.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype006.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype007.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype008.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype009.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype010.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype011.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype012.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype013.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype014.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype015.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype016.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype017.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype018.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype019.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype020.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.earnings.ptype021.Median)
    hold off
    title('Median earnings conditional on ptype')


    % Plot mean timel conditional on ptype
    figure(3)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype001.Mean)
    hold on
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype002.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype003.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype004.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype005.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype006.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype007.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype008.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype009.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype010.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype011.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype012.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype013.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype014.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype015.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype016.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype017.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype018.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype019.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype020.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype021.Mean)
    hold off
    title('Mean timel (time in human capital creation) conditional on ptype')
    % Plot median timel conditional on ptype
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype001.Median)
    hold on
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype002.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype003.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype004.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype005.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype006.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype007.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype008.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype009.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype010.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype011.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype012.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype013.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype014.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype015.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype016.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype017.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype018.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype019.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype020.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.timel.ptype021.Median)
    hold off
    title('Median timel (time in human capital creation) conditional on ptype')


    % Plot mean humancapital conditional on ptype
    figure(1)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype001.Mean)
    hold on
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype002.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype003.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype004.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype005.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype006.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype007.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype008.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype009.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype010.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype011.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype012.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype013.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype014.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype015.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype016.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype017.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype018.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype019.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype020.Mean)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype021.Mean)
    hold off
    title('Mean human capital conditional on ptype')
    % Plot median humancapital conditional on ptype
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype001.Median)
    hold on
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype002.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype003.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype004.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype005.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype006.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype007.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype008.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype009.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype010.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype011.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype012.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype013.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype014.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype015.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype016.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype017.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype018.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype019.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype020.Median)
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.humancapital.ptype021.Median)
    hold off
    title('Median human capital conditional on ptype')

    % So issue is coming from the human capital distribution
    figure(1);
    plot(h_grid,cumsum(StationaryDist.ptype001(:,1)))
    hold on
    for jj=2:N_j
        plot(h_grid,cumsum(StationaryDist.ptype001(:,jj)))
    end
    hold off
end


% Double-check: following is zero
cumsum(StationaryDist.ptweights)-cumsum(sum(jequaloneDist,1))'


% Table 2
fprintf('Table 2 of HVY2006 (alpha=%1.1f, parametric, accumulation starts at age %i) \n', Params.alpha, Params.agejshifter+1)
% This just reports the 'mean absolute deviation'
AbsOfMean=abs(log(TargetMoments.AgeConditionalStats.earnings.Mean./AgeConditionalStats.earnings.Mean));
AbsOfMedian=abs(log(TargetMoments.AgeConditionalStats.earnings.RatioMeanToMedian./AgeConditionalStats.earnings.RatioMeanToMedian));
AbsOfGini=abs(log(TargetMoments.AgeConditionalStats.earnings.Gini./AgeConditionalStats.earnings.Gini));
MeanAbsDeviation=(sum(AbsOfMean)+sum(AbsOfMedian)+sum(AbsOfGini))/(3*N_j);
fprintf('Mean absolute log deviation is %8.2f (%%) \n',100*MeanAbsDeviation)
% Note that mean absolute deviation is not actual measure used in the calibration exercise, that is the square, not absolute, so let's calculate that too
SqLogOfMean=log(TargetMoments.AgeConditionalStats.earnings.Mean./AgeConditionalStats.earnings.Mean).^2;
SqLogOfRatioMeanToMedian=log(TargetMoments.AgeConditionalStats.earnings.RatioMeanToMedian./AgeConditionalStats.earnings.RatioMeanToMedian).^2;
SqLogOfGini=log(TargetMoments.AgeConditionalStats.earnings.Gini./AgeConditionalStats.earnings.Gini).^2;
SumSqLogDeviation=(sum(SqLogOfMean)+sum(SqLogOfRatioMeanToMedian)+sum(SqLogOfGini))/(3*N_j);
fprintf('End of outputs relating to Table 2 \n')
fprintf('Not in Table 2, but might be interesting comparison as it is in the same metric as our calibration \n')
fprintf('Sum square log deviation is %8.2f (%%) \n',100*SumSqLogDeviation)
fprintf('This sum square log deviation is the sum across ages (rows of below) and mean-meanmedianratio-Gini (columns)')
[SqLogOfMean',SqLogOfRatioMeanToMedian',SqLogOfGini']


% Table 3

% Table 5 (Table 4 is same, but for non-parametric case)
% Note this is about ability and initial human capital; calibrated parameters were about joint log-normal distribution of this
FnsToEvaluate.ability=@(l,h,ability) ability;
FnsToEvaluate.humancapital=@(l,h) h;
AgeConditionalStats=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,[],simoptions);
fprintf('Table 5 of HVY2006: for panel A, alpha=%1.1f and accumulation starting at age=%i \n', Params.alpha, Params.agejshifter+1)
fprintf('Mean ability at initial age=%8.3f \n',AgeConditionalStats.ability.Mean(1))
fprintf('Coeff. of variation of ability at initial age=%8.3f \n',AgeConditionalStats.ability.StdDeviation(1)/AgeConditionalStats.ability.Mean(1)) % Coeff of variation=std dev/mean
fprintf('Skewness of ability at initial age=%8.3f \n',AgeConditionalStats.ability.RatioMeanToMedian(1)) % HVY2006 measure skewness as mean/median
fprintf('Mean human capital at initial age=%8.3f \n',AgeConditionalStats.humancapital.Mean(1))
fprintf('Coeff. of variation of  human capital at initial age=%8.3f \n',AgeConditionalStats.humancapital.StdDeviation(1)/AgeConditionalStats.humancapital.Mean(1)) % Coeff of variation=std dev/mean
fprintf('Skewness of  human capital at initial age=%8.3f \n',AgeConditionalStats.humancapital.Mean(1)/AgeConditionalStats.humancapital.Median(1)) % HVY2006 measure skewness as mean/median
% fprintf('Correlation of ability and human capital at initial age=%8.3f \n')
fprintf('End of outputs relating to Table 5 \n')

% Table 6
% Autocorrelation of earnings, in levels and in growth rates


save ./SavedOutput/HVY2006full.mat

% load  ./SavedOutput/HVY2006full.mat



%% Comment:
% What this code does is actually way more computationally complex than the
% exercise of HVY2006. The trick is to notice that none of the parameters
% being calibrated actually matter for the value and policy function. Hence
% HVY2006 solve the value/policy function once, and then do their
% calibration step just looping over the initial dist/stationary
% dist/life-cycle profile steps (ours includes the redundant value/policy
% fn step inside the iterations). Because the value/policy function step is 
% one of the harder parts computationally this makes a big difference. But 
% with 20 odd years of extra computing power we can anyway just
% solve the model even with this extra computation and it shows off
% how the calibration command in VFI Toolkit works.






