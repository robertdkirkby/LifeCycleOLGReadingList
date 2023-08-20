% Carroll (1997) - Buffer-Stock Saving and the Life-Cycle/Permanent Income Hypothesis
%
% This model has permanent shocks. In a model with permanent shocks what
% you should do is renormalize the value function problem so that the permanent
% shocks disappear from the state-space and just act like a stochastic discount
% factor (you can see how to code this at https://github.com/robertdkirkby/example_permamentincomeshock )
% But permanent shocks don't really help you get to the modern literature
% with age-dependent shock processes, so instead I just solve the whole
% value function problem with the permanent shock.
%
% The value function here does not use the 'gross wealth' setup. This code
% just directly solves the model in Equation 1 of C1997.
%
% Because the shock space here is large, the agent distribution takes a while to solve.

%% Basic setup

n_a=401; % Assets
n_z=51; % Permanent shock (P in C1997 notation)
n_e=12; % zero-income shock (V in C1997 notation) [working paper, Appendix I, makes it clear there were 11 points for Z, so one more for the zero gives us 12]

N_j=50; % Carroll (1997) makes no mention I can find of how many periods? 
% Graph axes suggest it is ages 25-65, so I am guessing 40 periods of working life.
% Some later graphs from a deterministic model (Figure 5), seem to also show 10 periods of retirement (to age 75)
% This is not mentioned in the model description (pg 5, eqn 1).
% But then Figures 6&7 are from the stocastic model and also appear to show 10 periods of retirement
% The deterministic models underlying Figure 5 all have retirement income=70% of income in final working year.

%% Parameters
% Paper does not report all parameter values, it does
% state "I will solve the model and present most of my results for a single
% baseline set of parameter values, but will also present a summary of
% results for alternative choices of all parameter values so that readers
% can determine how sensitive my results are to changes in that parameter".
% In the NBER working paper version https://www.nber.org/system/files/working_papers/w5788/w5788.pdf
% Table 1 on page 52 of the pdf appear to show all these alternative
% parametrizations and to mark the baseline with a cross. I use these as
% the parameters here. (They match the parameter p that is reported in the
% paper, so seems fine.)

% Preferences
Params.delta=0.04; % Carroll (1997) pg 7
Params.beta=1/(1+Params.delta); % Carroll (1997) pg 5
Params.gamma=2; % pg 8 says there is CRRA utility with parameter 2

Params.r=0; % R=1+r, Carroll (1997) equations use R, but the working paper anyway actually reports r, rather than R (and I prefer to use r). Carroll (1997) pg 7 does say that r=0
Params.g=0.02;

% Exogenous shocks
Params.p=0.005; % Paper says "p=0.5 percent per year", so this is p=0.005 as a fraction
Params.sigma_lnZ=0.1;
Params.mew_lnZ=(-(Params.sigma_lnZ^2)/2)/(1-Params.p); % "mew_lnZ is chosen to make E_t V_{t+1}=0"; 
Params.sigma_lnN=0.1;
Params.mew_lnN=-(Params.sigma_lnN^2)/2;

% Retirement
% Description of model on pg 5 equation 1 does not have any retirement
% But figures 6 & 7 suggest there is retirement, and it seems plausible to
% me the model did contain it (otherwise all the discussion of buffer-stock
% vs retirement as savigns motives would be a bit awkward). So I put in ten
% periods of retirement, and make incomes fall to 70% of final working
% period income
Params.Jr=41; % retire at 65 (74 is last year of life)
Params.pension=0.7; % fraction of final working period income

Params.agej=1:1:N_j;


%% Grids
a_grid=10*linspace(0,1,n_a)'.^3; 
% Note: normally the bottom of the grid is implictly enforcing a borrowing constraint
% But in this model the fact that there is the possibility of a zero-income
% event (that v=0) in every period of (remaining) lifetime means that
% no-one would ever want to borrow anyway (as there would be a risk you get
% v=0 every period and end up with -Inf utility). That is, there is a
% 'natural borrowing limit' of zero, and so the zero at the bottom of the
% asset grid is not actually acting like a borrowing limit as no-one would
% ever want to breach it anyway.

% First, the V shock which is i.i.d. (page 6 of C1997 gives exact specification)
% Consists of a first stage, which Carroll (1997) calls Z
% "In solving the model the, the lognormal distributions were truncated at 
% three standard deviations from the mean"
farmertodaoptions.nSigmas=3;
[lnZ_grid,pi_Z]=discretizeAR1_FarmerToda(Params.mew_lnZ,0,Params.sigma_lnZ,n_e-1,farmertodaoptions);
pi_Z=pi_Z(1,:)';
% [sum(lnZ_grid.*pi_Z),Params.mew_lnZ] % should be equal, they are
% Now add in the zero-income event, and use Z instead of lnZ
V_grid=[0;exp(lnZ_grid)];
pi_V=[Params.p; (1-Params.p)*pi_Z]; % There is no difference between pi_lnZ and pi_Z

% sum(V_grid.*pi_V) 
% Not exactly equal to one because of the discretization, so I renormalize the grid so that we get exactly 1
V_grid=V_grid./sum(V_grid.*pi_V);
sum(V_grid.*pi_V) % This should be 1 (it is :)

% Now the permanent shock P
% Carroll (1997) has that P_t=G_t P_{t-1} N_t
% G=(1+g) is the growth factor for permanent labor income
% ln(N) ~ TN(-(sigma_lnN^2)/2, sigma_lnN^2)
% "The choice of mean for ln(N) was [chosen] to make E_t N_{t+1}=1"
% Can therefore tell that there is typo in Carroll (1997) as it is missing a minus
% on the expression for the mean of ln(N) (corrected in the version above)
% Similarly, does the truncated lognormal at three standard deviation from the mean for N_t.
% If we log this process we get: lnP_t= g + lnP_t-1 + lnN_t
% Since N_t is log-normal, we just get nlN_t ~ (-(sigma_lnN^2)/2, sigma_lnN^2), which is nice and easy
% Note: TN refers to truncated normal, as it is max of +- three std devs
% We need to use the extension of Farmer-Toda to age-dependent parameters to handle the income growth g and that there are permanent shocks
kirkbyoptions.nSigmas=3; 
kirkbyoptions.initialj1mewz=0; % WHAT DOES CARROLL USE TO START?
[lnP_grid_J, pi_P_J,jequaloneDistP,otheroutputs] = discretizeLifeCycleAR1_Kirkby(Params.g*ones(1,N_j)-(Params.sigma_lnN^2)/2,ones(1,N_j),Params.sigma_lnN*ones(1,N_j),n_z,N_j,kirkbyoptions);
% Note: Carroll (1997) does what you should do in a model with exogneous
% labor and  permament shocks, namely renormalize and then solve. We just
% take the lazy (but computationally costly) option here, on the plus side
% it makes it easy to modify this shock process.
P_grid_J=exp(lnP_grid_J);
% jequaloneDistP is the initial distribution

% To handle retirement I want P to be constant once reaching retirement
for jj=Params.Jr:N_j
    P_grid_J(:,jj)=P_grid_J(:,Params.Jr-1);
end
for jj=Params.Jr-1:N_j
    pi_P_J(:,:,jj)=eye(n_z,n_z); % the N_j-th one is not actually used for anything, but setting it anyway
end

% Put into VFI Toolkit notation
vfoptions.n_e=n_e;
vfoptions.e_grid=V_grid;
vfoptions.pi_e=pi_V;
simoptions.n_e=vfoptions.n_e;
simoptions.e_grid=vfoptions.e_grid;
simoptions.pi_e=vfoptions.pi_e;
vfoptions.z_grid_J=P_grid_J;
simoptions.z_grid_J=vfoptions.z_grid_J;
z_grid=P_grid_J(:,1); % just a placeholder
vfoptions.pi_z_J=pi_P_J;
simoptions.pi_z_J=vfoptions.pi_z_J;
pi_z=pi_P_J(:,:,1); % just a placeholder

%%
DiscountFactorParamNames={'beta'};

ReturnFn=@(aprime,a,P,V,r,gamma,agej,Jr,pension) ...
    Carroll1997_ReturnFn(aprime,a,P,V,r,gamma,agej,Jr,pension);

%% Solve value function
tic;
[V, Policy]=ValueFnIter_Case1_FHorz(0,n_a,n_z,N_j, [], a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);
vftime=toc

%% Paper talks a lot about consumption and savings, plot a bunch of the policies for this a different ages
FnsToEvaluate.consumption=@(aprime,a,P,V,r,agej,Jr,pension) Carroll1997_ConsumptionFn(aprime,a,P,V,r,agej,Jr,pension);
FnsToEvaluate.income=@(aprime,a,P,V,agej,Jr,pension) (agej<Jr)*(V*P)+(agej>=Jr)*(pension*P); %note, because of pi_P_J, P will be constant in retirement.
FnsToEvaluate.assets=@(aprime,a,P,V) a;
FnsToEvaluate.aprime=@(aprime,a,P,V) aprime;

ValuesOnGrid=EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1(Policy, FnsToEvaluate, Params, [], 0, n_a, n_z, N_j, [], a_grid, z_grid, [],simoptions);

% plot consumption for: median value of V, three values of P (on three graphs), as function of assets for a few different ages
fig1=figure(1);
subplot(3,1,1); plot(a_grid, ValuesOnGrid.consumption(:,ceil(2*n_z/3),7,1),a_grid, ValuesOnGrid.consumption(:,ceil(2*n_z/3),7,11),a_grid, ValuesOnGrid.consumption(:,ceil(2*n_z/3),7,21),a_grid, ValuesOnGrid.consumption(:,ceil(2*n_z/3),7,31),a_grid, ValuesOnGrid.consumption(:,ceil(2*n_z/3),7,41))
title('Consumption policy function: high permanent shock')
subplot(3,1,2); plot(a_grid, ValuesOnGrid.consumption(:,ceil(n_z/2),7,1),a_grid, ValuesOnGrid.consumption(:,ceil(n_z/2),7,11),a_grid, ValuesOnGrid.consumption(:,ceil(n_z/2),7,21),a_grid, ValuesOnGrid.consumption(:,ceil(n_z/2),7,31),a_grid, ValuesOnGrid.consumption(:,ceil(n_z/2),7,41))
title('Consumption policy function: median permanent shock')
subplot(3,1,3); plot(a_grid, ValuesOnGrid.consumption(:,ceil(1*n_z/3),7,1),a_grid, ValuesOnGrid.consumption(:,ceil(1*n_z/3),7,11),a_grid, ValuesOnGrid.consumption(:,ceil(1*n_z/3),7,21),a_grid, ValuesOnGrid.consumption(:,ceil(1*n_z/3),7,31),a_grid, ValuesOnGrid.consumption(:,ceil(1*n_z/3),7,41))
legend('Age 25','Age 35','Age 45','Age 55','Age 65')
title('Consumption policy function: low permanent shock')
% Comment: that curve near the zero assets is precautionary savings in action :)

% plot next-period assets for: median value of V, three values of P (on three graphs), as function of assets for a few different ages
fig2=figure(2);
subplot(3,1,1); plot(a_grid, ValuesOnGrid.aprime(:,ceil(2*n_z/3),7,1),a_grid, ValuesOnGrid.aprime(:,ceil(2*n_z/3),7,11),a_grid, ValuesOnGrid.aprime(:,ceil(2*n_z/3),7,21),a_grid, ValuesOnGrid.aprime(:,ceil(2*n_z/3),7,31),a_grid, ValuesOnGrid.aprime(:,ceil(2*n_z/3),7,41))
title('Next-period assets policy function: high permanent shock')
subplot(3,1,2); plot(a_grid, ValuesOnGrid.aprime(:,ceil(n_z/2),7,1),a_grid, ValuesOnGrid.aprime(:,ceil(n_z/2),7,11),a_grid, ValuesOnGrid.aprime(:,ceil(n_z/2),7,21),a_grid, ValuesOnGrid.aprime(:,ceil(n_z/2),7,31),a_grid, ValuesOnGrid.aprime(:,ceil(n_z/2),7,41))
title('Next-period assets policy function: median permanent shock')
subplot(3,1,3); plot(a_grid, ValuesOnGrid.aprime(:,ceil(1*n_z/3),7,1),a_grid, ValuesOnGrid.aprime(:,ceil(1*n_z/3),7,11),a_grid, ValuesOnGrid.aprime(:,ceil(1*n_z/3),7,21),a_grid, ValuesOnGrid.aprime(:,ceil(1*n_z/3),7,31),a_grid, ValuesOnGrid.aprime(:,ceil(1*n_z/3),7,41))
legend('Age 25','Age 35','Age 45','Age 55','Age 65')
title('Next-period assets policy function: low permanent shock')
% Comment: that curve near the zero assets is precautionary savings in action :)
% Less obvious in next-period assets than it was in consumption, you probably need to zoom in on bottom-right of graph to really see it



%% Define where agents are born
% I make them all born with zero assets, and median value of V shock
% Distribution over P shock is what was set above during the
% discretization process for P shock
jequaloneDist=zeros(n_a,n_z,n_e);
jequaloneDist(1,:,7)=jequaloneDistP';

% I am not going to calculate any aggregates so this won't matter but 
AgeWeightParamNames={'ageweight'};
Params.ageweight=ones(1,N_j)/N_j;

%% Simulate agent distribution
simoptions.parallel=5; % this is slower than the default (=2), but reduces the memory use
simoptions.verbose=1;
tic;
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy,0,n_a,n_z,N_j,pi_z,Params,simoptions);
agentdisttime=toc

%% Let's recreate a few graphs from Carroll (1997)
% Compute life-cycle profiles
AgeConditionalStats=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,[],Params,0,n_a,n_z,N_j,[],a_grid,z_grid,simoptions);

fig3=figure(3);
plot(24+Params.agej,AgeConditionalStats.consumption.Mean, 24+Params.agej,AgeConditionalStats.income.Mean, 24+Params.agej,AgeConditionalStats.assets.Mean)
legend('Consumption','Income','Assets')
title('Mean Life-cycle profiles')
% Note: the permanent shocks have a drift/growth-rate, hence why incomes grow




