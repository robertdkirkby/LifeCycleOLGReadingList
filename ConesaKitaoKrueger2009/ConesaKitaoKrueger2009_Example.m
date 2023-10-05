% Example implementing model of Conesa, Kitao & Krueger (2009) - Taxing Capital? Not a Bad Idea After All
% This example simply evaluates two different settings for the tax parameters, it does not attempt
% to find the optimal (CKK2009 original codes evaluate 55 differrent settings for the tax parameters and pick the best)

% A line I needed for running on the Server
addpath(genpath('./MatlabToolkits/'))

% One decision variable: labor supply (l)
% One endogenous state variable: assets (k)
% One stochastic exogenous state variable (eta, which represents idiosyncratic labor productivity)
% Age
% Permanent agent types

Params.J=81; % Number of period in life-cycle (represents ages 21 to 100)
Params.M=2; % Number of permanent types of agent.

% Grid sizes to use
n_l=21; % Labor supply
n_k=501; % Assets
n_eta=7; % Idiosyncratic labor productivity shocks
N_j=Params.J; % Number of periods in finite horizon.
N_i=Params.M; % Number of permanent types of agent.

% Remark: CKK2009, pg 33, say they discretize eta as AR(1) using 7 states.
% They do not appear to specify how they discretize. I have used Tauchen method with tauchen q=3. (Maybe they use Tauchen-Hussey?? Probably, Conesa & Krueger (1999) used Tauchen-Hussey)
% Note: equation (10) is 'incorrect' in sense that the representative firm ignores the payroll tax when deciding labor demand.

% Different model options
Model_AltUtilityFn=0; % 0 is non-seperable (utility of cons and leisure), 1 is seperable (1 forms the basis for Table 5)

% Remark: the parameters alpha and alpha_i are completely unrelated.



%% Set model parameters

% Preference parameters
Params.gamma_c=1.8; % Consumption share (CKK2009 call this gamma)
Params.gamma=4; % Risk aversion (CKK2009 call this sigma)
Params.beta=1.001; % Note that once this is combined with the conditional probability of survival it gives numbers that are typically less than 1.
% When using utility fn that is seperable in consumption and leisure (Model_AltUtilityFn=1)
if Model_AltUtilityFn==1
    % CKK2009, pg 43 gives
    Params.beta=0.972;
    Params.gamma_c=2;
    Params.gamma_l=3;
    Params.chi=1.92;
    % Not used
    Params.gamma=0;
else
    % Not used
    Params.gamma_l=0;
    Params.chi=0;
end
Params.Model_AltUtilityFn=Model_AltUtilityFn; % Need to pass it to the return fn later

% %Technology parameters
Params.A=1; % Aggregate productivity scale factor
Params.alpha=0.36; % Capital share
Params.delta=0.0833; % Depreciation rate

%People become economically active (j=1) at age 20, retire at 65, and max age is 100
Params.J=N_j; %=100-20+1
Params.Jr=46; % =65-20+1
Params.I_j=[1:1:Params.Jr,zeros(1,Params.J-Params.Jr)]; % Indicator for retirement

% Population growth of 1.1%
Params.n=0.011;

% Government tax and spending
% Income taxes are progressive and depend on three parameters
Params.tau_kappa0=0.258;
Params.tau_kappa1=0.768; % CKK2009 look for the optimal in range of 6 to 8
% tau_kappa2 is set as part of general eqm; is used to ensure government budget balance
Params.tau_ss=0.124; % CKK2009, Table 1; this is the social security payroll tax rate (half of this is paid by employer and half by employee)
% Social Security benefits SS are set as part of general equilibrium
% Social Security taxes are levied on labor income up to ybar which is set to 2.5 times average income
Params.ybardivYtarget=2.5; % ybar itself is set as part of general eqm
Params.tau_c=0.05; % CKK2009, Table 1; this is the tax rate on consumption
% G is set below as part of general eqm; and then kept constant across tax experiments
Params.GdivYtarget=0.17;
% tau_k is set below as finding the optimal tau_k is core purpose of the paper
% We do however want to first solve the model without capital tax as this
% is used by CKK2009 as a reference case when calculating the compensative
% equivalent variation (gain in welfare) acheived by setting the optimal
% tax rates.
Params.tau_k=0;

%Probability of surviving, conditional on being alive are taken from Bell & Miller (2002) (CKK2009, Table 1 on pg 32)
%F. C. Bell and M. L. Miller (2002), Life Tables for the United States Social Security Area 1900-2100, Actuarial Study No. 116, Office of the Chief Actuary
%https://www.ssa.gov/oact/NOTES/pdf_studies/study116.pdf
%Table 8 — Period Probabilities of Death Within One Year (qx) at Selected Exact Ages, by Sex and Calendar Year (Cont.)
%The raw data from there is
%          Sex and Exact Age
%    |  Male                                                Female
%Year| [0 30 60 65 70 100]                                  [0 30 60 65 70 100]
%2005| [0.00642,0.00151,0.01209,0.02032,0.03054,0.36873]    [0.00538,0.00068,0.00769,0.01289,0.01972,0.32090]
%I just take the numbers for Males, and then set my actual values based on a linear interpolation of the data.
dj_temp=interp1([0,30,60,65,70,100],[0.00642,0.00151,0.01209,0.02032,0.03054,0.36873],0:1:100,'linear');
Params.sj=1-dj_temp(20:100);
Params.sj(1)=1;
Params.sj(end)=0;
% I use linear interpolation to fill in the numbers inbetween those reported by Bell & Miller (2002).
% I have aditionally imposed that the prob of death at age 20 be zero and that prob of death at age 100 is one. 
% Note: CKK2009 do not specify that they use the year 2005 numbers, nor if they just use the numbers for males, 
% nor how they interpolate these.


%Age profile of productivity (based on lifetime profile of earnings, Hansen 1993), Epsilon_j
%First the raw data of Hansen 1993:
%Table II. Weights assigned to age-sex
%         Males       Females
%Age      Weight      Weight
%16-19    0.56        0.52
%20-24    0.78        0.69
%25-34    1.14        0.89
%35-44    1.37        0.90
%45-54    1.39        0.87
%55-64    1.33        0.84
%65 +     0.89        0.66
%14-17    0.56        0.52
%14-19    0.56        0.52
%18-24    0.78        0.69
%25-44    1.24        0.89
%45-64    1.37        0.86
%epsilon_j is set based on the first entries of the data for males
epsilon_j=zeros(Params.J,1);
epsilon_j(1:4)=0.78; epsilon_j(5:14)=1.14; epsilon_j(15:24)=1.37; 
epsilon_j(25:34)=1.39; epsilon_j(35:44)=1.33; epsilon_j(45:65)=0.89; % I1989 pg 316 refers to 'interpolated to inbetween years', given that Hansen (1993) provides numbers for all years in the form of 'buckets' it is unclear exactly what this means.
Params.epsilon_j=epsilon_j; %/mean(epsilon_j(1:44)); % I1998 pg 316, "normalized to average unity over the 44 working years".
% Idle note: this is same process used by Imrohoroglu (1989) and Conesa & Krueger (1999).
% Imrohoroglu bans retirees from working. One way to do this is just to
% set epsilon_j=0 for retirees. Rather than do this via epsilon_j it is
% done by I_j which in practice is multiplied by epsilon_j.
Params.I_j=[ones(Params.Jr-1,1);zeros(Params.J-Params.Jr+1,1)];
% Note that because CKK2009 uses 80 periods the last 15 entries of Params.epsilon_j are actually zero anyway.

% Permanent types of ability (CKK2009, pg 33)
Params.sigmasq_alpha=0.14; % CKK2009, Table 1 on pg 32
Params.sigma_alpha=sqrt(Params.sigmasq_alpha);
Params.alpha_i=[exp(-Params.sigma_alpha);exp(Params.sigma_alpha)];
Params.prob_alpha_i=0.5; % CKK2008, pg 33 set population weights for each permanent type to half.

% Idiosyncratic shocks to labour productivity
Params.rho_eta=0.98; % CKK2009 call this rho.
Params.sigmasq_eta=0.0289;

%%
l_grid=linspace(0,1,n_l)';

k_grid=28*(linspace(0,1,n_k).^3)'; % I originally use 20 as max, but a number of points stayed at that level, so increased to 25, then again to 28
% CKK2009: use grid for assets from 0 to 75 (Their GRID.f90)
% Note: 'top51gridpointmass' below shows zero, so noone tries to go near
% k=28 (noone reaches any of the top 51 grid points on k (top 51 out of n_k points))

% CKK2009, pg 33. log(eta) is AR(1) with persistence parameter rho_eta and
% unconditional variance sigmasq_eta. They discretize with seven states.
% They do not explain how they discretize, have assumed Tauchen method with Tauchen q equal to 3.
Params.tauchenq=3;
[eta_grid, pi_eta]=discretizeAR1_Tauchen(0,Params.rho_eta,sqrt(Params.sigmasq_eta*(1-Params.rho_eta^2)),n_eta,Params.tauchenq); %[states, transmatrix], transmatix is (z,zprime)
eta_grid=exp(eta_grid);

%% Get into format for VFI Toolkit
d_grid=l_grid;
a_grid=k_grid;
z_grid=eta_grid;
pi_z=pi_eta;

n_d=n_l;
n_a=n_k;
n_z=n_eta;

%% General equilibrium, and initial parameter values for the general eqm parameters

% Note: since we use Cobb-Douglas prodn fn and only look at stationary eqm
% we don't actually need w here as it is just a function of r in any case.
GEPriceParamNames={'r','w','G','SS','ybar','tau_kappa2','Tr_beq'};
% % Originally I used the following, but the GE did not converge with enough accuracy
% Params.r=0.06; % interest rate on assets
% Params.w=1; % wages
% Params.G=0.17; % Government spending
% Params.SS=0.4; % Social security benefits
% Params.ybar=2.5*1.4; % ybar is multiple-of-avg-income up to which social security payroll taxes are paid (target is 2.5*Y, the 1.4 is a guess for Y)
% Params.tau_kappa2=0.122; % set to balance government budget deficit
% Params.Tr_beq=0.1; % Accidental bequests

% Using CMA-ES algorithm (heteroagentoptions.fminalgo=4) instead of the
% default (heteroagentoptions.fminalgo=1, matlab fminsearch()) we get
% convergence, following is a new initial guess that is essentially
% equilbirium that it finds.
Params.r=0.1265; % interest rate on assets
Params.w=0.8673; % wages
Params.G=0.2519; % Government spending
Params.SS=0.8757; % Social security benefits
Params.ybar=3.7046; % ybar is multiple-of-avg-income up to which social security payroll taxes are paid (target is 2.5*Y, the 1.4 is a guess for Y)
Params.tau_kappa2=2.8288; % set to balance government budget deficit
Params.Tr_beq=0.0689; % Accidental bequests

%% Now, create the return function 
DiscountFactorParamNames={'beta','sj'};

ReturnFn=@(l,kprime,k,eta,r,w,tau_k, tau_ss,ybar,tau_c,SS,chi,gamma_c,gamma_l,gamma,epsilon_j,alpha_i,I_j,A,alpha,delta, Tr_beq, tau_kappa0, tau_kappa1, tau_kappa2,Model_AltUtilityFn) ConesaKitaoKrueger2009_ReturnFn(l,kprime,k,eta,r,w,tau_k, tau_ss,ybar,tau_c,SS,chi,gamma_c,gamma_l,gamma,epsilon_j,alpha_i,I_j,A,alpha,delta, Tr_beq, tau_kappa0, tau_kappa1, tau_kappa2,Model_AltUtilityFn)

%% Now solve the value function iteration problem, just to check that things are working before we go to General Equilbrium

fprintf('Test ValueFnIter \n')
tic;
[V, Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,N_i, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames);
toc

% max(max(max(max(Policy))))<n_a % Double check that never try to leave top of asset grid.
% sum(sum(sum(sum(Policy==n_a))))

%% Distribution across the two permanent types of agent
Params.alpha_i_dist=[Params.prob_alpha_i; 1-Params.prob_alpha_i];
PTypeDistParamNames={'alpha_i_dist'};

%% Calculate the population distribution across the ages (population at time
% t is made stationary by dividing it by \Pi_{i=1}^{t} (1+n_{i}) (product of all
% growth rates up till current period))
% CKK2009 uses mewj, but has different notation for sj (calls it psi_j). The VFI Toolkit needs to give
% it a name so that it can be automatically used when calculating model outputs.
Params.mewj=ones(Params.J,1);
for jj=2:Params.J
    Params.mewj(jj)=Params.mewj(jj-1)*(1/(1+Params.n))*Params.sj(jj-1);
end
Params.mewj=Params.mewj/sum(Params.mewj); % normalize to measure one
Params.mewj=Params.mewj'; % Age weights must be a row vector.

AgeWeightsParamNames={'mewj'}; % Many finite horizon models apply different weights to different 'ages'; eg., due to survival rates or population growth rates.


%% Initial distribution of agents at birth (j=1)
% Note that this is the same for both agent types (is independent of agent type)
jequaloneDist=zeros(n_a,n_z); 
jequaloneDist(1,1+floor(n_z/2))=1; % CKK2009, pg 35 "all newborn households start with zero assets and average labor producitivity"

%% Test
fprintf('Test StationaryDist \n')
tic;
StationaryDist=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightsParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_z,N_j,N_i,pi_z,Params);
toc


%% Set up the General Equilibrium conditions (on assets/interest rate, assuming a representative firm with Cobb-Douglas production function)

% Steady State Aggregates
% Aggregate assets K
FnsToEvaluate.K = @(d,aprime,a,z) a;
% Aggregate effective labour supply (in efficiency units), CKK2009 calls this N
FnsToEvaluate.N = @(d,aprime,a,z,I_j,alpha_i,epsilon_j) I_j*alpha_i*epsilon_j*z*d;
% Tr, accidental bequest transfers % The divided by (1+n) is due to population growth and because this is Tr_{t+1}
FnsToEvaluate.Tr = @(d,aprime,a,z,sj,n) (1-sj)*aprime/(1+n); 
% Consumption
FnsToEvaluate.C = @(d,aprime,a,z,r,w,tau_k, tau_ss,ybar,tau_c,SS,epsilon_j,alpha_i,I_j,Tr_beq, tau_kappa0, tau_kappa1, tau_kappa2) ConesaKitaoKrueger2009_ConsumptionFn(d,aprime,a,z,r,w,tau_k, tau_ss,ybar,tau_c,SS,epsilon_j,alpha_i,I_j,Tr_beq, tau_kappa0, tau_kappa1, tau_kappa2);
% Labor income tax revenue
FnsToEvaluate.TaxRevenue = @(d,aprime,a,z,r,tau_ss,ybar,epsilon_j,alpha_i,I_j,A,alpha,delta, tau_kappa0, tau_kappa1, tau_kappa2) ConesaKitaoKrueger2009_LaborIncomeTaxRevenueFn(d,aprime,a,z,r,tau_ss,ybar,epsilon_j,alpha_i,I_j,A,alpha,delta, tau_kappa0, tau_kappa1, tau_kappa2);
% Mass of retirees (note that this could be calculated directly as is exogenous)
FnsToEvaluate.MassRetired = @(d,aprime,a,z,I_j) 1-I_j;
% Total social security tax revenue
FnsToEvaluate.SSrevenue = @(d,aprime,a,z,tau_ss,ybar,I_j,alpha_i,epsilon_j) tau_ss*min(I_j*alpha_i*epsilon_j*z*d,ybar);

% General Equilibrium Equations
% Recall that GEPriceParamNames={'r','w','G','SS','ybar','tau_kappa2','Tr_beq'};
% interest rate equals marginal product of capital (net of depreciation) (eqn 9)
GeneralEqmEqns.CapitalMarket = @(r,K,N,A,alpha,delta) r-(A*(alpha)*(K^(alpha-1))*(N^(1-alpha))-delta); % Rate of return on assets is related to Marginal Product of Capital
GeneralEqmEqns.LaborMarket = @(w,K,N,A,alpha) w-(A*(1-alpha)*(K^(alpha))*(N^(-alpha))); % wage is related to Marginal Product of (effective) Labor
% G as fraction of Y is equal to GdivYtarget
GeneralEqmEqns.GovSizeTarget = @(G,N,K,GdivYtarget,A,alpha) G-GdivYtarget*(A*(K^(alpha))*(N^(1-alpha)));
% SS benefits are equal to the amount raised by the social security payroll taxes (eqn 11)
GeneralEqmEqns.SSbalance = @(SS,MassRetired,SSrevenue) SS*MassRetired-SSrevenue;
% ybar divided by Y is equal to ybardivYtarget
GeneralEqmEqns.ybarTarget = @(ybar,K,N,ybardivYtarget,A,alpha) ybar-ybardivYtarget*(A*(K^(alpha))*(N^(1-alpha)));
% Government budget balance (role of tau_kappa2 is via TaxRevenue)
GeneralEqmEqns.GovBudgetBalance = @(G,tau_k,r,K,Tr,TaxRevenue,tau_c,C) G-(tau_k*r*(K+Tr)+TaxRevenue+tau_c*C); % Government budget balance
% Accidental bequests (eqn 12)
GeneralEqmEqns.AccBeqBalance = @(Tr_beq,Tr) Tr_beq-Tr;

%% Test
disp('Test AggVars')
AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1_PType(StationaryDist, Policy, FnsToEvaluate, Params,n_d,n_a,n_z,N_j,N_i,d_grid, a_grid, z_grid);

%% Calculate the general equilibrium without capital tax
fprintf('GE without capital tax \n')

% heteroagentoptions.fminalgo=4 % 4 is CMA-ES algorithm; the default is 1 which is fminsearch() but that doesn't seem to work for CKK2009 model

heteroagentoptions.fminalgo=5;
heteroagentoptions.fminalgo5.howtoupdate={...
    'CapitalMarket','r',0,0.1;...  % CaptialMarket GE condition will be positive if r is too big, so subtract
    'LaborMarket','w',0,0.1;... % LaborMarket GE condition will be positive if w is too big, so subtract
    'GovSizeTarget','G',0,0.1;... % GovSizeTarget GE condition will be positive if G is too big, so subtract
    'SSbalance','SS',0,0.1;... % SSbalance GE condition will be positive if SS is too big, so subtract
    'ybarTarget','ybar',0,0.1;... % ybarTarget GE condition will be positive if ybar is too big, so subtract
    'GovBudgetBalance','tau_kappa2',1,0.1;... % GovBudgetBalance GE condition will be negative if tau_kappa2 is too big, so add
    'AccBeqBalance','Tr_beq',0,0.1;... % AccBeqBalance GE condition will be positive if Tr_beq is too big, so subtract
    };

heteroagentoptions.verbose=1;
tic;
[p_eqm,~,GeneralEqmConditions]=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d, n_a, n_z, N_j, N_i, [], pi_z, d_grid, a_grid, z_grid,jequaloneDist, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightsParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptions);
GEtime=toc
Params.r=p_eqm.r;
Params.w=p_eqm.w;
Params.G=p_eqm.G;
Params.SS=p_eqm.SS;
Params.ybar=p_eqm.ybar;
Params.tau_kappa2=p_eqm.tau_kappa2;
Params.Tr_beq=p_eqm.Tr_beq;

% We need to reload some of this later
save ./CKK2009_notaxeqm.mat p_eqm GeneralEqmConditions Params GEtime


% Comment: Matlab fminsearch() was unable to find the general eqm. Had to
% switch to CMA-ES algorithm which is slower but much more robust. From the
% look of CKK2009 codes they use a shooting algorithm, what they do is use a guess
% for prices, use this to create a new guess for prices (e.g. newr=marginal
% product of capital net depreciation), and then update prices as
% 0.8*'guess for prices'+0.2*'new guess for prices'. I have since then
% switched my codes over to fminalgo=5 which uses the same kind of
% shooting algorithms to find the general equilibrium as they did.

%%
% Average Hours Worked
FnsToEvaluate_Extra.HoursWorked = @(d,aprime,a,z) d; % Hours worked

% Calculate some aggregate variables in economy without capital taxes
[V, Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,N_i, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames);
StationaryDist=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightsParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_z,N_j,N_i,pi_z,Params);
AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1_PType(StationaryDist, Policy, FnsToEvaluate, Params, n_d, n_a, n_z,N_j,N_i, d_grid, a_grid, z_grid);
AggVarsExtra=EvalFnOnAgentDist_AggVars_FHorz_Case1_PType(StationaryDist, Policy, FnsToEvaluate_Extra, Params, n_d, n_a, n_z,N_j,N_i, d_grid, a_grid, z_grid);
AvgHoursWorked_0=AggVarsExtra.HoursWorked.Mean;
N_0=AggVars.N.Mean;
K_0=AggVars.K.Mean;
Y_0=Params.A*(K_0^Params.alpha)*(N_0^(1-Params.alpha));
C_0=AggVars.C.Mean;
% Calculate the welfare without capital tax (needed to later calculate the CEV; pg 36 of CKK2009)
Wcl_0=sum(Params.(PTypeDistParamNames{1}).*[V.ptype001(1,1+floor(n_z/2),1);V.ptype002(1,1+floor(n_z/2),1)]);
% We will also need the following for CEV decomposition calculations
FnsToEvaluate_Decisions.c = @(d,aprime,a,z,r,w,tau_k, tau_ss,ybar,tau_c,SS,epsilon_j,alpha_i,I_j,Tr_beq, tau_kappa0, tau_kappa1, tau_kappa2) ConesaKitaoKrueger2009_ConsumptionFn(d,aprime,a,z,r,w,tau_k, tau_ss,ybar,tau_c,SS,epsilon_j,alpha_i,I_j,Tr_beq, tau_kappa0, tau_kappa1, tau_kappa2);
FnsToEvaluate_Decisions.l=@(d,aprime,a,z) d; % hours worked (not leisure)
ValuesOnGrid=EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1_PType(StationaryDist, Policy, FnsToEvaluate_Decisions, Params,n_d,n_a,n_z,N_j,N_i,d_grid, a_grid, z_grid);
Policy_c0.ptype001=ValuesOnGrid.c.ptype001;
Policy_c0.ptype002=ValuesOnGrid.c.ptype002;
Policy_l0.ptype001=ValuesOnGrid.l.ptype001;
Policy_l0.ptype002=ValuesOnGrid.l.ptype002;
Policy_0=Policy;

V_0=V; % I use this for a double-check later on
StationaryDist_0=StationaryDist; % I use this for a double-check later on

% Life cycle profiles for graphs
% Assets K
FnsToEvaluate_LifeCycleProfiles.assets = @(d,aprime,a,z) a;
% Hours worked
FnsToEvaluate_LifeCycleProfiles.hoursworked = @(d,aprime,a,z) d;
% Consumption
FnsToEvaluate_LifeCycleProfiles.consumption = @(d,aprime,a,z,r,w,tau_k, tau_ss,ybar,tau_c,SS,epsilon_j,alpha_i,I_j,Tr_beq, tau_kappa0, tau_kappa1, tau_kappa2) ConesaKitaoKrueger2009_ConsumptionFn(d,aprime,a,z,r,w,tau_k, tau_ss,ybar,tau_c,SS,epsilon_j,alpha_i,I_j,Tr_beq, tau_kappa0, tau_kappa1, tau_kappa2);
% Total income tax revenue
FnsToEvaluate_LifeCycleProfiles.taxrevenue = @(d,aprime,a,z,r,tau_k,tau_ss,ybar,epsilon_j,alpha_i,I_j,A,alpha,delta,Tr_beq, tau_kappa0, tau_kappa1, tau_kappa2) ConesaKitaoKrueger2009_TotalIncomeTaxRevenueFn(d,aprime,a,z,r,tau_k,tau_ss,ybar,epsilon_j,alpha_i,I_j,A,alpha,delta,Tr_beq, tau_kappa0, tau_kappa1, tau_kappa2);
AgeConditionalStats_0=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate_LifeCycleProfiles,Params,n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,z_grid);


% Want to check that the asset grid has appropriate maximum
PolicyValues=PolicyInd2Val_Case1_FHorz_PType(Policy,n_d,n_a,n_z,N_j,d_grid,a_grid);
% plot(shiftdim(PolicyValues(2,:,:),1)-a_grid*ones(1,N_z))
% This plot should be less than zero over to the right-hand end of plot

fprintf('Following are to check that max value of grid on assets is large enough: \n')
fprintf('  Number of points that choose the max assets is %i for PType 1 and %i for PType 2 \n', sum(sum(sum(sum(Policy_0.ptype001(2,:,:,:)==n_a)))),  sum(sum(sum(sum(Policy_0.ptype002(2,:,:,:)==n_a)))) )
fprintf('  Mass of agents on top 50  asset points is %8.12f for PType 1 and %8.12f for PType 2 \n', sum(sum(sum(sum(StationaryDist.ptype001(end-49:end,:,:))))),  sum(sum(sum(sum(StationaryDist.ptype002(end-49:end,:,:))))) )
fprintf('  Mass of agents on top 100 asset points is %8.12f for PType 1 and %8.12f for PType 2 \n', sum(sum(sum(sum(StationaryDist.ptype001(end-99:end,:,:))))),  sum(sum(sum(sum(StationaryDist.ptype002(end-99:end,:,:))))) )
fprintf('  Mass of agents on top 200 asset points is %8.12f for PType 1 and %8.12f for PType 2 \n', sum(sum(sum(sum(StationaryDist.ptype001(end-199:end,:,:))))),  sum(sum(sum(sum(StationaryDist.ptype002(end-199:end,:,:))))) )
fprintf('  Mass of agents on top 300 asset points is %8.12f for PType 1 and %8.12f for PType 2 \n', sum(sum(sum(sum(StationaryDist.ptype001(end-299:end,:,:))))),  sum(sum(sum(sum(StationaryDist.ptype002(end-299:end,:,:))))) )
fprintf('  Mass of agents on top 500 asset points is %8.12f for PType 1 and %8.12f for PType 2 \n', sum(sum(sum(sum(StationaryDist.ptype001(end-499:end,:,:))))),  sum(sum(sum(sum(StationaryDist.ptype002(end-499:end,:,:))))) )
% Note to self: seems okay, the top 50 asset points are empty for both types.

save ./SavedOutput/CKK2009_notaxeqm.mat p_eqm GeneralEqmConditions Params V Policy_0 StationaryDist AggVars AggVarsExtra AgeConditionalStats_0 Policy_c0 Policy_l0 PolicyValues a_grid d_grid


%% In this example we will not attempt to compute the optimal tax. Instead we just try out two different values  for the tax parameters and see which is best
% CKK2009 original codes actually do this for 55 different vectors of the
% tax parameters. In an ideal world, you would set it up as a constrained
% maximiztion (choose tax parameters to maximize social welfare function
% subject to model general equilibrium)

% We will try two values of tau_k
tau_k_vec=[0.3,0.4];

% Set an initial guess for GE params
Params.r=0.1042;
Params.w=0.9236;
Params.G=0.2683;
Params.SS=0.8762;
Params.ybar=3.9452;
Params.tau_kappa2=2.82;
Params.Tr_beq=0.0771;
% Set up general equilibrium conditions
heteroagentoptions.fminalgo=5;
heteroagentoptions.fminalgo5.howtoupdate={...
    'CapitalMarket','r',0,0.3;...  % CaptialMarket GE condition will be positive if r is too big, so subtract
    'LaborMarket','w',0,0.3;... % LaborMarket GE condition will be positive if w is too big, so subtract
    'GovSizeTarget','G',0,0.3;... % GovSizeTarget GE condition will be positive if G is too big, so subtract
    'SSbalance','SS',0,0.3;... % SSbalance GE condition will be positive if SS is too big, so subtract
    'ybarTarget','ybar',0,0.3;... % ybarTarget GE condition will be positive if ybar is too big, so subtract
    'GovBudgetBalance','tau_kappa2',1,0.3;... % GovBudgetBalance GE condition will be negative if tau_kappa2 is too big, so add
    'AccBeqBalance','Tr_beq',0,0.3;... % AccBeqBalance GE condition will be positive if Tr_beq is too big, so subtract
    };
heteroagentoptions.toleranceGEprices_percent=10^(-4);
heteroagentoptions.verbose=1;

%% First, try out tau_k=0.3
Params.tau_k=tau_k_vec(1);
% Calculate the general equilibrium and associated level of welfare
[p_eqm1,~,GeneralEqmConditions]=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d, n_a, n_z, N_j, N_i, [], pi_z, d_grid, a_grid, z_grid,jequaloneDist, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightsParamNames, PTypeDistParamNames, GEPriceParamNames, heteroagentoptions);
Params.r=p_eqm1.r;
Params.w=p_eqm1.w;
Params.G=p_eqm1.G;
Params.SS=p_eqm1.SS;
Params.ybar=p_eqm1.ybar;
Params.tau_kappa2=p_eqm1.tau_kappa2;
Params.Tr_beq=p_eqm1.Tr_beq;

[V, ~]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,N_i, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames);
% Calculate the welfare
welfare_tauk1=sum(Params.(PTypeDistParamNames{1}).*[V.ptype001(1,1+floor(n_z/2),1);V.ptype002(1,1+floor(n_z/2),1)]);

%% Second, try out tau_k=0.4
Params.tau_k=tau_k_vec(2);
% Calculate the general equilibrium and associated level of welfare
[p_eqm2,~,GeneralEqmConditions]=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d, n_a, n_z, N_j, N_i, [], pi_z, d_grid, a_grid, z_grid,jequaloneDist, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightsParamNames, PTypeDistParamNames, GEPriceParamNames, heteroagentoptions);
Params.r=p_eqm2.r;
Params.w=p_eqm2.w;
Params.G=p_eqm2.G;
Params.SS=p_eqm2.SS;
Params.ybar=p_eqm2.ybar;
Params.tau_kappa2=p_eqm2.tau_kappa2;
Params.Tr_beq=p_eqm2.Tr_beq;

[V, ~]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,N_i, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames);
% Calculate the welfare
welfare_tauk2=sum(Params.(PTypeDistParamNames{1}).*[V.ptype001(1,1+floor(n_z/2),1);V.ptype002(1,1+floor(n_z/2),1)]);


%% Which of our two tax rates is the optimal (welfare maximizing)
[opt_tau_k,optindex]=max([welfare_tauk1,welfare_tauk2]);
if optindex==1
    p_eqm_opt=p_eqm1;
elseif optindex==2
    p_eqm_opt=p_eqm2;
end

%% Table 2: Calculate a bunch of things about aggregate variables in economy with optimal taxation (as in the better of the two tax rates that we tried)
clear heteroagentoptions
heteroagentoptions.fminalgo=1 % 4 is CMA-ES algorithm; the default is 1 which is fminsearch() but that doesn't seem to work for CKK2009 model
% Using fminalgo=1 here as starting from what is essentially the GE anyway (just want to recalculate it as a double-check)

% Use what we found for the optimal tax rate
Params.tau_k=tau_k_vec(optindex);
Params.r=p_eqm_opt.r;
Params.w=p_eqm_opt.w;
Params.G=p_eqm_opt.G;
Params.SS=p_eqm_opt.SS;
Params.ybar=p_eqm_opt.ybar;
Params.tau_kappa2=p_eqm_opt.tau_kappa2;
Params.Tr_beq=p_eqm_opt.Tr_beq;
% Probably unnecessary, but just resolve the general eqm
[p_eqm,~,GeneralEqmConditions]=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d, n_a, n_z, N_j, N_i, [], pi_z, d_grid, a_grid, z_grid,jequaloneDist, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightsParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptions);
Params.r=p_eqm.r;
Params.w=p_eqm.w;
Params.G=p_eqm.G;
Params.SS=p_eqm.SS;
Params.ybar=p_eqm.ybar;
Params.tau_kappa2=p_eqm.tau_kappa2;
Params.Tr_beq=p_eqm.Tr_beq;

[V, Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,N_i, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames);
StationaryDist=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightsParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_z,N_j,N_i,pi_z,Params);
AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1_PType(StationaryDist, Policy, FnsToEvaluate, Params, n_d, n_a, n_z,N_j,N_i, d_grid, a_grid, z_grid);
AggVarsExtra=EvalFnOnAgentDist_AggVars_FHorz_Case1_PType(StationaryDist, Policy, FnsToEvaluate_Extra, Params, n_d, n_a, n_z,N_j,N_i, d_grid, a_grid, z_grid);

Wcl_optimal=sum(Params.(PTypeDistParamNames{1}).*[V.ptype001(1,1+floor(n_z/2),1); V.ptype002(1,1+floor(n_z/2),1)]); % (1,4,1,:) is no assets (1), average labor productivity (1+floor(n_z/2)), age one (1), both permanent types (:)

% Average Hours Worked
% Total Labor Supply, N
% Capital Stock, K
% Output, Y
% Aggregate Consumption, C
AvgHoursWorked_optimal=AggVarsExtra.HoursWorked.Mean;
N_optimal=AggVars.N.Mean;
K_optimal=AggVars.K.Mean;
Y_optimal=Params.A*(K_optimal^Params.alpha)*(N_optimal^(1-Params.alpha));
C_optimal=AggVars.C.Mean;

CEV=(Wcl_optimal/Wcl_0)^(1/(Params.gamma_c*(1-Params.gamma)))-1;
% CEV_Old=(Wcl_optimal/Wcl_0)^((1/Params.gamma_c)*(1-Params.gamma))-1; % This is how I first misread eqn (27) on page 36 of CKK2009. Their parentheses are unclear.

% We will also need the following for CEV decomposition calculations
ValuesOnGrid=EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1_PType(StationaryDist, Policy, FnsToEvaluate_Decisions, Params,n_d,n_a,n_z,N_j,N_i,d_grid, a_grid, z_grid);
Policy_coptimal.ptype001=ValuesOnGrid.c.ptype001;
Policy_coptimal.ptype002=ValuesOnGrid.c.ptype002;
Policy_loptimal.ptype001=ValuesOnGrid.l.ptype001;
Policy_loptimal.ptype002=ValuesOnGrid.l.ptype002;
Policy_optimal=Policy;
% I use the following as a double-check on CEV decompositions
V_optimal=V;

% Some life-cycle profiles for Figure 1
AgeConditionalStats_optimal=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate_LifeCycleProfiles,Params,n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,z_grid);

% Calculate the percentage changes
AvgHoursWorked_pch=100*(AvgHoursWorked_optimal-AvgHoursWorked_0)/AvgHoursWorked_0;
N_pch=100*(N_optimal-N_0)/N_0;
K_pch=100*(K_optimal-K_0)/K_0;
Y_pch=100*(Y_optimal-Y_0)/Y_0;
C_pch=100*(C_optimal-C_0)/C_0;

% Table 2
FID = fopen('./ConesaKitaoKrueger2009_Table2.tex', 'w');
fprintf(FID, 'Changes in Aggregate Variables in the Optimal Tax System \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lc} \n \\hline \\hline \n');
fprintf(FID, ' Variable & Change in percent  \\\\ \\hline \n');
fprintf(FID, ' Average Hours Worked & %8.2f \\\\ \n',AvgHoursWorked_pch);
fprintf(FID, ' Total Labor Supply, N & %8.2f \\\\ \n',N_pch);
fprintf(FID, ' Capital Stock, K & %8.2f \\\\ \n',K_pch);
fprintf(FID, ' Output, Y & %8.2f \\\\ \n',Y_pch);
fprintf(FID, ' Aggregate Consumption, C & %8.2f \\\\ \n',C_pch);
fprintf(FID, ' CEV & %8.2f \\\\ \n',100*CEV);
fprintf(FID, '\\hline \n \\end{tabular*} \n');
% fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
% fprintf(FID, 'Note: . \n');
% fprintf(FID, '}} \\end{minipage}');
fclose(FID);

% Check if anyone is trying to leave the top of the asset grid
typenames=fieldnames(StationaryDist);
top51gridpointmass_type1=sum(sum(StationaryDist.(typenames{1}),2),3); % sum over eta and age
top51gridpointmass_type1=sum(top51gridpointmass_type1(end-50:end)); % top 51 grid points on assets
top51gridpointmass_type2=sum(sum(StationaryDist.(typenames{2}),2),3); % sum over eta and age
top51gridpointmass_type2=sum(top51gridpointmass_type2(end-50:end)); % top 51 grid points on assets

% save ./SavedOutput/CKK2009_optimalcapitaltaxeqm.mat Params p_eqm GeneralEqmConditions V Policy_optimal StationaryDist AggVars AggVarsExtra AgeConditionalStats_optimal Policy_coptimal Policy_loptimal CEV Wcl_optimal Wcl_0 top51gridpointmass_type1 top51gridpointmass_type2

%% Table 3: Calculate a few decompositions of the CEV
% See Footnote 17 on pg 36 of CKK2009 for the formulae that define the following concepts.

% Need to use the initial Params
load ./SavedOutput/CKK2009_notaxeqm.mat Params

%%
% As a double-check, report the V_0 based welfare measures. We should get the same thing when we calculate 'W_c0CEV_Cl0' for CEV_C=0 (as it is just
% evaluating the same value fn, but using consumption and leisure policy functions)
W_c0l0_directfromV0=sum(Params.(PTypeDistParamNames{1}).*[V_0.ptype001(1,1+floor(n_z/2),1); V_0.ptype002(1,1+floor(n_z/2),1)]); % (1,4,1,:) is no assets (1), average labor productivity (1+floor(n_z/2)), age one (1), both permanent types (:).

% As a double-check, report the V_optimal based welfare measures. We should get the same thing when we calculate 'W_cstarlstar' as part of finding CEV_L (as it is just
% evaluating the same value fn, but using consumption and leisure policy functions)
W_cstarlstar_directfromVoptimal=sum(Params.(PTypeDistParamNames{1}).*[V_optimal.ptype001(1,1+floor(n_z/2),1); V_optimal.ptype002(1,1+floor(n_z/2),1)]); % (1,4,1,:) is no assets (1), average labor productivity (1+floor(n_z/2)), age one (1), both permanent types (:).

fprintf('Calculating W_c0l0_directfromV0 we get %8.8f, we should find the same we look at W_c0CEV_Cl0 for CEV_C=0 \n', W_c0l0_directfromV0)
% Consumption
% Total
chooseCEV_C=@(CEV_C) chooseCEV_C_function(CEV_C, Policy_optimal, Policy_0, Policy_coptimal, Policy_loptimal, Policy_c0, Policy_l0,n_d,n_a,n_z,N_j,pi_z, Params,PTypeDistParamNames,DiscountFactorParamNames); %,F0);
fprintf('Calculating CEV_C: \n')
[CEV_C,thisshouldbezero_C]=fminsearch(chooseCEV_C,0);
fprintf(' \n')
save ./SavedOutput/CKK2009_CEV_C.mat CEV_C thisshouldbezero_C
% Level
CEV_CL=C_optimal/C_0-1; % See Footnote 17 on pg 36 of CKK2009
% Distribution
CEV_CD=CEV_C-CEV_CL;

fprintf('Calculating W_cstarlstar_directfromVoptimal we get %8.8f, we should find the same we look at W_cstarlstar as part of finding CEV_L \n', W_cstarlstar_directfromVoptimal)
% Leisure
% Total
chooseCEV_L=@(CEV_L) chooseCEV_L_function(CEV_L, Policy_optimal, Policy_0, Policy_coptimal, Policy_loptimal, Policy_c0, Policy_l0,n_d,n_a,n_z,N_j,pi_z, Params,PTypeDistParamNames,DiscountFactorParamNames);
fprintf('Calculating CEV_L: \n')
[CEV_L,thisshouldbezero_L]=fminsearch(chooseCEV_L,0);
fprintf(' \n')
save ./SavedOutput/CKK2009_CEV_L.mat CEV_L thisshouldbezero_L
% Level
CEV_LL=(1-AvgHoursWorked_optimal)/(1-AvgHoursWorked_0)-1; % See Footnote 17 on pg 36 of CKK2009; note that 'leisure=1-hours worked' (and expectations operator is linear, so E(1-l)=1-E(l)).
% Distribution
CEV_LD=CEV_L-CEV_LL;

fprintf('We have CEV=%8.4f, CEV_C=%8.4f, CEV_L=%8.4f \n',CEV, CEV_C, CEV_L)
fprintf('When calculating CEV_C=%8.4f we had thisshouldbezero_C=%8.4f \n',CEV_C, thisshouldbezero_C)
fprintf('When calculating CEV_L=%8.4f we had thisshouldbezero_L=%8.4f \n',CEV_L, thisshouldbezero_L)

% Note: I don't really need to calculate CEV_L above, I could just use following
% formula to get it. I do it so that following formula can then be used as a double-check of results.
CEV_test_thisshouldbezero=(1+CEV)-(1+CEV_C)*(1+CEV_L)

% Table 3
FID = fopen('.ConesaKitaoKrueger2009_Table3.tex', 'w');
fprintf(FID, 'Decomposition of Welfare \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}llr} \n \\hline \\hline \n');
fprintf(FID, ' \\multicolumn{2}{l}{Total change (in percent)} & %8.2f  \\\\ \\hline \n', 100*CEV);
fprintf(FID, ' \\multirow{3}{*}{Consumption} & Total & %8.2f \\\\ \n', 100*CEV_C); 
fprintf(FID, ' & Level & %8.2f \\\\ \n', 100*CEV_CL); 
fprintf(FID, ' & Distribution & %8.2f \\\\ \n', 100*CEV_CD); 
fprintf(FID, ' \\multirow{3}{*}{Leisure} & Total & %8.2f \\\\ \n', 100*CEV_L); 
fprintf(FID, ' & Level & %8.2f \\\\ \n', 100*CEV_LL); 
fprintf(FID, ' & Distribution & %8.2f \\\\ \n', 100*CEV_LD); 
fprintf(FID, '\\hline \n \\end{tabular*} \n');
% fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
% fprintf(FID, 'Note: . \n');
% fprintf(FID, '}} \\end{minipage}');
fclose(FID);
% Note: 'in percent', hence all the 100*

%% Life-cycle profiles for Figure 1

save ./SavedOutput/CKK.mat

figure(1)
% Assets
subplot(2,2,1); plot(20:1:100,AgeConditionalStats_0.assets.ptype001.Mean,20:1:100,AgeConditionalStats_optimal.assets.ptype001.Mean)
hold on
plot(20:1:100,AgeConditionalStats_0.assets.ptype002.Mean,20:1:100,AgeConditionalStats_optimal.assets.ptype002.Mean)
hold off
title('Assets')
% Work hours
subplot(2,2,2); plot(20:1:100,AgeConditionalStats_0.hoursworked.ptype001.Mean,20:1:100,AgeConditionalStats_optimal.hoursworked.ptype001.Mean)
hold on
plot(20:1:100,AgeConditionalStats_0.hoursworked.ptype002.Mean,20:1:100,AgeConditionalStats_optimal.hoursworked.ptype002.Mean)
hold off
title('Hours worked')
legend('PType 1, tax0','PType 1, taxopt','PType 2, tax0','PType 2, taxopt')
% Consumption
subplot(2,2,3); plot(20:1:100,AgeConditionalStats_0.consumption.ptype001.Mean,20:1:100,AgeConditionalStats_optimal.consumption.ptype001.Mean)
hold on
plot(20:1:100,AgeConditionalStats_0.consumption.ptype002.Mean,20:1:100,AgeConditionalStats_optimal.consumption.ptype002.Mean)
hold off
title('Consumption')
% Total income taxes paid
subplot(2,2,4); plot(20:1:100,AgeConditionalStats_0.taxrevenue.ptype001.Mean,20:1:100,AgeConditionalStats_optimal.taxrevenue.ptype001.Mean)
hold on
plot(20:1:100,AgeConditionalStats_0.taxrevenue.ptype002.Mean,20:1:100,AgeConditionalStats_optimal.taxrevenue.ptype002.Mean)
hold off
title('Income Tax Revenue')
% saveas(gcf,'./SavedOutput/Graphs/CKK2009_Figure1.png')



