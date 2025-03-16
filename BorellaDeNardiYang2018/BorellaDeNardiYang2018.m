% Borella, De Nardi & Yang (2018) - The Aggregate Implications of Gender and Marriage
% We solve Economy 4 which is a model with three kinds of households: single female, single male, and married couple.
% Because moving between these is possible, and the moves are exogenous, we
% need to think of one agent who changes 'type' based on an exogenous state.

% Two small notational changes from BDNY2018. The deterministic
% age-dependent component of ages is called kappa_j (they call it e) and
% the markov (AR(1)) shock is called z (they call it epsilon).

% Everything is nice and smooth. 

% Two small errors/typos in BDNY2018 paper: 
% 1) pg 11, eqn 2 tells us that n is 0 to 1, but when used in eqn 5 it should 
% be hours (so 5840 times as much). Effectively, eqn 5 (and 13) is missing 
% a constant parameter than multiplies n. I will refer to this missing
% parameter as 'potential hours'.
% 2) pg 21, the multivariate normal for the initial distribution has
% delta_a on left and right hand sides, it should only appear on the left
% (should not appear on the right).

% The paper is missing a few parameters but the authors were very good/nice about 
% sending them in an clean well-documented format when I requested them by 
% email. You can find them all below, specifically, the conditional
% survival probabilities, the parametrization of the initial distribution,
% the 'potential hours' that multiplies n in eqn 5 (and 14),
% and kappa_j (the kappa_j are shown in Figure 13, but they sent me the exact
% numbers).

% Reproduces Figures 9, 10, 11, 14 of BDNY2018 (the parts relating to model, those figures also contain data which is not reproduced here)
% Create additional new figures 1,2,3 to show some other parts of the model to help illustrate what is going on.

CreateFigures=1


%%
n_d=[51,51]; % Two decision variables, endogenous labor supply for female and male, respectively
n_a=701; % Endogenous state, asset holdings
n_z=[7,7,4]; % Three exogenous states, female labor productivity, male labor productivity, household-type
N_j=76; 
Params.agejshifter=24; % Ages 25 to 100.

% A quick explanation of the third exogenous state which takes four values.
% This is the 'four' household types. We have single males, single females,
% and married couples. Because the married couples can transition to being
% single male/female, we cannot set them up at permanent types. So we set
% up an exogenous state to control these transitions. Notice that these
% transition probabilities occur based on the conditional survival
% probabilities, and that the conditional survival probabilities depend on
% this state. As a result, instead of treating conditional survival
% probabilities as another discount factor (which is standard) we instead
% treat the conditional survival probabilities as transition probabilities
% of the exogenous state. We therefore need to add a fourth grid point to
% the exogenous state which signifies death. This is not particully
% computationally efficient, but it is easy enough (VFI Toolkit has no
% feature to permit the conditional survival probabilities to depend on the
% exogenous state which would be the better way to handle it.)

% The third exogenous state is called 'htype', and htype==1 is married
% couple, htype==2 is single female, htype==3 is single male, and hype==4
% is deceased.

%% Parameters
% Economy 4, so most parameters come from Table 9
% [some from Table 3 and Table 4]
% Paper is missing a few parameter values, BDNY were kind enough to send them by email. 
% Big thanks to them! Specifically, kappa_j (deterministic age-dependent component of 
% earnings), s_j (conditional survival probabilities), potential hours, and the 
% parameters that determine  the initial age distribution are those that are from email.

% Age j
Params.agej=1:1:N_j;
% actual age is agej+agejshifter
Params.age=Params.agej+Params.agejshifter;
% Retirement age
Params.Jr=66-Params.agejshifter; % Retire at age 66


% Preferences
Params.beta=0.961; % Discount factor
Params.gamma=2; % risk aversion coefficient
Params.omega=0.496; % Consumption weight

% Labor force participation fixed cost (only used for working age)
Params.participationcost_men=0.316;
Params.participationcost_women_single=0.376;
Params.participationcost_women_married=exp(0.0009.*Params.agej.^2-0.0236*Params.agej-1.0708)./(1+exp(0.0009.*Params.agej.^2-0.0236*Params.agej-1.0708));
% Participation cost for women: this formula comes from pg 25 of BDNY2018.
% Figure 14 plots it and the formula here does reproduce the figure.
% [I originally used age rather than agej in formula, but Fig 14 makes it clear that is should be agej]
if CreateFigures==1
    figure(14)
    plot(Params.age(1:Params.Jr-1),Params.participationcost_men*ones(1,Params.Jr-1))
    hold on
    plot(Params.age(1:Params.Jr-1),Params.participationcost_women_single*ones(1,Params.Jr-1))
    plot(Params.age(1:Params.Jr-1),Params.participationcost_women_married(1:Params.Jr-1))
    hold off
    ylabel('Participation Cost')
    xlabel('Age')
    legend('Men','Single Women','Married Women')
end

% Hours worked factor
Params.potentialhours=5840; % BDNY2018 by email tell me that the 'total potential hours' which multiplies n is5840 (=365 days*16 hrs per day)

% Prices
Params.r=0.04; % interest rate

% Retirement and Pensions
Params.tau_ss=0.038; % Social security tax rate on employees

% Earnings process, deterministic component
Params.kappa_j_men=[14.28697096, 14.92876666, 15.50469684, 16.01569213, 16.46435926, 16.85445477, 17.19077667, 17.47871754, 17.72400614, 17.93263075, 18.11042972, 18.26303369, 18.39575626, 18.51342232, 18.62027424, 18.72015966, 18.81615859, 18.91073972, 19.00573926, 19.10217112, 19.20062825, 19.30071226, 19.40137617, 19.50086844, 19.5965984, 19.68513928, 19.76237742, 19.82327991, 19.86221136, 19.87262189, 19.84732029, 19.77866833, 19.65858228, 19.47902022, 19.2315277, 18.90867856, 18.50342777, 18.01015922, 17.42479166, 16.74579195, 15.97365815, zeros(1,N_j-Params.Jr+1)]; % fill in zeros for retirement ages
Params.kappa_j_women=[10.23342169, 10.33615559, 10.43110283, 10.52023978, 10.605471, 10.6885492, 10.77108253, 10.85453538, 10.94021008, 11.0292408, 11.12258805, 11.22103341, 11.32520757, 11.43549369, 11.55209724, 11.67506067, 11.80415986, 11.93894804, 12.07874142, 12.22256703, 12.36923078, 12.51715371, 12.66452781, 12.8091787, 12.94859127, 13.07991699, 13.19993625, 13.30513432, 13.3916985, 13.45548729, 13.49216302, 13.49726402, 13.46622906, 13.3946316, 13.27811246, 13.11277571, 12.89515939, 12.6225513, 12.29306153, 11.90600625, 11.46174536, zeros(1,N_j-Params.Jr+1)]; % fill in zeros for retirement ages

% Earnings processes (from Table 4)
% Persistence of markov shocks
Params.rho_z_men=0.973;
Params.rho_z_women=0.963;
% Params.rho_z_menwomen=0.973;
Params.rho_z_married=0.972;
% Variance of markov innovations
Params.sigmasq_upsilon_men=0.016;
Params.sigmasq_upsilon_women=0.014;
% Params.sigmasq_upsilon_menwomen=0.021;
Params.sigmasq_upsilon_married=0.013;
% Initial variance
Params.sigmasq_init_z_men=0.128;
Params.sigmasq_init_z_women=0.128;
% Params.sigmasq_init_z_menwomen=0.128;
Params.sigmasq_init_z_married=0.128;

% Conditional Survival Probabilities [Just used to construct transition probabilities for htype]
% BDNY2018 Appendix A states "We model the probability of being alive at
% time t as a logit function" from HRS biennial data. They do state that s_j=1 for working ages
% I emailed authors and they sent me the following parameter values
Params.s_j_men=[ones(1,Params.Jr-1), 0.9708540312, 0.9654954844, 0.9601461089, 0.9549150054, 0.9498879234, 0.9451082213, 0.940578516, 0.9362672177, 0.9321064622, 0.928034319, 0.9239414182, 0.919726484, 0.9152684335, 0.9104409934, 0.9051246504, 0.8991750989, 0.8924468522, 0.884799902, 0.8760654515, 0.8661334404, 0.8548043058, 0.8419632417, 0.8275050364, 0.8113131369, 0.7934129312, 0.7737980124, 0.7526544571, 0.7301095265, 0.7065543154, 0.682371187, 0.6581025282, 0.6341793149, 0.6114001461, 0.5902426503, 0.5713174184];
Params.s_j_women=[ones(1,Params.Jr-1), 0.9867173292, 0.9839537557, 0.9810895031, 0.9781775739, 0.9752641812, 0.9723774046, 0.969525577, 0.966698489, 0.963864254, 0.9609913506, 0.9580178755, 0.9548816988, 0.9515039178, 0.947796448, 0.9436688848, 0.9390081383, 0.933691744, 0.9275915251, 0.9205495885, 0.9124369293, 0.9030503063, 0.8922328875, 0.8798241982, 0.86564477, 0.8496186806, 0.8316466875, 0.8117952836, 0.7901019966, 0.7668593618, 0.7423955654, 0.7172294449, 0.6918221763, 0.667022471, 0.6434068853, 0.6216936057];
Params.s_j_married_men=[ones(1,Params.Jr-1), 0.9857772492, 0.9826486379, 0.979351984, 0.975942647, 0.9724705996, 0.968967363, 0.9654432869, 0.9618879751, 0.9582658633, 0.9545406942, 0.9506398489, 0.9464884755, 0.9419901481, 0.9370354902, 0.9315104472, 0.9252723856, 0.9181662883, 0.9100301624, 0.9006674505, 0.8899198529, 0.8775450903, 0.863366277, 0.8472150339, 0.8289168868, 0.8084411651, 0.785749516, 0.7610170302, 0.7343962364, 0.706334129, 0.677301214, 0.6479487248, 0.6188217937, 0.5908262567, 0.5645208063, 0.5405648483];
Params.s_j_married_women=[ones(1,Params.Jr-1), 0.9935989401, 0.992044741, 0.9903554811, 0.9885517802, 0.9866536883, 0.9846739092, 0.9826155039, 0.980471028, 0.9782190319, 0.975836103, 0.973276773, 0.970490844, 0.9674117874, 0.9639603055, 0.960048177, 0.9555625464, 0.9503736337, 0.9443364884, 0.937270943, 0.92900752, 0.9193032442, 0.9079407056, 0.8946906353, 0.8793009575, 0.861614884, 0.8414655991, 0.8188646896, 0.7938284563, 0.7666627684, 0.7377541306, 0.7077193293, 0.6771431241, 0.6470122624, 0.6180256858, 0.5909978175];
% Do a plot of the survival probabilities, just to see them
if CreateFigures==1
    figure(1)
    plot(Params.age,Params.s_j_men)
    hold on
    plot(Params.age,Params.s_j_women)
    plot(Params.age,Params.s_j_married_men)
    plot(Params.age,Params.s_j_married_women)
    hold off
    ylabel('Age Conditional Survival Probabilities')
    xlabel('Age')
    legend('Single Men','Single Women','Married Men','Married Women')
end


%% Grids
a_grid=500000*linspace(0,1,n_a)'.^3; 
% I just had to guess a max, then adjust it till it was reasonable
% BDNY2018 are essentially using this grid in USD (can be nice when you go to data)

n1_grid=linspace(0,1,n_d(1))'; % can tell from size of 'participationcost' parameters that model is on 0 to 1 (not actual hours worked which is what figures show)
n2_grid=linspace(0,1,n_d(2))';
d_grid=[n1_grid; n2_grid];

% Discretization is an AR(1), but because there is an initial distribution
% (variance) which is not the asymptotic/unconditional variance, we are
% better (in sense grid will be more efficient) using a life-cycle AR(1) method.
% Because male and female exogenous shocks are independent, we will create them separately and then join the two
% Set the initial distribution
kirkbyoptions_men.initialj1sigmaz=sqrt(Params.sigmasq_init_z_men); % age j=1 is a normal distribution with mean zero and this standard deviation [you can also set the mean, but default is zero]
[z_grid_J_men, pi_z_J_men,jequaloneDistz_men,~] = discretizeLifeCycleAR1_KFTT(0,Params.rho_z_men,sqrt(Params.sigmasq_upsilon_men),n_z(2),N_j,kirkbyoptions_men);
kirkbyoptions_women.initialj1sigmaz=sqrt(Params.sigmasq_init_z_women); % age j=1 is a normal distribution with mean zero and this standard deviation [you can also set the mean, but default is zero]
[z_grid_J_women, pi_z_J_women,jequaloneDistz_women,~] = discretizeLifeCycleAR1_KFTT(0,Params.rho_z_women,sqrt(Params.sigmasq_upsilon_women),n_z(1),N_j,kirkbyoptions_women);
% Note, this is currently ln(z). We want the actual z
z_grid_J_men=exp(z_grid_J_men);
z_grid_J_women=exp(z_grid_J_women);
% Set up the htype
htype_grid=(1:1:4)';
pi_htype_J=zeros(4,4,N_j);% Build from the conditional survival probabilities
for jj=1:N_j
    % Married couples
    pi_htype_J(1,1,jj)=Params.s_j_married_men(jj)*Params.s_j_married_women(jj);
    pi_htype_J(1,2,jj)=1-Params.s_j_married_men(jj);
    pi_htype_J(1,3,jj)=1-Params.s_j_married_women(jj);
    % Single women
    pi_htype_J(2,2,jj)=Params.s_j_married_women(jj);
    pi_htype_J(2,4,jj)=1-Params.s_j_married_women(jj);
    % Single men
    pi_htype_J(3,3,jj)=Params.s_j_married_men(jj);
    pi_htype_J(3,4,jj)=1-Params.s_j_married_men(jj);
    % Note: since V for deceased is zero, what we put as transtion
    % probability from deceased to deceased (1 or 0) anyway irrelevant.
    pi_htype_J(4,4,jj)=1;
end

% Put together the exogenous shocks
z_grid_J=zeros(sum(n_z),N_j);
pi_z_J=zeros(prod(n_z),prod(n_z),N_j);
for jj=1:N_j
    z_grid_J(:,jj)=[z_grid_J_women(:,jj);z_grid_J_men(:,jj);htype_grid];
    pi_z_J(:,:,jj)=kron(pi_htype_J(:,:,jj),kron(pi_z_J_men(:,:,jj),pi_z_J_women(:,:,jj))); % note, kron in reverse order
end

%% Return Function and Discount Parameters
DiscountFactorParamNames={'beta'}; 
% Normally we would put the survival probabilities into the discount
% factor. But here we instead have the survival probabilities acting
% through the exogenous state transition probabilities as we need 

ReturnFn=@(n_f,n_m,aprime,a,z_f,z_m,htype,r,omega,gamma,tau_ss,participationcost_men,participationcost_women_single,participationcost_women_married,kappa_j_women,kappa_j_men,potentialhours)...
    BorellaDeNardiYang2018_ReturnFn(n_f,n_m,aprime,a,z_f,z_m,htype,r,omega,gamma,tau_ss,participationcost_men,participationcost_women_single,participationcost_women_married,kappa_j_women,kappa_j_men,potentialhours);

%% Solve value function
vfoptions=struct(); % use defaults
vfoptions.verbose=1;
vfoptions.lowmemory=1;
tic;
[V,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid_J,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
toc

%% Plot some policy functions to take a look
% I can't be bothered just now :)

%% You could take a look at value functions to see if better/worse off being single/married
% Not really important in this model, as all changes in martial status are completely exogenous anyway.



%% Age j=1 agent distribution
% BDNY2018 Appendix A
% Model joint distribution of initial assets and wage shocks at age 25 as a
% joint log-normal in the logarithm of the asset and wage shock
% (ln(a+delta_a) )    ( mu_a                  )
% (  ln z        ) ~ N(    mu_z       , Sigma_s)
% 
% delta_a is a shift parameter for assets to have only positive values to be able to take logs.
% For Economy 4 with couples it is joint normal of the two ln(z) for female
% and male and the assets (so tri-variate instead of bi-variate).
%
% NOTE: VERSION IN BDNY2018 ERRONEOUSLY PUTS delta_a ON RIGHT HAND SIDE AS WELL (mu_a+delta_a instead of mu_a as here)


% BDNY2018 don't report the parameters of this tri-variate log-normal
% distribution (except the initial variance of ln(z) which is in Table 4)
% I emailed and they kindly sent me the following.
% Joint Distribution at age 25														
% Couples				Variance/Covariance Matrix			
% 	                      mu	    delta_a			            log(wealth)	Prod. Shock	Prod. Shock Spouse
% log(wealth)	        11.35902	69155.39		log(wealth)	0.1578558	0.018727	0.0153577
% Prod. Shock	        0.0231509	0		      Prod. Shock	0.018727	0.091125	0.0153023
% Prod. Shock Spouse	0.0013973	0		 Prod. Shock Spouse	0.0153577	0.0153023	0.0777305

% Note that we already have grids on assets, and on productivity for female
% and male. So we just want to assign mass on these grids to fit this
% distribution.
% What we use is to allocate based on the joint-normal cdf (to the grids in logs).

% Let's first just create the mean and var-covar from the values in email.
delta_a=69155.39;
% Note that in codes, the ordering is assets, productivity shock female,
% productivity shock male. So have to swap female and male order from what 
% BDNY sent in the email. (I have assumed 'Spouse' refers to female)
jequal1_mu=[11.35902, 0.0013973, 0.0231509];
jequal1_sigma=[0.1578558, 0.0153577,0.018727; 0.0153577, 0.0777305, 0.0153023; 0.018727, 0.0153023, 0.091125];

% Get log(producitivity) grids for age j=1.
logz_f_jequal1grid=log(z_grid_J_women(:,1));
logz_m_jequal1grid=log(z_grid_J_men(:,1));

% Initial dist for on assets is on log(assets+delta_a)
tic;
jequaloneDist=MVNormal_ProbabilitiesOnGrid([log(a_grid+delta_a); logz_f_jequal1grid; logz_m_jequal1grid],jequal1_mu, jequal1_sigma, [n_a,n_z(1),n_z(2)]);
initdisttime=toc

% Take a look over the asset grid
if CreateFigures==1
    figure(2)
    plot(a_grid,cumsum(sum(sum(sum(jequaloneDist,4),3),2)))
    title('Initial dist: marginal cdf over assets')
    % Almost everyone starts with 0 to 200,000 USD
end

% Comment: Evaluating a Multivariate-Normal cdf is not
% computationally cheap/fast (for anything but a small grid).
% An obvious alternative to get code faster here would be to only put the
% initial dist onto the lower 1/2 (or somesuch) of the asset distribution. Following two
% lines demonstrate how you would do this
% jequaloneDist=zeros([n_a,n_z(1),n_z(2),'gpuArray');
% jequaloneDist(1:floor(n_a/2),:,:)=MVNormal_ProbabilitiesOnGrid([log(a_grid(1:floor(n_a/2))+delta_a); logz_f_jequal1grid; logz_m_jequal1grid],jequal1_mu, jequal1_sigma, [floor(n_a/2),n_z(1),n_z(2)]);

% Distribution is currently over [n_a,n_z(1),n_z(2)]
% Needs to be over [n_a,n_z(1),n_z(2), n_z(3)]

% To do this we have to define relative mass of the household types.
% BDNY2018 did this, as required for Fig 11

% Table 1 of BDNY2018 says 0.43 married, 0.07 single female, and 0.07
% single male (at age 25). Paper does not explictly say that these (initial) weights
% are were used for model, but BDNY2018 confirmed this to me by email.
mass_hhtypes=[0.43,0.07,0.07,0]; % married, single female, single male, deceased (note, these masses will change over time, this is just for j=1) 
mass_hhtypes=mass_hhtypes/sum(mass_hhtypes); % normalize total mass to one
jequaloneDist=jequaloneDist.*shiftdim(mass_hhtypes,-2);



%% Agent distribution
AgeWeightParamNames={'mewj'};
Params.mewj=(1:1:N_j)/N_j;
% BDNY2018 don't define this as it is not needed since only model output they look at is life-cycle profiles, 
% but toolkit requires it (will have no impact on results as long as you are only looking at life-cycle profiles are these
% are anyway conditional on age).

simoptions=struct(); % use default simoptions
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,pi_z_J,Params,simoptions);


%% Plot some life-cycle profiles
% We only want age-conditional means, by using simoptions.whichstats we can specify this and make things faster
simoptions.whichstats=[1,0,0,0,0,0,0]; % http://discourse.vfitoolkit.com/t/how-to-cut-runtimes-for-life-cycle-profiles-or-allstats-when-you-dont-need-all-the-stats-simoptions-whichstats/253


% Start with some basics, just look fraction of households of each type
FnsToEvaluate_hhtype.marriedhh=@(n_f,n_m,aprime,a,z_f,z_m,htype) (htype==1);
FnsToEvaluate_hhtype.singlefemalehh=@(n_f,n_m,aprime,a,z_f,z_m,htype) (htype==2);
FnsToEvaluate_hhtype.singlemalehh=@(n_f,n_m,aprime,a,z_f,z_m,htype) (htype==3);
FnsToEvaluate_hhtype.deceasedhh=@(n_f,n_m,aprime,a,z_f,z_m,htype) (htype==4);
AgeConditionalStats_hhtype=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate_hhtype,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid_J,simoptions);

if CreateFigures==1
    figure(3)
    plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats_hhtype.marriedhh.Mean./(1-AgeConditionalStats_hhtype.deceasedhh.Mean)) % fraction married
    hold on
    plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats_hhtype.singlefemalehh.Mean./(1-AgeConditionalStats_hhtype.deceasedhh.Mean)) % fraction single female
    plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats_hhtype.singlemalehh.Mean./(1-AgeConditionalStats_hhtype.deceasedhh.Mean)) % fraction single male
    hold off
    title('Fraction of (living) households of each type')
    legend('married','single female','single male')
end

% Note that if we just compute standard life-cycles then their
% interpretation is problematic (as, e.g., our model in a single female
% household imputes a male labor supply of zero, which is correct but awkward)

% Hence we will do each household type one by one, and compute stats conditional on this type.
FnsToEvaluate_married.participation_men=@(n_f,n_m,aprime,a,z_f,z_m,htype) (n_m>0);
FnsToEvaluate_married.participation_women=@(n_f,n_m,aprime,a,z_f,z_m,htype) (n_f>0);
FnsToEvaluate_married.hours_men=@(n_f,n_m,aprime,a,z_f,z_m,htype,potentialhours) potentialhours*n_m;
FnsToEvaluate_married.hours_women=@(n_f,n_m,aprime,a,z_f,z_m,htype,potentialhours) potentialhours*n_f;
FnsToEvaluate_married.earnings_men=@(n_f,n_m,aprime,a,z_f,z_m,htype,kappa_j_men,potentialhours) kappa_j_men*z_m*potentialhours*n_m; % BDNY2018 refer to these as e, epsilon, n
FnsToEvaluate_married.earnings_women=@(n_f,n_m,aprime,a,z_f,z_m,htype,kappa_j_women,potentialhours) kappa_j_women*z_f*potentialhours*n_f; % BDNY2018 refer to these as e, epsilon, n
FnsToEvaluate_married.assets=@(n_f,n_m,aprime,a,z_f,z_m,htype) a;
simoptions.SampleRestrictionFn=@(n_f,n_m,aprime,a,z_f,z_m,htype) (htype==1);
AgeConditionalStats_married=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate_married,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid_J,simoptions);

FnsToEvaluate_sfemale.participation=@(n_f,n_m,aprime,a,z_f,z_m,htype) (n_f>0);
FnsToEvaluate_sfemale.hours=@(n_f,n_m,aprime,a,z_f,z_m,htype,potentialhours) potentialhours*n_f;
FnsToEvaluate_sfemale.earnings=@(n_f,n_m,aprime,a,z_f,z_m,htype,kappa_j_women,potentialhours) kappa_j_women*z_f*potentialhours*n_f; % BDNY2018 refer to these as e, epsilon, n
FnsToEvaluate_sfemale.assets=@(n_f,n_m,aprime,a,z_f,z_m,htype) a;
simoptions.SampleRestrictionFn=@(n_f,n_m,aprime,a,z_f,z_m,htype) (htype==2);
AgeConditionalStats_sfemale=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate_sfemale,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid_J,simoptions);

FnsToEvaluate_smale.participation=@(n_f,n_m,aprime,a,z_f,z_m,htype) (n_m>0);
FnsToEvaluate_smale.hours=@(n_f,n_m,aprime,a,z_f,z_m,htype,potentialhours) potentialhours*n_m;
FnsToEvaluate_smale.earnings=@(n_f,n_m,aprime,a,z_f,z_m,htype,kappa_j_men,potentialhours) kappa_j_men*z_m*potentialhours*n_m; % BDNY2018 refer to these as e, epsilon, n
FnsToEvaluate_smale.assets=@(n_f,n_m,aprime,a,z_f,z_m,htype) a;
simoptions.SampleRestrictionFn=@(n_f,n_m,aprime,a,z_f,z_m,htype) (htype==3);
AgeConditionalStats_smale=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate_smale,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid_J,simoptions);


% Finally, a few things across all (living) households
FnsToEvaluate.participation=@(n_f,n_m,aprime,a,z_f,z_m,htype,potentialhours) ((n_f>0)+(n_m>0))/2*(htype==1)+(n_f>0)*(htype==2)+(n_m>0)*(htype==3); % note: divided by 2 for married HH (did this as that is what matches Fig 11 of BDNY2018)
FnsToEvaluate.hours=@(n_f,n_m,aprime,a,z_f,z_m,htype,potentialhours) potentialhours*(n_f+n_m)/2*(htype==1)+potentialhours*n_f*(htype==2)+potentialhours*n_m*(htype==3); % note: divided by 2 for married HH (did this as that is what matches Fig 11 of BDNY2018)
FnsToEvaluate.earnings=@(n_f,n_m,aprime,a,z_f,z_m,htype,kappa_j_women,kappa_j_men,potentialhours) (kappa_j_women*z_f*potentialhours*n_f+kappa_j_men*z_m*potentialhours*n_m)/2*(htype==1)+kappa_j_women*z_f*potentialhours*n_f*(htype==2)+kappa_j_men*z_m*potentialhours*n_m*(htype==3); % note: divided by 2 for married HH (did this as that is what matches Fig 11 of BDNY2018)
FnsToEvaluate.assets=@(n_f,n_m,aprime,a,z_f,z_m,htype) a/2*(htype==1)+a*(htype==2)+a*(htype==3); % would normally just be a, but need to do complicate it so I can divide by 2 for married couple
simoptions.SampleRestrictionFn=@(n_f,n_m,aprime,a,z_f,z_m,htype) (htype~=4);
AgeConditionalStats=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid_J,simoptions);


%% Figures of BDNY2018
% Let's do Figures 9, 10, 11 of BDNY2018
% These are the figures relating to Economy 4 (which is the model being done here)
% We only plot the 'half' of each Figure which is the model (they also plot data)

% To make it easy to compare to those in paper, I enforce similar xlim and ylim (ranges for the x and y axes)

% Figures 9 and 10 of BDNY2018
if CreateFigures==1
    figure(9)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:Params.Jr-1), AgeConditionalStats_smale.participation.Mean(1:Params.Jr-1), Params.agejshifter+(1:1:Params.Jr-1), AgeConditionalStats_sfemale.participation.Mean(1:Params.Jr-1), Params.agejshifter+(1:1:Params.Jr-1), AgeConditionalStats_married.participation_men.Mean(1:Params.Jr-1), Params.agejshifter+(1:1:Params.Jr-1), AgeConditionalStats_married.participation_women.Mean(1:Params.Jr-1))
    legend('Single Men','Single Women','Married Men','Married Women')
    title('Participation')
    ylim([0.2,1])
    xlim([25,65])
    subplot(2,1,2); plot(Params.agejshifter+(1:1:Params.Jr-1), AgeConditionalStats_smale.hours.Mean(1:Params.Jr-1), Params.agejshifter+(1:1:Params.Jr-1), AgeConditionalStats_sfemale.hours.Mean(1:Params.Jr-1), Params.agejshifter+(1:1:Params.Jr-1), AgeConditionalStats_married.hours_men.Mean(1:Params.Jr-1), Params.agejshifter+(1:1:Params.Jr-1), AgeConditionalStats_married.hours_women.Mean(1:Params.Jr-1))
    legend('Single Men','Single Women','Married Men','Married Women')
    title('Hours')
    % ylim([500,2500])
    xlim([25,65])

    figure(10)
    subplot(2,1,1); plot(Params.agejshifter+(1:1:Params.Jr-1), AgeConditionalStats_smale.earnings.Mean(1:Params.Jr-1), Params.agejshifter+(1:1:Params.Jr-1), AgeConditionalStats_sfemale.earnings.Mean(1:Params.Jr-1), Params.agejshifter+(1:1:Params.Jr-1), AgeConditionalStats_married.earnings_men.Mean(1:Params.Jr-1), Params.agejshifter+(1:1:Params.Jr-1), AgeConditionalStats_married.earnings_women.Mean(1:Params.Jr-1))
    legend('Single Men','Single Women','Married Men','Married Women')
    title('Labor Income')
    xlim([25,65])
    subplot(2,1,2); plot(Params.agejshifter+(1:1:N_j), AgeConditionalStats_smale.assets.Mean, Params.agejshifter+(1:1:N_j), AgeConditionalStats_sfemale.assets.Mean, Params.agejshifter+(1:1:N_j), AgeConditionalStats_married.assets.Mean)
    legend('Single Men','Single Women','Couples')
    title('Assets')
    xlim([25,100])
end

% Figure 11 of BDNY2018
if CreateFigures==1
    % Note: this Figure is 'per person', so all married couple numbers are divided by 2
    figure(11)
    subplot(2,2,1); plot(Params.agejshifter+(1:1:Params.Jr-1), AgeConditionalStats.participation.Mean(1:Params.Jr-1))
    legend('Economy 4')
    title('Participation')
    ylim([0.2,1])
    xlim([25,65])
    subplot(2,2,2); plot(Params.agejshifter+(1:1:Params.Jr-1), AgeConditionalStats.hours.Mean(1:Params.Jr-1))
    title('Hours')
    ylim([500,2500])
    xlim([25,65])
    subplot(2,2,3); plot(Params.agejshifter+(1:1:Params.Jr-1), AgeConditionalStats.earnings.Mean(1:Params.Jr-1))
    title('Labor Income')
    ylim([8000,65000])
    xlim([25,65])
    subplot(2,2,4); plot(Params.agejshifter+(1:1:N_j), AgeConditionalStats.assets.Mean)
    title('Assets')
    ylim([0,500000])
    xlim([25,100])
end

% If you run with CreateFigures=0, you can just load the following and create all the graphs later
% (I ran the main code on a headless server, then created graphs on my desktop)
save BDNY2018.mat
% load BDNY2018.mat

