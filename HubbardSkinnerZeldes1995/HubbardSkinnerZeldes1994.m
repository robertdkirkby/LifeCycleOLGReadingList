%% This code implements the baseline model of
% Hubbard, Skinner & Zeldes (1994) - The importance of precautionary motives in explaining individual and aggregate saving
% The model detailed in this Carnegie-Rochester Conference Series on Public Policy paper appears to also be exactly the same as that used in
% Hubbard, Skinner & Zeldes (1994) - Expanding the Life-Cycle Model: Precautionary Saving and Public Policy
% Hubbard, Skinner & Zeldes (1995) - Precautionary Saving and Social Insurance
% published in AER:P&P and JPE respecitively. Although I do not replicate all the results of these other two.

% If you want the codes to run faster, reduce simoptions.numbersims below
% (the number of household simulations used when creating panel data simulations)

%% State space variables of the model
% Age
% One endogenous variable: assets
% Two stochastic exogenous variables: income shock, medical expense shock
% Three permanent types: 
%  pt1: high-school dropout
%  pt2: high-school
%  pt3: university

%% Declare some parameter values

Params.J=80; % Ages 21 to 100 inclusive.
WorkingAgeVec=21:1:65;
RetiredAgeVec=66:100;
Params.age=[WorkingAgeVec,RetiredAgeVec];

%
n_a=501;
maxa=5*10^5;
n_z=[21,21]; % 21,21 % income, medical. HSZ1994 use [9,9] (HSV1994, pg 112)
Names_i={'NoHighSchool','HighSchool','College'}; % Number of fixed types
N_j=Params.J; % Number of periods in finite horizon
Params.q=3; % For tauchen method. They used q=2.5: pg 112 of pdf, "range between 2.5 standard deviations (of the unconditional distribution) above and below..."

Params.gamma=3;
% Gamma plays three roles:      
        % gamma is coefficient of relative risk aversion
        % 1/gamma is the intertemporal elasticity of substitution of consumption
        % gamma+1 is the coefficient of relative prudence (Kimball, 1990)

Params.delta=0.03; % rate at which agents discount future
Params.beta=1/(1+Params.delta);
Params.r=0.03; % interest rate on assets

% Add a parameter which just indicates the permanent type (for use when simulating panel data)
Params.hhtype=[1,2,3];

% Mortality rates (probability of survival)
% Table A.1
Params.dj=[0.00058, 0.00061, 0.00062, 0.00064, 0.00065, 0.00067, 0.00069, 0.00070, 0.00072, 0.00075, 0.00078, 0.00082, 0.00086, 0.00091, 0.00098, 0.00105, 0.00115, 0.00128, 0.00144, 0.00161, 0.00180, 0.00200,...
    0.00221, 0.00242, 0.00266, 0.00292, 0.00320, 0.00349, 0.00380, 0.00413, 0.00450, 0.00490, 0.00533,...
    0.00533, 0.00581, 0.00632, 0.00689, 0.00749, 0.00811, 0.00878, 0.00952, 0.01033, 0.01121, 0.01223,...
    0.01332, 0.01455, 0.01590, 0.01730, 0.01874, 0.02028, 0.02203, 0.02404, 0.02623, 0.02863, 0.03128,...
    0.03432, 0.03778, 0.04166, 0.04597, 0.05078, 0.05615, 0.06214, 0.06885, 0.07631, 0.08455, 0.09352,...
    0.10323, 0.11367, 0.12484, 0.13677, 0.14938, 0.16289, 0.17721, 0.19234, 0.20828, 0.22418, 0.23980, 0.25495, 0.26937, 0.28284]; % conditional probability of death
Params.sj=1-Params.dj; % Conditional survival probabilities. Act as a discount rate.

% Wage (Earnings) incomes, W
% Table A.2
Params.DeterministicWj_working.NoHighSchool=[10993-196*WorkingAgeVec+2167*((WorkingAgeVec.^2)/100)-2655*((WorkingAgeVec.^3)/10000)];
Params.DeterministicWj_working.HighSchool=[-11833+1421*WorkingAgeVec-137*((WorkingAgeVec.^2)/100)-2186*((WorkingAgeVec.^3)/10000)]; 
Params.DeterministicWj_working.College=[72270-5579*WorkingAgeVec+19200*((WorkingAgeVec.^2)/100)-18076*((WorkingAgeVec.^3)/10000)];
% Table A.3
Params.DeterministicWj_retired.NoHighSchool=[17861-133*RetiredAgeVec]; 
Params.DeterministicWj_retired.HighSchool=[29733-245*RetiredAgeVec]; 
Params.DeterministicWj_retired.College=[48123-429*RetiredAgeVec];
% Is not mentioned in paper, but based on Figures 3a-c is clear that
% earnings from ages 90+ are simply frozen at their age 89 level.
Params.DeterministicWj_retired.NoHighSchool(25:35)=Params.DeterministicWj_retired.NoHighSchool(24)*ones(size(Params.DeterministicWj_retired.NoHighSchool(25:35)));
Params.DeterministicWj_retired.HighSchool(25:35)=Params.DeterministicWj_retired.HighSchool(24)*ones(size(Params.DeterministicWj_retired.HighSchool(25:35)));
Params.DeterministicWj_retired.College(25:35)=Params.DeterministicWj_retired.College(24)*ones(size(Params.DeterministicWj_retired.College(25:35)));
Params.DeterministicWj.NoHighSchool=[Params.DeterministicWj_working.NoHighSchool, Params.DeterministicWj_retired.NoHighSchool];
Params.DeterministicWj.HighSchool=[Params.DeterministicWj_working.HighSchool, Params.DeterministicWj_retired.HighSchool];
Params.DeterministicWj.College=[Params.DeterministicWj_working.College, Params.DeterministicWj_retired.College];
% % Compare to Fig A.1(a-c): they won't be exactly the same as drop the year fixed effects, but eyeballing suggests they are fine.
% plot(21:1:100, Params.DeterministicWj.NoHighSchool, 21:1:100, Params.DeterministicWj.HighSchool, 21:1:100, Params.DeterministicWj.College)
% Stochastic Wj (Table A.4): u_it, an AR(1) in regressions in paper
Params.w_rho=[0.955; 0.946; 0.955];
Params.w_sigmasqepsilon=[0.033; 0.025; 0.016];
Params.w_sigmasqu=Params.w_sigmasqepsilon./(1-Params.w_rho.^2);
Params.w_sigmasqupsilon=[0.040; 0.021; 0.014]; % Estimated from PSID but not used in model.
% Bottom of pg 103 (45th pg): Wit=exp(log(DeterministicWj)-0.5*sigmasqu+uit)
% "This specification ensures that when we compare the certainty case with
% the earnings uncertainty case, we hold the age-conditional means of earnings constant."

% Medical Expenses, M
% Table A.5
Params.DeterministicMj_working=[5.373+0.073*WorkingAgeVec-0.753*((WorkingAgeVec.^2)/1000); 4.749+0.106*WorkingAgeVec-1.084*((WorkingAgeVec.^2)/1000); 4.816+0.109*WorkingAgeVec-1.090*((WorkingAgeVec.^2)/1000)];
Params.DeterministicMj_retired=[14.441-0.200*RetiredAgeVec+1.288*((RetiredAgeVec.^2)/1000); 11.371-0.101*RetiredAgeVec+0.540*((RetiredAgeVec.^2)/1000); 9.553-0.054*RetiredAgeVec-0.297*((RetiredAgeVec.^2)/1000)];
Params.DeterministicMj=[Params.DeterministicMj_working, Params.DeterministicMj_retired]; % is in logs
Params.m_rho=0.901*ones(3,1);
Params.m_sigmasqepsilon=[0.175; 0.156; 0.153];
Params.m_sigmasqmew=Params.m_sigmasqepsilon./(1-Params.m_rho.^2);
Params.m_sigmasqomega=[0.220; 0.220; 0.220];
% Consumption Floor
Params.Cbar=7000; % (middle of pg. 111)

%% Create grids and discretize exogenous shocks
% Stochastic Wj (Table A.4): u_it, an AR(1) in regressions in paper
[z1_grid.NoHighSchool,pi_z1.NoHighSchool]=discretizeAR1_Tauchen(0,Params.w_rho(1),sqrt(Params.w_sigmasqepsilon(1)),n_z(1),Params.q);
[z1_grid.HighSchool,pi_z1.HighSchool]=discretizeAR1_Tauchen(0,Params.w_rho(2),sqrt(Params.w_sigmasqepsilon(2)),n_z(1),Params.q);
[z1_grid.College,pi_z1.College]=discretizeAR1_Tauchen(0,Params.w_rho(3),sqrt(Params.w_sigmasqepsilon(3)),n_z(1),Params.q);
% Bottom of pg 103 (45th pg): Wit=exp(log(DeterministicWj)-0.5*sigmasqu+uit)
% "This specification ensures that when we compare the certainty case with
% the earnings uncertainty case, we hold the age-conditional means of earnings constant."

% Medical Expenses, M
[z2_grid.NoHighSchool,pi_z2.NoHighSchool]=discretizeAR1_Tauchen(0,Params.m_rho(1),sqrt(Params.m_sigmasqepsilon(1)),n_z(2),Params.q);
[z2_grid.HighSchool,pi_z2.HighSchool]=discretizeAR1_Tauchen(0,Params.m_rho(2),sqrt(Params.m_sigmasqepsilon(2)),n_z(2),Params.q);
[z2_grid.College,pi_z2.College]=discretizeAR1_Tauchen(0,Params.m_rho(3),sqrt(Params.m_sigmasqepsilon(3)),n_z(2),Params.q);

%% Grids
z_grid.NoHighSchool=[z1_grid.NoHighSchool; z2_grid.NoHighSchool];
z_grid.HighSchool=[z1_grid.HighSchool; z2_grid.HighSchool];
z_grid.College=[z1_grid.College; z2_grid.College];

pi_z.NoHighSchool=kron(pi_z2.NoHighSchool,pi_z1.NoHighSchool); % note, kron() in reverse order
pi_z.HighSchool=kron(pi_z2.HighSchool,pi_z1.HighSchool); % note, kron() in reverse order
pi_z.College=kron(pi_z2.College,pi_z1.College); % note, kron() in reverse order

a_grid=linspace(0,maxa,n_a)'; % Could probably do better by adding more grid near zero

%% Now, create the return function 
DiscountFactorParamNames={'beta','sj'};

ReturnFn=@(aprime,a,W_z1,M_z2,age,gamma,r,Cbar,DeterministicWj, w_sigmasqu, DeterministicMj, m_sigmasqmew)...
    HubbardSkinnerZeldes1994_ReturnFn(aprime,a,W_z1,M_z2,age,gamma,r,Cbar,DeterministicWj, w_sigmasqu, DeterministicMj, m_sigmasqmew);

%% Now solve the value function iteration problem
vfoptions.verbose=1;
tic;
[V, Policy]=ValueFnIter_Case1_FHorz_PType(0,n_a,n_z,N_j,Names_i, 0, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames,vfoptions);
toc


%% Create some outputs to replicate those of HubbardSkinnerZeldes1994

% Note: The annual growth rate of the population is assumed to be 1%.
Params.n=0.01;

simoptions.parallel=1;
simoptions.numbersims=10^5; % Number of simulations on which panel data (and life-cycle profile) results will be based.

PTypeDistParamNames={'ptypeweights'};
Params.ptypeweights=[0.22,0.56,0.22]'; % Hubbard, Skinner & Zeldes (1994) do not appear to report these 
                             % weights/probabilities anywhere in either of the two papers. 
                             % I have 'ballparked' them based on alternative sources for 1984 as
                             % fractions of US population (pg 1 of http://www.russellsage.org/sites/all/files/chartbook/Educational%20Attainment%20and%20Achievement.pdf )


%% Generate the kinds of outputs that are reported in Hubbard, Skinner & Zeldes (1994) in Tables 1,2,3
% 
% % Note: The annual growth rate of the population is assumed to be 1%.
% Params.n=0.01;
% 
% simoptions.numbersims=10^3; % Number of simulations on which panel data (and life-cycle profile) results will be based.

%% Create an initial distribution from which agents are drawn/born.

% I have been unable to find any mention in the paper of how the initial
% distribution from which agents are born is determined. I therefore simply
% assume that they are all born with zero assets, the stationary distribution 
% on income shocks, and the mean shock value on medical expenses.
% (Paper does specify that all 'newborn' are age 21; ie. first period.)
% More tricky is what weight to attach to each of the permanent types, this
% remains unclear. HubbardSkinnerZeldes1994 do not appear to report these
% numbers in the paper at all.

z1staty.NoHighSchool=ones(size(z1_grid.NoHighSchool))/length(z1_grid.NoHighSchool);
for ii=1:1000
    z1staty.NoHighSchool=pi_z1.NoHighSchool'*z1staty.NoHighSchool;
end
z1staty.HighSchool=ones(size(z1_grid.HighSchool))/length(z1_grid.HighSchool);
for ii=1:1000
    z1staty.HighSchool=pi_z1.HighSchool'*z1staty.HighSchool;
end
z1staty.College=ones(size(z1_grid.College))/length(z1_grid.College);
for ii=1:1000
    z1staty.College=pi_z1.College'*z1staty.College;
end

% PTypeDist=[0.25,0.25,0.5]';
PTypeDist=Params.(PTypeDistParamNames{1});
InitialDist.NoHighSchool=zeros([n_a,n_z]); 
InitialDist.NoHighSchool(1,:,ceil(n_z(2)/2))=z1staty.NoHighSchool'*PTypeDist(1);
InitialDist.HighSchool=zeros([n_a,n_z]); 
InitialDist.HighSchool(1,:,ceil(n_z(2)/2))=z1staty.HighSchool'*PTypeDist(2);
InitialDist.College=zeros([n_a,n_z]); 
InitialDist.College(1,:,ceil(n_z(2)/2))=z1staty.College'*PTypeDist(3);

%% Life-cycle profiles

FnsToEvaluate.Assets = @(aprime,a,z1,z2) a; 
FnsToEvaluate.Earnings= @(aprime,a,z1,z2,DeterministicWj,w_sigmasqu) exp(log(DeterministicWj)-0.5*w_sigmasqu+z1);
FnsToEvaluate.age = @(aprime,a,z1,z2,age) age;
FnsToEvaluate.Income = @(aprime,a,z1,z2,r,DeterministicWj,w_sigmasqu) r*a+exp(log(DeterministicWj)-0.5*w_sigmasqu+z1);  
FnsToEvaluate.Consumption= @(aprime,a,z1,z2,r,Cbar,DeterministicWj,w_sigmasqu, DeterministicMj) HubbardSkinnerZeldes1994_ConsumptionFn(aprime,a,z1,z2,r,Cbar,DeterministicWj,w_sigmasqu, DeterministicMj);
FnsToEvaluate.TR = @(aprime,a,z1,z2,r,Cbar,DeterministicWj,w_sigmasqu, DeterministicMj) HubbardSkinnerZeldes1994_TRFn(aprime,a,z1,z2,r,Cbar,DeterministicWj,w_sigmasqu, DeterministicMj);
FnsToEvaluate.AssetIncomeRatio = @(aprime,a,z1,z2,r,DeterministicWj,w_sigmasqu) a/(r*a+exp(log(DeterministicWj)-0.5*w_sigmasqu+z1));  
FnsToEvaluate.SavingsRate = @(aprime,a,z1,z2,r,DeterministicWj,w_sigmasqu) (aprime-a)/(r*a+exp(log(DeterministicWj)-0.5*w_sigmasqu+z1));  
FnsToEvaluate.GrossSavings = @(aprime,a,z1,z2) aprime-a;  
FnsToEvaluate.ptype = @(aprime,a,z1,z2,hhtype) hhtype;

simoptions.lifecyclepercentiles=0; % Just mean and median, no percentiles.
SimLifeCycleProfiles=SimLifeCycleProfiles_FHorz_Case1_PType(InitialDist,PTypeDistParamNames,Policy, FnsToEvaluate,Params,0,n_a,n_z,N_j,Names_i,0,a_grid,z_grid,pi_z, simoptions);

% % Figure: Assets
% figure(1)
% plot(Params.age,SimLifeCycleProfiles.NoHighSchool(1,:,1),Params.age,SimLifeCycleProfiles.HighSchool(1,:,1),Params.age,SimLifeCycleProfiles.College(1,:,1))
% % Figure: Earnings
% figure(2)
% plot(Params.age,SimLifeCycleProfiles.NoHighSchool(2,:,1),Params.age,SimLifeCycleProfiles.HighSchool(2,:,1),Params.age,SimLifeCycleProfiles.College(2,:,1))


% Simulate Panel Data
% Same variables as we used for the life-cycle profiles.
SimPanelValues=SimPanelValues_FHorz_Case1_PType(InitialDist,PTypeDistParamNames,Policy, FnsToEvaluate,Params,0,n_a,n_z,N_j,Names_i,0,a_grid,z_grid,pi_z, simoptions);


%% Table 1: Asset-Income ratio and Savings Rate (aggregate and conditional on fixed-type)
% Taking into account population growth n and (conditional) survival rates sj.
ageweights=cumprod(Params.sj).*cumprod((1/(1+Params.n))*ones(1,Params.J)); % First part is based on survival, second part on population growth
ageweights=ageweights./sum(ageweights); % Normalize to sum to 1

% At first, not clear if Table 1 reports the mean of the ratio, or the ratio of the mean (analagously for savings rate)
% From Appendix B it is clear that Table 1 reports the ratio of the mean.

% Following is based on Appendix B: Constructing aggregate consumption, earnings, and assets.

% First, calculate sample means for consumption, earnings, and assets all conditional on age.
% Note that these are what we already calculated for the life-cycle profiles, so rather 
% than get them from panel data simulation we can just get them from there. 
% [The life-cycle profiles command is anyway internally based on a simulated panel data].
AssetsByAgeAndFixedType=[SimLifeCycleProfiles.NoHighSchool.Assets.Mean;SimLifeCycleProfiles.HighSchool.Assets.Mean;SimLifeCycleProfiles.College.Assets.Mean];
EarningsByAgeAndFixedType=[SimLifeCycleProfiles.NoHighSchool.Earnings.Mean;SimLifeCycleProfiles.HighSchool.Earnings.Mean;SimLifeCycleProfiles.College.Earnings.Mean];
IncomeByAgeAndFixedType=[SimLifeCycleProfiles.NoHighSchool.Income.Mean;SimLifeCycleProfiles.HighSchool.Income.Mean;SimLifeCycleProfiles.College.Income.Mean];
ConsumptionByAgeAndFixedType=[SimLifeCycleProfiles.NoHighSchool.Consumption.Mean;SimLifeCycleProfiles.HighSchool.Consumption.Mean;SimLifeCycleProfiles.College.Consumption.Mean];
TransfersByAgeAndFixedType=[SimLifeCycleProfiles.NoHighSchool.TR.Mean;SimLifeCycleProfiles.HighSchool.TR.Mean;SimLifeCycleProfiles.College.TR.Mean];
GrossSavingsByAgeAndFixedType=[SimLifeCycleProfiles.NoHighSchool.GrossSavings.Mean; SimLifeCycleProfiles.HighSchool.GrossSavings.Mean; SimLifeCycleProfiles.College.GrossSavings.Mean];

AssetsByFixedType=sum(AssetsByAgeAndFixedType.*ageweights,2);
Assets_Agg=sum(PTypeDist.*AssetsByFixedType);
EarningsByFixedType=sum(EarningsByAgeAndFixedType.*ageweights,2);
Earnings_Agg=sum(PTypeDist.*EarningsByFixedType);
IncomeByFixedType=sum(IncomeByAgeAndFixedType.*ageweights,2);
Income_Agg=sum(PTypeDist.*IncomeByFixedType);
ConsumptionByFixedType=sum(ConsumptionByAgeAndFixedType.*ageweights,2);
Consumption_Agg=sum(PTypeDist.*ConsumptionByFixedType);
TransfersByFixedType=sum(TransfersByAgeAndFixedType.*ageweights,2);
Transfers_Agg=sum(PTypeDist.*TransfersByFixedType);
GrossSavingsByFixedType=sum(GrossSavingsByAgeAndFixedType.*ageweights,2);
GrossSavings_Agg=sum(PTypeDist.*GrossSavingsByFixedType);

DisposableIncomeByFixedType=IncomeByFixedType+TransfersByFixedType;
DisposableIncome_Agg=Income_Agg+Transfers_Agg;

AssetEarningsRatioByFixedType=AssetsByFixedType./EarningsByFixedType;
AssetEarningsRatio_Agg=Assets_Agg/Earnings_Agg;
AssetIncomeRatioByFixedType=AssetsByFixedType./IncomeByFixedType;
AssetIncomeRatio_Agg=Assets_Agg/Income_Agg;
AssetDispIncomeRatioByFixedType=AssetsByFixedType./DisposableIncomeByFixedType;
AssetDispIncomeRatio_Agg=Assets_Agg/DisposableIncome_Agg;


% SavingsRateByFixedType=GrossSavingsByFixedType./EarningsByFixedType
% SavingsRate_Agg=GrossSavings_Agg/Earnings_Agg
SavingsRateByFixedType=GrossSavingsByFixedType./IncomeByFixedType;
SavingsRate_Agg=GrossSavings_Agg/Income_Agg;
% Income here is after tax but before transfers.

Table1row=[Params.delta, Params.gamma, AssetIncomeRatioByFixedType', AssetIncomeRatio_Agg, SavingsRateByFixedType', SavingsRate_Agg];
% Note: VFI Toolkit could alternatively do this based on stationary agent
% distribution of this model. But will here follow the simulated panel data 
% approach used by Hubbard, Skinner & Zeldes (1994).


%% Table 2: creates the numbers relevant for Table 2
% Conditional probabilities of having consumption approximately equal to income.

HighCorrIncomeCons=zeros(6,3);
NumberOfHHs=zeros(6,3);
LowAssets_NumberOfHHs=zeros(6,3);
LowAssets_HighCorrIncomeCons=zeros(6,3);

% AverageIncomeByAgeBin=[mean(SimLifeCycleProfiles(4,1:9,1)),mean(SimLifeCycleProfiles(4,10:19,1)),mean(SimLifeCycleProfiles(4,20:29,1)),mean(SimLifeCycleProfiles(4,30:39,1)),mean(SimLifeCycleProfiles(4,40:49,1)),mean(SimLifeCycleProfiles(4,50:end,1))];

for ii=1:simoptions.numbersims
    for jj=1:80
        cons=SimPanelValues.Consumption(jj,ii);
        income=SimPanelValues.Income(jj,ii);
        assets=SimPanelValues.Assets(jj,ii);
        age=SimPanelValues.age(jj,ii);
        hhtype=SimPanelValues.ptype(jj,ii);
        
        if age<=29
            agebin=1;
        elseif age<=39
            agebin=2;
        elseif age<=49
            agebin=3;
        elseif age<=59
            agebin=4;
        elseif age<=69
            agebin=5;
        else
            agebin=0;
        end
        
        % Do the needed calculations for each agebin
        if agebin>0
            NumberOfHHs(agebin,hhtype)=NumberOfHHs(agebin,hhtype)+1;
            if (cons/income)>=0.95 && (cons/income)<=1.05
                HighCorrIncomeCons(agebin,hhtype)=HighCorrIncomeCons(agebin,hhtype)+1;
            end
            
            if assets<income % Is not quite clear from paper which HHs in simulated data are considered to satisfy the 'Initial Assets <0.5*Average Income' criterion.
                LowAssets_NumberOfHHs(agebin,hhtype)=LowAssets_NumberOfHHs(agebin,hhtype)+1;
                if (cons/income)>=0.95 && (cons/income)<=1.05
                    LowAssets_HighCorrIncomeCons(agebin,hhtype)=LowAssets_HighCorrIncomeCons(agebin,hhtype)+1;
                end
            end
        end
        % Also do the Total in the sixth row
        NumberOfHHs(6,hhtype)=NumberOfHHs(6,hhtype)+1;
        if (cons/income)>=0.95 && (cons/income)<=1.05
            HighCorrIncomeCons(6,hhtype)=HighCorrIncomeCons(6,hhtype)+1;
        end
        if assets<income % Is not quite clear from paper which HHs in simulated data are considered to satisfy the 'Initial Assets <0.5*Average Income' criterion.
                LowAssets_NumberOfHHs(6,hhtype)=LowAssets_NumberOfHHs(6,hhtype)+1;
                if (cons/income)>=0.95 && (cons/income)<=1.05
                    LowAssets_HighCorrIncomeCons(6,hhtype)=LowAssets_HighCorrIncomeCons(6,hhtype)+1;
                end
            end
    end
end

Table2=[HighCorrIncomeCons./NumberOfHHs;LowAssets_HighCorrIncomeCons./LowAssets_NumberOfHHs];

%% Table 3 
% Campbell-Mankiw-Lusardi Euler Eqns
% My impression from paper is that these regressions from (pooled) simulated panel
% do not correct for the 1% population growth that they assume (ie. that it
% is ignored here). I therefore follow this.
% According to Hubbard, Skinner & Zeldes (1994) the original Cambell-Mankiw
% regressions were on aggregate consumption and income, but in this model
% both of those are constants so this regression would not be possible.
% Presumably it has therefore been done (simulated) microdata instead.

% cons=SimPanelValues(5,:,:);
% income=SimPanelValues(4,:,:);

vecsize=[(N_j-1)*simoptions.numbersims,1];

DeltaC=reshape(SimPanelValues.Consumption(2:end,:)-SimPanelValues.Consumption(1:(end-1),:),vecsize);
DeltaY=reshape(SimPanelValues.Income(2:end,:)-SimPanelValues.Income(1:(end-1),:),vecsize);
DeltalnC=reshape(log(SimPanelValues.Consumption(2:end,:))-log(SimPanelValues.Consumption(1:(end-1),:)),vecsize);
DeltalnY=reshape(log(SimPanelValues.Income(2:end,:))-log(SimPanelValues.Income(1:(end-1),:)),vecsize);
age=reshape(SimPanelValues.age(2:end,:),vecsize);
constant=ones(vecsize);

% All the regressions appear to be two-stage least squared, using lags as
% instruments for the change in income terms.

% First regression
y2=DeltaY(4:end);
X2=[constant(4:end),DeltaC(2:end-2),DeltaC(1:end-3),DeltaY(2:end-2),DeltaY(1:end-3),age(4:end),age(4:end).^2];
[b2,bint,r,rint,stats]=regress(y2,X2);
DeltaYfitted=X2*b2;
y=DeltaC(4:end);
X=[constant(4:end),DeltaYfitted];
[bcolumn1,bintcolumn1,r,rint,stats]=regress(y,X);
% Second regression
y2=DeltalnY(4:end);
X2=[constant(4:end),DeltalnC(2:end-2),DeltalnC(1:end-3),DeltalnY(2:end-2),DeltalnY(1:end-3),age(4:end),age(4:end).^2];
[b2,bint,r,rint,stats]=regress(y2,X2);
DeltalnYfitted=X2*b2;
y=DeltalnC(4:end);
X=[constant(4:end),DeltalnYfitted];
[bcolumn2,bintcolumn2,r,rint,stats]=regress(y,X);
% Third regression
y2=DeltalnY(4:end);
X2=[constant(4:end),DeltalnC(3:end-1),DeltalnC(2:end-2),DeltalnC(1:end-3),DeltalnY(3:end-1),DeltalnY(2:end-2),DeltalnY(1:end-3),age(4:end),age(4:end).^2];
[b2,bint,r,rint,stats]=regress(y2,X2);
DeltalnYfitted=X2*b2;
y=DeltalnC(4:end);
X=[constant(4:end),DeltalnYfitted];
[bcolumn3,bintcolumn3,r,rint,stats]=regress(y,X);
% Fourth regression
y2=DeltalnY(4:end);
X2=[constant(4:end),DeltalnC(3:end-1),DeltalnC(2:end-2),DeltalnC(1:end-3),DeltalnY(3:end-1),DeltalnY(2:end-2),DeltalnY(1:end-3),age(4:end),age(4:end).^2];
[b2,bint,r,rint,stats]=regress(y2,X2);
DeltalnYfitted=X2*b2;
y=DeltalnC(4:end);
X=[constant(4:end),DeltalnYfitted,age(4:end),age(4:end).^2];
[bcolumn4,bintcolumn4,r,rint,stats]=regress(y,X);

% These regressions ignore the surivival probabilies. I suspect that Hubbard, Skinner & Zeldes (1994) do not ignore these.

Table3=nan(8,4);
Table3(1,1)=bcolumn1(2);
Table3(2,1)=bcolumn1(2)/((bintcolumn1(2,2)-bcolumn1(2))/1.96); % Coefficient estimate divided by standard error. Since matlab just gives the 95% confidence intervals (which correspond to plus and minus 1.96 std errors) we have to calculate the standard error in the coefficient.
Table3(3,2)=bcolumn2(2);
Table3(4,2)=bcolumn2(2)/((bintcolumn2(2,2)-bcolumn2(2))/1.96);
Table3(3,3)=bcolumn3(2);
Table3(4,3)=bcolumn3(2)/((bintcolumn3(2,2)-bcolumn3(2))/1.96);
Table3(3,4)=bcolumn4(2);
Table3(4,4)=bcolumn4(2)/((bintcolumn4(2,2)-bcolumn4(2))/1.96);
Table3(5,4)=bcolumn4(3);
Table3(6,4)=bcolumn4(3)/((bintcolumn4(3,2)-bcolumn4(3))/1.96);
Table3(7,4)=bcolumn4(4);
Table3(8,4)=bcolumn4(4)/((bintcolumn4(4,2)-bcolumn4(4))/1.96);


%% Generate the Figures, these are related to those reported in HSZ1994

% Figure 1
figure(1)
plot(Params.age(1:(end-19)), LifeCycProfiles.NoHighSchool.Assets.Mean(1:(end-19))/1000)
hold on
plot(Params.age(1:(end-19)), LifeCycProfiles.HighSchool.Assets.Mean(1:(end-19))/1000)
plot(Params.age(1:(end-19)), LifeCycProfiles.College.Assets.Mean(1:(end-19))/1000)
hold off
title({'Average Assets by Age';'All Certain, $1 Floor'})
legend('No High School Degree', 'High School Degree', 'College')
xlabel('Age')
ylabel('Thousands')

% Figure 3a
figure(5)
plot(Params.age, LifeCycProfiles.NoHighSchool.Earnings.Mean/1000) % Earnings
hold on
plot(Params.age, LifeCycProfiles.NoHighSchool.Consumption.Mean/1000) % Consumption
hold off
title({'Average Consumption and Earnings by Age';'No High School Degree'})
legend('Earnings','Consumption')
xlabel('Age')
ylabel('Thousands')

