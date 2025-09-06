% Cocco, Gomes & Maenhout (2005) - Consumption and Portfolio Choice over the Life Cycle
% This code solves a slightly modified version of this paper (persistent
% shocks instead of permanent shocks, renormalization of the EZ preferences)
%
% Househould chooses a combination of safe and risky asset to save for retirement.
% With exogenous labor risky income can be thought of as a fixed asset.
%
% Labor income risk decreases as people age, which suggests they should switch into 
% riskier assets as they age (their labor income becomes more like a safe asset)
%
% The v shock is markov, so VFI Toolkit calls it an 'z' shock
% The e shock is iid, so VFI Toolkit calls it an 'e' shock
% The eta shock is iid and occurs 'between' periods, so VFI Toolkit calls it an 'u' shock
% To use 'u' shocks we need what VFI Toolkit calls 'riskyasset' (in a standard endogenous state we can choose aprime directly, 
% in riskyasset we choose d, and then some iid shock that occurs between periods determines aprime)
%
% This is not the actual model of CGM2005 as there is a minor difference since they use a permanent (unit
% root) shock, whereas I use a markov (a discretized AR(1) process).
% (Actually I ended up using an autocorrelation of 1, so it should be CGM2005 anyway)
%
% Note that using permanent shock v could be treated in a simpler manner (you
% can essentially divide everything by v, see their paper). I don't do this
% as since the paper was published the empirical evidence has become clear
% that it is better modelled as highly persistent, in which case we would
% need the shock to be treated as I do here. Makes it easier to modify :)
%
% I cheat in one place. In retirement CGM2005 freeze the final period
% transitiory iid shock to it's value in the final working period. Because
% this would require me to treat it as a markov, I instead just freeze it
% as zero valued.
%
% I'm guessing that CGM2005 treat the conditional surivival probability in their EZ
% specification as a discount factor? It does not appear to be discussed in article.

n_d=[51,250]; % Fraction of assets invested in risky asset; Savings of total assets.
n_a=250; % Number of grid points for safe asset and risky asset holdings respectively

% Grid sizes for the labor income shocks
n_z=21; % Persistent shocks: AR(1)
n_e=3; % Transitory shocks: i.i.d.
% Grid sizes for the stochastic component of the excess returns to the risky asset
n_u=5; % i.i.d

vfoptions.riskyasset=1; % riskyasset aprime(d,u)
simoptions.riskyasset=1;

vfoptions.exoticpreferences='EpsteinZin'; % Use Epstein-Zin preferences
vfoptions.EZutils=0; % EZ in consumption-units
vfoptions.EZriskaversion='gamma'; % risk-aversion parameter for EZ
vfoptions.EZeis='psi'; % elasticity of intertemporal substitution for EZ

Names_i={'College','HighSchool','NoHighSchool'};

Params.agejshifter.College=21; % College are 'born' at age 22, so this shifts 'agej' to 'age in years'
Params.agejshifter.HighSchool=19; % Non-college are 'born' at age 20, so this shifts 'agej' to 'age in years'
Params.agejshifter.NoHighSchool=19; % Non-college are 'born' at age 20, so this shifts 'agej' to 'age in years'

Params.J.College=100-Params.agejshifter.College; % Age ends at 100
Params.J.HighSchool=100-Params.agejshifter.HighSchool; % Age ends at 100
Params.J.NoHighSchool=100-Params.agejshifter.NoHighSchool; % Age ends at 100

N_j.College=Params.J.College;
N_j.HighSchool=Params.J.HighSchool;
N_j.NoHighSchool=Params.J.NoHighSchool;

%% To speed up the use of riskyasset we use 'refine_d', which requires us to set the decision variables in a specific order
vfoptions.refine_d=[0,1,1]; % tell the code how many d1, d2, and d3 there are
% Idea is to distinguish three categories of decision variable:
%  d1: decision is in the ReturnFn but not in aprimeFn
%  d2: decision is in the aprimeFn but not in ReturnFn
%  d3: decision is in both ReturnFn and in aprimeFn
% Note: ReturnFn must use inputs (d1,d3,..) 
%       aprimeFn must use inputs (d2,d3,..)
% n_d must be set up as n_d=[n_d1, n_d2, n_d3]
% d_grid must be set up as d_grid=[d1_grid; d2_grid; d3_grid];
% It is possible to solve models without any d1, as is the case here.
simoptions.refine_d=vfoptions.refine_d;

%% Parameters

% Preferences
Params.gamma=10; % CES utility (risk aversion; this number is insanely large; most estimates are more like 1.5 to 4)
% CGM2005 need such a huge risk aversion because the income process is not
% very risky and so (present value of lifetime-earnings) acts as if it were an almost riskless asset (and so agents want
% to hold little further actual risk free assets and instead hold mostly risky assets)
Params.psi=0.5; % Intertemporal Elasticity of substitution

Params.beta=0.96;

Params.Jr.College=65-Params.agejshifter.College+1; % Retire at at 65yrs old
Params.Jr.HighSchool=65-Params.agejshifter.HighSchool+1; % Retire at at 65yrs old
Params.Jr.NoHighSchool=65-Params.agejshifter.NoHighSchool+1; % Retire at at 65yrs old
Params.agej.College=1:1:Params.J.College;
Params.agej.HighSchool=1:1:Params.J.HighSchool;
Params.agej.NoHighSchool=1:1:Params.J.NoHighSchool;

%% Survival probabilities (CGM2005 call this pj, I use nomenclature sj, the survival hazard rate conditional on j)
% CGM2005 take them from the 'mortality tables of the National Center for Health Statistics'
% When I searched for this it was not obvious where to find a 'national'
% rate, there where large numbers of rates but at more disaggregate levels.
% I thus instead just use them from 
% F. C. Bell and M. L. Miller (2005), Life Tables for the United States Social Security Area 1900-2100, Actuarial Study No. 120, Office of the Chief Actuary
% http://www.socialsecurity.gov/oact/NOTES/s2000s.html
% Table 8 — Period Probabilities of Death Within One Year (qx) at Selected Exact Ages, by Sex and Calendar Year (Cont.)
% The raw data from there is
%          Sex and Exact Age
%    |  Male                                                Female
% Year| [0 30 60 65 70 100]                                  [0 30 60 65 70 100]
% 2010| [0.00587,0.00116,0.01086,0.01753,0.02785,0.39134]    [0.00495,0.00060,0.00734,0.01201,0.01912,0.34031]
% I just take the numbers for Males, and then set my actual values based on a linear interpolation of the data.
dj_temp=interp1([0,30,60,65,70,100],[0.00587,0.00116,0.01086,0.01753,0.02785,0.39134],0:1:100,'linear');

Params.sj.College=1-dj_temp((Params.agej.College(1)+Params.agejshifter.College):100);
Params.sj.College(1)=1;
Params.sj.College(end)=0;
Params.sj.HighSchool=1-dj_temp((Params.agej.HighSchool(1)+Params.agejshifter.HighSchool):100);
Params.sj.HighSchool(1)=1;
Params.sj.HighSchool(end)=0;
Params.sj.NoHighSchool=1-dj_temp((Params.agej.NoHighSchool(1)+Params.agejshifter.NoHighSchool):100);
Params.sj.NoHighSchool(1)=1;
Params.sj.NoHighSchool(end)=0;
% I use linear interpolation to fill in the numbers inbetween those reported by Bell & Miller (2005).
% I have additionally imposed that the prob of death at age 20 (or 22 if college) be zero and that prob of death at age 100 is one.

%% Labor income is exogenous
% it consists of 
% i) a deterministic function of age
% ii) a persisitent stochastic component v, modelled as AR(1)
% iii) a transitory stochastic component u, modelled as iid

% The deterministic function of the age is one of the following three
Params.PolynomialCoeffs.College=[-4.3148, 0.3194, -0.0577, 0.0033];
Params.PolynomialCoeffs.HighSchool=[-2.1700, 0.1682, -0.0323, 0.0020];
Params.PolynomialCoeffs.NoHighSchool=[-2.1361, 0.1684, -0.0353, 0.0023];
% Use this
Params.kappa_ij.College=Params.PolynomialCoeffs.College(1)+Params.PolynomialCoeffs.College(2)*(Params.agej.College+Params.agejshifter.College)+Params.PolynomialCoeffs.College(3)*((Params.agej.College+Params.agejshifter.College).^2)/10+Params.PolynomialCoeffs.College(4)*((Params.agej.College+Params.agejshifter.College).^3)/100;
Params.kappa_ij.College(Params.Jr.College:end)=zeros(Params.J.College-Params.Jr.College+1,1); % Make these zero, they will be filled in with actual retirement income based on lambda later in the code
Params.kappa_ij.HighSchool=Params.PolynomialCoeffs.HighSchool(1)+Params.PolynomialCoeffs.HighSchool(2)*(Params.agej.HighSchool+Params.agejshifter.HighSchool)+Params.PolynomialCoeffs.HighSchool(3)*((Params.agej.HighSchool+Params.agejshifter.HighSchool).^2)/10+Params.PolynomialCoeffs.HighSchool(4)*((Params.agej.HighSchool+Params.agejshifter.HighSchool).^3)/100;
Params.kappa_ij.HighSchool(Params.Jr.HighSchool:end)=zeros(Params.J.HighSchool-Params.Jr.HighSchool+1,1); % Make these zero, they will be filled in with actual retirement income based on lambda later in the code
Params.kappa_ij.NoHighSchool=Params.PolynomialCoeffs.NoHighSchool(1)+Params.PolynomialCoeffs.NoHighSchool(2)*(Params.agej.NoHighSchool+Params.agejshifter.NoHighSchool)+Params.PolynomialCoeffs.NoHighSchool(3)*((Params.agej.NoHighSchool+Params.agejshifter.NoHighSchool).^2)/10+Params.PolynomialCoeffs.NoHighSchool(4)*((Params.agej.NoHighSchool+Params.agejshifter.NoHighSchool).^3)/100;
Params.kappa_ij.NoHighSchool(Params.Jr.NoHighSchool:end)=zeros(Params.J.NoHighSchool-Params.Jr.NoHighSchool+1,1); % Make these zero, they will be filled in with actual retirement income based on lambda later in the code

% My kappa_ij is the log of what CGM2005 call fj, hence we need to take logs now
Params.kappa_ij.College(1:Params.Jr.College-1)=log(Params.kappa_ij.College(1:Params.Jr.College-1)); % Note that retirement is zero so avoid taking log of it, will be filled in later
Params.kappa_ij.HighSchool(1:Params.Jr.HighSchool-1)=log(Params.kappa_ij.HighSchool(1:Params.Jr.HighSchool-1)); % Note that retirement is zero so avoid taking log of it, will be filled in later
Params.kappa_ij.NoHighSchool(1:Params.Jr.NoHighSchool-1)=log(Params.kappa_ij.NoHighSchool(1:Params.Jr.NoHighSchool-1)); % Note that retirement is zero so avoid taking log of it, will be filled in later

% The replacement rate for pensions is lambda (because of how codes work, we want the logs of lambda)
Params.lambda.College=log(0.938873);
Params.lambda.HighSchool=log(0.68212);
Params.lambda.NoHighSchool=log(0.88983);
% And kappa_ij in retirement is just a frozen version of kappa_ij, but reduced by lambda (just add lambda to kappa_ij (both are in logs))
Params.kappa_ij.College(Params.Jr.College:end)=Params.kappa_ij.College(Params.Jr.College-1)+Params.lambda.College;
Params.kappa_ij.HighSchool(Params.Jr.HighSchool:end)=Params.kappa_ij.HighSchool(Params.Jr.HighSchool-1)+Params.lambda.HighSchool;
Params.kappa_ij.NoHighSchool(Params.Jr.NoHighSchool:end)=Params.kappa_ij.NoHighSchool(Params.Jr.NoHighSchool-1)+Params.lambda.NoHighSchool;



% I discretize this as an AR(1) with rho_z=1 and using the sigma_epsilonz
% that CGM2005 estimated. Note that really for permanent shocks we should
% renormalize the model, what is being solved here is harder
% computationally but means you can easily switch to a different earnings
% process.
Params.rho_z=1;
Params.sigma_epsilonz.College=sqrt(0.0169);
Params.sigma_epsilonz.HighSchool=sqrt(0.0106);
Params.sigma_epsilonz.NoHighSchool=sqrt(0.0105);
% Following comment uses notation of CGM2005 (epsilonz is what they call u)
% % Note that the sum of two normally distributed variables is itself normal,
% % so the only way to get from sigma_u to sigma_w and sigma_zeta is
% % arbitrary (it could be non-arbitrary if you have the original panel data,
% % but as all we have is sigma_u it must be arbitrary)
% % https://en.wikipedia.org/wiki/Sum_of_normally_distributed_random_variables
kirkbyoptions.nSigmas=2;
kirkbyoptions.initialj1mewz=0; % WHAT DO CGM2005 USE TO START?
[z_grid_J.College, pi_z_J.College,jequaloneDistP,otheroutputs] = discretizeLifeCycleAR1_KFTT(zeros(1,N_j.College),Params.rho_z*ones(1,N_j.College),Params.sigma_epsilonz.College*ones(1,N_j.College),n_z,N_j.College,kirkbyoptions);
[z_grid_J.HighSchool, pi_z_J.HighSchool,jequaloneDistP,otheroutputs] = discretizeLifeCycleAR1_KFTT(zeros(1,N_j.HighSchool),Params.rho_z*ones(1,N_j.HighSchool),Params.sigma_epsilonz.HighSchool*ones(1,N_j.HighSchool),n_z,N_j.HighSchool,kirkbyoptions);
[z_grid_J.NoHighSchool, pi_z_J.NoHighSchool,jequaloneDistP,otheroutputs] = discretizeLifeCycleAR1_KFTT(zeros(1,N_j.NoHighSchool),Params.rho_z*ones(1,N_j.NoHighSchool),Params.sigma_epsilonz.NoHighSchool*ones(1,N_j.NoHighSchool),n_z,N_j.NoHighSchool,kirkbyoptions);

% Transitory income shock, e
Params.sigma_e.College=sqrt(0.0584);
Params.sigma_e.HighSchool=sqrt(0.0738);
Params.sigma_e.NoHighSchool=sqrt(0.1056);
Params.rho_e=0; % iid
farmertodaoptions.nSigmas=1;
farmertodaoptions.method='even';
[e_grid.College, pi_e.College]=discretizeAR1_FarmerToda(0,Params.rho_e,Params.sigma_e.College,n_e,farmertodaoptions);
[e_grid.HighSchool, pi_e.HighSchool]=discretizeAR1_FarmerToda(0,Params.rho_e,Params.sigma_e.HighSchool,n_e,farmertodaoptions);
[e_grid.NoHighSchool, pi_e.NoHighSchool]=discretizeAR1_FarmerToda(0,Params.rho_e,Params.sigma_e.NoHighSchool,n_e,farmertodaoptions);
pi_e.College=pi_e.College(1,:)'; % This is iid
pi_e.HighSchool=pi_e.HighSchool(1,:)'; % This is iid
pi_e.NoHighSchool=pi_e.NoHighSchool(1,:)'; % This is iid
% To make e equal zero in retirement, I set it up as age-dependent
e_grid_J.College=[e_grid.College.*ones(1,Params.Jr.College-1),zeros(n_e,N_j.College-Params.Jr.College+1)];
e_grid_J.HighSchool=[e_grid.HighSchool.*ones(1,Params.Jr.HighSchool-1),zeros(n_e,N_j.HighSchool-Params.Jr.HighSchool+1)];
e_grid_J.NoHighSchool=[e_grid.NoHighSchool.*ones(1,Params.Jr.NoHighSchool-1),zeros(n_e,N_j.NoHighSchool-Params.Jr.NoHighSchool+1)];
pi_e_J.College=[pi_e.College.*ones(1,Params.Jr.College-1),[0;1;0].*ones(1,N_j.College-Params.Jr.College+1)];
pi_e_J.HighSchool=[pi_e.HighSchool.*ones(1,Params.Jr.HighSchool-1),[0;1;0].*ones(1,N_j.HighSchool-Params.Jr.HighSchool+1)];
pi_e_J.NoHighSchool=[pi_e.NoHighSchool.*ones(1,Params.Jr.NoHighSchool-1),[0;1;0].*ones(1,N_j.NoHighSchool-Params.Jr.NoHighSchool+1)];

%% Retirement income
% Is a constant fraction of last-working period income
% This allows a simple computational trick for a model with permanent
% income shocks in which you are doing the renormalization
% Here we just have to force it.
% Freeze z at the retirement value, and then just apply a scaling factor to reduce it
for jj=Params.Jr.College:N_j.College
    z_grid_J.College(:,jj)=z_grid_J.College(:,Params.Jr.College-1);
    pi_z_J.College(:,:,jj)=eye(n_z,n_z);
end
for jj=Params.Jr.HighSchool:N_j.HighSchool
    z_grid_J.HighSchool(:,jj)=z_grid_J.HighSchool(:,Params.Jr.HighSchool-1);
    pi_z_J.HighSchool(:,:,jj)=eye(n_z,n_z);
end
for jj=Params.Jr.NoHighSchool:N_j.NoHighSchool
    z_grid_J.NoHighSchool(:,jj)=z_grid_J.NoHighSchool(:,Params.Jr.NoHighSchool-1);
    pi_z_J.NoHighSchool(:,:,jj)=eye(n_z,n_z);
end
% Eqn (5) of CGM2005 on pg 496 seems to say they also freeze v_Jr. I cannot
% really do this as it is a pain in the ass (means I would have to treate v
% as a markov). Instead I just overwrite it as zero.
% The replacement rate is lambda, and this done via the return fn

%% Share of assets invested in the risky asset
riskyshare_grid=linspace(0,1,n_d(1))'; % Share of assets, from 0 to 1

%% Grids for the total assets
a_grid=300*linspace(0,1,n_a)'; % This grid has to be equally spaced (is hardcoded into Phiaprime)

% Set up d for VFI Toolkit (is the two decision variables)
d_grid=[riskyshare_grid; a_grid]; % Note: this does not have to be a_grid, I just chose to use same grid for savings as for assets

%% Asset returns
% Rather than treat the excess mean and stardard deviation of the risky asset
% returns as two separate concepts, I just stick them together as
% properties of the u shocks.

Params.Rf=1.02; % (gross) return to the risk free asset Rf=1+rf

% Excess returns to the risky asset consist of a deterministic constant and a shock
% Rr-Rf=mew+eta
% I reformulate this to just be
% Rr-Rf=u, u~N(mew,sigma_u^2)

% u is the stochastic component of the excess returns to the risky asset
Params.mew=0.04; % Mean of the excess returns to the risky asset
% Note: Table 4 reports mew as 0.06, (mew-1 in their notation), but the
% text on pg 500 makes clear it should be 0.04 (which also fits the
% description of the number in Table 4, just not the formula/expression)
Params.sigma_u=0.157; % Standard deviation of innovations to the risky asset
Params.rho_u=0; % Asset return risk component is modeled as iid
[u_grid, pi_u]=discretizeAR1_FarmerToda(Params.mew,Params.rho_u,Params.sigma_u,n_u);
pi_u=pi_u(1,:)'; % This is iid

%% Put in form for VFI Toolkit
% z are markov
% e are iid
% u are iid and occur between time periods

% z and u are already set up appropriately. Just need to put e into the vfoptions and simoptions.
vfoptions.n_e=n_e;
vfoptions.e_grid=e_grid_J;
vfoptions.pi_e=pi_e_J;
simoptions.n_e=vfoptions.n_e;
simoptions.e_grid=vfoptions.e_grid;
simoptions.pi_e=vfoptions.pi_e;

%% Discount factor and return function
DiscountFactorParamNames={'beta','sj'};

ReturnFn=@(savings,a,z,e,kappa_ij)...
    CoccoGomesMaenhout2005_ReturnFn(savings,a,z,e,kappa_ij);

%% Define aprime function used for riskyasset (value of next period assets, determined by this period decision, and u shock)

% riskyasset: aprime_val=aprimeFn(d,u)
aprimeFn=@(riskyshare,savings,u, Rf) CoccoGomesMaenhout2005_aprimeFn(riskyshare,savings, u, Rf); % Will return the value of aprime
% Note that u is risky asset return and effectively includes both the (excess) mean and standard deviation of risky assets

% Put the risky asset into vfoptions and simoptions
vfoptions.aprimeFn=aprimeFn;
vfoptions.n_u=n_u;
vfoptions.u_grid=u_grid;
vfoptions.pi_u=pi_u;
simoptions.aprimeFn=aprimeFn;
simoptions.n_u=n_u;
simoptions.u_grid=u_grid;
simoptions.pi_u=pi_u;
% Because a_grid and d_grid are involved in risky assets, but are not
% normally needed for agent distriubiton simulation, we have to also
% include these in simoptions
simoptions.a_grid=a_grid;
simoptions.d_grid=d_grid;

%% Solve for the value function and policy fn
vfoptions.verbose=1;
tic;
[V, Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i, d_grid, a_grid, z_grid_J, pi_z_J, ReturnFn, Params, DiscountFactorParamNames,vfoptions);
toc

%% Take a look at the policy
fig4=figure(4);
plot(a_grid,a_grid(Policy.HighSchool(2,:,11,2,end))) 
hold on
plot(a_grid,a_grid(Policy.HighSchool(2,:,11,2,20-Params.agejshifter.HighSchool))) % age 20
plot(a_grid,a_grid(Policy.HighSchool(2,:,11,2,30-Params.agejshifter.HighSchool))) % age 30
plot(a_grid,a_grid(Policy.HighSchool(2,:,11,2,40-Params.agejshifter.HighSchool))) % age 40
plot(a_grid,a_grid(Policy.HighSchool(2,:,11,2,50-Params.agejshifter.HighSchool))) % age 50
plot(a_grid,a_grid(Policy.HighSchool(2,:,11,2,60-Params.agejshifter.HighSchool))) % age 60
plot(a_grid,a_grid(Policy.HighSchool(2,:,11,2,70-Params.agejshifter.HighSchool))) % age 70
hold off

%% Make sure that u shocks seem appropriate in terms of where they send you for aprime
% aprimeFnParamNames in same fashion
l_d=length(n_d);
l_u=length(n_u);
temp=getAnonymousFnInputNames(aprimeFn);
if length(temp)>(l_d+l_u)
    aprimeFnParamNames={temp{l_d+l_u+1:end}}; % the first inputs will always be (d,u)
else
    aprimeFnParamNames={};
end
aprimeFnParamsVec=CreateVectorFromParams(Params, aprimeFnParamNames,N_j);
[aprimeIndex,aprimeProbs]=CreateaprimeFnMatrix_RiskyAsset(aprimeFn, n_d, n_a, n_u, d_grid, a_grid, u_grid, aprimeFnParamsVec,0); % Note, is actually aprime_grid (but a_grid is anyway same for all ages)
% Note: aprimeIndex is [N_d,N_u], whereas aprimeProbs is [N_d,N_u]
% plot the first, mid, and last u shock
PolicyKron=KronPolicyIndexes_FHorz_Case2(Policy.HighSchool, n_d, n_a, n_z,N_j.HighSchool,simoptions.n_e); % Just pretend e is another z
N_d=prod(n_d); N_u=prod(n_u);

fig5=figure(5);
subplot(3,1,1); plot(a_grid,a_grid(aprimeIndex(PolicyKron(:,11,2,20-Params.agejshifter.HighSchool)))) % age 20
title('Next period state for the different u shocks')
hold on
subplot(3,1,1); plot(a_grid,a_grid(aprimeIndex(PolicyKron(:,11,2,20-Params.agejshifter.HighSchool)+N_d*1))) % age 20
subplot(3,1,1); plot(a_grid,a_grid(aprimeIndex(PolicyKron(:,11,2,20-Params.agejshifter.HighSchool)+N_d*2))) % age 20
subplot(3,1,1); plot(a_grid,a_grid(aprimeIndex(PolicyKron(:,11,2,20-Params.agejshifter.HighSchool)+N_d*3))) % age 20
subplot(3,1,1); plot(a_grid,a_grid(aprimeIndex(PolicyKron(:,11,2,20-Params.agejshifter.HighSchool)+N_d*(N_u-1)))) % age 20
hold off
subplot(3,1,2); plot(a_grid,riskyshare_grid(Policy.HighSchool(1,:,11,2,20-Params.agejshifter.HighSchool)))
title('Corresponding risky asset share')
ylim([0,1])
subplot(3,1,3); plot(a_grid,a_grid(Policy.HighSchool(2,:,11,2,20-Params.agejshifter.HighSchool)))
title('Corresponding asset savings')


%% Mass of each PType
% Because CGM2005 do not compute any aggregate (just statistics conditional
% on the PType) they do not actually need to define this
Params.PTypeMasses=[0.3,0.4,0.3];
PTypeDistParamNames={'PTypeMasses'};

%% If we wanted to compute the stationary distribution
jequaloneDist=zeros([n_a,n_z,vfoptions.n_e]);
jequaloneDist(1,floor(n_z/2)+1,floor(n_e/2)+1)=1; % start everyone with zero assets and mid-point for all shocks
Params.AgeWeights.College=ones(N_j.College,1)/N_j.College;
Params.AgeWeights.HighSchool=ones(N_j.HighSchool,1)/N_j.HighSchool;
Params.AgeWeights.NoHighSchool=ones(N_j.NoHighSchool,1)/N_j.NoHighSchool;
AgeWeightsParamNames={'AgeWeights'};

StationaryDist=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightsParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_z,N_j,Names_i,pi_z_J, Params,simoptions);

% Take a look at distribution of agents over assets (check that they are not at top of grid)
fig1=figure(1);
plot(a_grid,StationaryDist.ptweights(1)*cumsum(sum(sum(sum(sum(StationaryDist.College,5),4),3),2)));
hold on
plot(a_grid,StationaryDist.ptweights(2)*cumsum(sum(sum(sum(sum(StationaryDist.HighSchool,5),4),3),2)));
plot(a_grid,StationaryDist.ptweights(3)*cumsum(sum(sum(sum(sum(StationaryDist.NoHighSchool,5),4),3),2)));
hold off
title('Stationary Dist of agents, cdf over assets')

% % Check amount of mass in top 50 grid points
% temp1=StationaryDist.College(end-50:end,:,:,:); sum(temp1(:))
% temp2=StationaryDist.HighSchool(end-50:end,:,:,:); sum(temp2(:))
% temp3=StationaryDist.NoHighSchool(end-50:end,:,:,:); sum(temp3(:))
% % Check amount of mass in top 20 grid points
% temp1=StationaryDist.College(end-20:end,:,:,:); sum(temp1(:))
% temp2=StationaryDist.HighSchool(end-20:end,:,:,:); sum(temp2(:))
% temp3=StationaryDist.NoHighSchool(end-20:end,:,:,:); sum(temp3(:))

%% For computing aggregates and life cycle profiles
FnsToEvaluate.Assets = @(d1,d2,a,v,e) a;
FnsToEvaluate.ShareRisky = @(d1,d2,a,v,e) d1; % Note: this is correct at individual level, but is not the same as the aggregate share of risky assets
FnsToEvaluate.Income = @(d1,d2,a,v,e, kappa_ij) exp(kappa_ij+v+e);
FnsToEvaluate.AssetIncomeRatio = @(d1,d2,a,v,e, kappa_ij) a/exp(kappa_ij+v+e);
FnsToEvaluate.SavingsRate = @(d1,d2,a,v,e, kappa_ij) (d2-a)/exp(kappa_ij+v+e);
FnsToEvaluate.Consumption = @(d1,d2,a,v,e, kappa_ij) CoccoGomesMaenhout2005_Consumption(d1,d2,a,v,e, kappa_ij);
FnsToEvaluate.CashOnHand = @(d1,d2,a,v,e, kappa_ij) a+exp(kappa_ij+v+e); % Essentially, assets plus earnings [This is only correct for working age]

AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1_PType(StationaryDist, Policy, FnsToEvaluate, Params, n_d, n_a, n_z,N_j, Names_i, d_grid, a_grid, z_grid_J, simoptions);

simoptions.lifecyclepercentiles=0; % Just mean and median, no percentiles.
simoptions.agejshifter=Params.agejshifter; % Note, this is 21,19,19, but internally will become 2,0,0
AgeConditionalStats=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy, FnsToEvaluate, Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid_J, simoptions);

% Plot life-cycle profile for mean of Risky Assets as a share of Total Assets
% plot(AgeConditionalStats.ShareRisky.Mean)


%% Figure 2 of CGM2005
% CGM2005 plot everything in terms of 'cash on hand'. While this makes
% sense for some theory, and is easy given the way they solve the problem,
% it is a real pain in the ass to do given how VFI Toolkit operate. So rather 
% than reproduce these I just plot them as in terms of assets.

% For consumption I need to do
ValuesOnGrid=EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1_PType(Policy, FnsToEvaluate, Params,n_d,n_a,n_z,N_j,Names_i,d_grid, a_grid, z_grid_J, simoptions);

% Note, CGM2005 just plot the 'HighSchool' type
% I plot for the median value of each schock
fig2=figure(2);
subplot(3,1,1); plot(a_grid,riskyshare_grid(Policy.HighSchool(1,:,11,2,end))) 
xlabel('assets')
title('Risky share in final period of life of High-schooler')
subplot(3,1,2); plot(a_grid,riskyshare_grid(Policy.HighSchool(1,:,11,2,20-Params.agejshifter.HighSchool))) % age 20
hold on
plot(a_grid,riskyshare_grid(Policy.HighSchool(1,:,11,2,30-Params.agejshifter.HighSchool))) % age 30
plot(a_grid,riskyshare_grid(Policy.HighSchool(1,:,11,2,55-Params.agejshifter.HighSchool))) % age 55
plot(a_grid,riskyshare_grid(Policy.HighSchool(1,:,11,2,75-Params.agejshifter.HighSchool))) % age 75
hold off
xlabel('assets')
title('Risky share of high-schooler')
legend('Age 20','Age 30','Age 55','Age 75')
ylim([0,1])
subplot(3,1,3); plot(a_grid,ValuesOnGrid.Consumption.HighSchool(:,11,2,20-Params.agejshifter.HighSchool)) % age 20
hold on
plot(a_grid,ValuesOnGrid.Consumption.HighSchool(:,11,2,30-Params.agejshifter.HighSchool)) % age 30
plot(a_grid,ValuesOnGrid.Consumption.HighSchool(:,11,2,55-Params.agejshifter.HighSchool)) % age 55
plot(a_grid,ValuesOnGrid.Consumption.HighSchool(:,11,2,75-Params.agejshifter.HighSchool)) % age 75
hold off
xlabel('assets')
title('Consumption of high-schooler')
legend('Age 20','Age 30','Age 55','Age 75')
% I only plot half of the first panel, the other half is a complete-markets result that I do not attempt to reproduce
% I could have ploted ValuesOnGrid.RiskyShare.HighSchool instead of riskyshare_grid(Policy.HighSchool. I just do the alternative to show how it is done

%% Figure 3 of CGM2005
% Note, these are grouped across permanent types
fig3=figure(3);
subplot(3,1,1); plot(Params.agejshifter.HighSchool+(1:1:N_j.HighSchool),AgeConditionalStats.Consumption.Mean)
hold on
plot(Params.agejshifter.HighSchool+(1:1:N_j.HighSchool),AgeConditionalStats.Income.Mean)
plot(Params.agejshifter.HighSchool+(1:1:N_j.HighSchool),AgeConditionalStats.Assets.Mean)
hold off
legend('Consumption','Income','Wealth')
xlabel('Age')
title('Life-cycle profiles (age-conditional mean)')
subplot(3,1,2);
% I don't bother to do this calculation
subplot(3,1,3); plot(Params.agejshifter.HighSchool+(1:1:N_j.HighSchool),AgeConditionalStats.ShareRisky.Mean)
hold on
plot(Params.agejshifter.HighSchool+(1:1:N_j.HighSchool),AgeConditionalStats.ShareRisky.QuantileCutoffs(2,:)) % 5th percentile
plot(Params.agejshifter.HighSchool+(1:1:N_j.HighSchool),AgeConditionalStats.ShareRisky.QuantileCutoffs(20,:)) % 95th percentile
hold off
legend('Mean','5th percentile','95th percentile')
xlabel('Age')
title('Risky share')
ylim([0,1])



%% Figure 9 of CGM2005
% Plot of risky share for the different permanent types
fig9=figure(9);
plot(Params.agejshifter.NoHighSchool+(1:1:N_j.NoHighSchool),AgeConditionalStats.ShareRisky.NoHighSchool.Mean)
hold on
plot(Params.agejshifter.HighSchool+(1:1:N_j.HighSchool),AgeConditionalStats.ShareRisky.HighSchool.Mean)
plot(Params.agejshifter.College+(1:1:N_j.College),AgeConditionalStats.ShareRisky.College.Mean)
hold on
legend('No High School','High School','College')
title('Life-cycle profile of Risky Share')
ylim([0,1])




