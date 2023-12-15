% Attanasio, Low & Sanchez-Marcos (2008) - Explaining Changes in Female Labor Supply in a Life-Cycle Model

% Two notational changes: call shocks to earnings z (ALSM2008 call them
% upsilon, eqns 10 and 11), and I call h_m what ALSM2008 would call exp(h_m) [because that
% way it is analagous to h_f]

% Note: In principle when using unit-root processes you can renormalize the
% whole problem and eliminate them from the state space. This is a better
% idea in sense it is compuationally faster/less demanding but is not done here. 
% The advantage of the approach here is that you can easily switch from unit-root 
% to a more realistic age-depedent VAR(1) process (e.g., autocorrelation
% and shocks could be age-dependent).

% Only Figure 10 is from the paper, the rest are extra things I created to
% help understand the model. The way households make decisions on female
% labor force participation (which is what ALSM2008 paper is all about) is
% all sensible, but the model looks really odd in terms of assets and
% especially the cdf over assets. [Inside return function are some commented out
% lines that I used to check that if I force something with assets and participation
% then model gives what you would expect, and it worked fine.]

% The way the model and parameter values are done is this paper is a mess.
% Required constantly leafing back and forth to figure out what the
% parameter values are. (E.g., see below for h_f and G(a_t), both of which were 
% left unclear by the paper and had to be found in the original codes which are 
% happily easily available.) Please put all your parameter values (including
% general eqm and alternative calibrations) into a table in an appendix :)

% 't' in ALSM2008 is the actual age (not model period; this is the likely interpretation of paper, but wasn't obvious)
% How is y_m0 set? (I guess it is just normalized to 1)
% (Table 3 of ALSM2008 says y_f0/y_m0 equals 0.64, but this is not enough to pin down both of them)
% No mention what initial h_f is (their codes make is look like it is 1)
% Some model outputs depend on the weights assigned to each age, but ALSM2008 never mentioned any weights (seems they equal weight them)

% Eqn (11) of ALSM2008 is incorrect. 
% First, it specifies the Covariance matrix of the innovations to the permanent shocks, but puts the
% 'correlation' in the off-diagonals. This is incorrect, the off-diagonals should be the covariance, 
% which is the correlation divided by the product of the standard deviations.
% Second, it says that zeta is N(Mew,SigmaSq), and has Mew being negative.
% This leads z (which they call upsilon) to drift negative (and so exp(z)
% gradually decreases with age). This is probably (but not certainly) an error.
% Authors appear to have gotten confused with the log() and exp(), as the expression
% for the Mew is what you use in that situation. ALSM2008 use the expression for the 
% mean of a random variable U which is distributed N(Mew,sigma) that will give us that 
% exp(U) is mean zero. But that is not how it gets used here as here it is z which is 
% in logs, not the innovations to z, the zetas, which are in levels. (Maybe
% they just use it as an approximation for the z in logs, but not clear)
% I guess that maybe this is a typo in the paper, but I go with it anyway. 
% The original codes cannot tell us if it is
% incorrect, as they just use a grid on these shocks that comes from having
% discretized the process, but contain nothing about the discretization
% that lead to these grids (as far as I can see).

% This code does not quite produce the same results as ALSM2008. Essentially, ALSM2008
% report that 47% of mothers with children under 3 years old
% participate/work (Table 3), whereas here it is almost zero.
% Most likely explanation is that ALSM2008 is rather vague on model and calibration, so
% I think these codes might get parametrization slightly off (potentially
% it is numerical error, but issues around parametrization seem most likely cause).
% E.g Look at Figure 10, both this code and ALSM2008 have same shapes over age, but
% the fraction of mothers working is way lower here than in paper.

n_d=2; % Female labor force participation
n_a=[501,151]; % Number of grid points for asset, and female labor force history

% Grid sizes for the labor productivity shocks
n_z=[11,11]; % Permament component (unit root), for female and male respectively (has to be this order because of how the innovations to them are defined by ALSM2008)

N_j=50;
Params.agejshifter=22; % So initial age is 23

Names_i={'nochild','youngmother','oldermother'};
% There are three types of household, the first never has children, the second which 
% have children when young and the third which have children when old
% [They differ in the parameters 'e' and 'childcarecosts']

% A line I needed for running on the Server
addpath(genpath('./MatlabToolkits/'))

CreateFigures=0


%% Parameters
Params.J=N_j;

Params.Jr=41; % 40 working years (Table 2 of ALSM2008)

Params.agej=1:1:Params.J; % the model period
Params.age=Params.agejshifter+Params.agej; % the age

% Preferences
Params.beta=0.98; % time discount factor
Params.gamma=1.5; % risk aversion (called nu in ALSM2008 codes)
Params.phi1=0.038; % scales down utility of consumption when female works
Params.phi2=0.002; % disutility of female working

% Asset returns
Params.R=1.015;

% Earnings
% Male earnings are a deterministic polynomial
Params.alpha1=0.034;
Params.alpha2=-0.00033;
% Table 2 of ALMS2008 says theses are parameters for 'regression of log wage on age and age^2
Params.h_m=Params.alpha1*(Params.age)+Params.alpha2*(Params.age.^2); % Note: is a polynomial in age [eqn 9, pg 1528 of ALSM2008]
% Because I use y_m=y_m0*h_m*z_m in return function, eqn (9) makes clear I
% need h_m to be exponential of this value
Params.h_m=exp(Params.h_m);
% ALSM2008 report that y_f0/y_m0=0.64 (Offered gender wage gap, Table 3 of ALSM2008)
% I am guessing they just normalize y_m0 to 1????
Params.y_m0=1;
Params.y_f0=0.64*Params.y_m0;
% Female human capital is set up later
% ALSM2008 codes allow for 'aggrowth' of wages, but this is set to zero.

% Exogenous processes is a joint-random walk
% z_f(j)=z_f(j-1)+zeta_f
% z_m(j)=z_m(j-1)+zeta_m
% Rewrite as a VAR(1)
% z_f(j)=mu_zeta_f + rho_z_f*z_f(j-1)+0**z_m(j-1)+zeta_f
% z_m(j)=mu_zeta_m + 0*z_f(j-1)+rho_z_m*z_m(j-1)+zeta_m
% (zeta_f,zeta_m)~N(0, sigma_zeta^2)
% mu_zeta=[-sigma_zeta_f^2/2; -sigma_zeta_m^2/2]
% sigma_zeta^2=[sigma_zeta_f^2, rho_zeta_f_zeta_m*(sigma_zeta_f*sigma_zeta_m); rho_zeta_f_zeta_m*(sigma_zeta_f*sigma_zeta_m), sigma_zeta_m^2] % Paper contains a typo in this (lists the off-diagonals as correlation when it should be the covariance)
% Note: ALSM2008 have mu_zeta as the mean of zeta (whereas I write it as the constant in a VAR(1))
Params.rho_zeta_f=1;
Params.rho_zeta_m=1;
Params.sigma_zeta_f=0.1342; % Table of paper reports as 0.13, in codes I found it is sqrt(0.018) [actually codes say the variance is 0.18] so I used the more accurate value here [because I was having issues with Var-Covar matrix failing to be positive semi-definite]
Params.sigma_zeta_m=0.1342;
Params.rho_zeta_f_zeta_m=0.25; % Values of these three parameters are from Table 2
Params.mu_zeta=[-Params.sigma_zeta_f^2/2; -Params.sigma_zeta_m^2/2];
% Above line is what ALSM2008 say,but seems likely incorrect (see comment
% at start), I tried out the following but didn't really do anything)
% Params.mu_zeta=[0; 0];

Params.sigma_zeta_sq=[Params.sigma_zeta_f^2, Params.rho_zeta_f_zeta_m*(Params.sigma_zeta_f*Params.sigma_zeta_m); Params.rho_zeta_f_zeta_m*(Params.sigma_zeta_f*Params.sigma_zeta_m), Params.sigma_zeta_m^2];
% The VAR(1) is thus
Params.Mew_z=Params.mu_zeta;
Params.Rho_z=[1,0; 0,1];
Params.SigmaSq_zeta=Params.sigma_zeta_sq;
% Note: This would be a non-stationary VAR(1), but this doesn't matter as
% model is finite-horizon and we use life-cycle VAR(1) discretization
% method.



%% Children (just used to calculate the consumption equivalence scale and the childcare cost, which are two age-dependent parameters)
% Number of children of each age (which is entirely determinsitic; bottom pg 1526, top of 1527)
% The first type of woman never has children
% The second type of woman has her first child at 24, with a second child arriving two years after the first
birthfirstchild_youngmother=24;
% The third type of woman has her first child at 29, with a second child arriving two years after the first
birthfirstchild_oldermother=29;

% Number of children will influence two age-dependent parameters
% equivalence scale for consumption which depends on age and number of children
Params.e.nochild=1.67*ones(1,N_j);
Params.e.youngmother=zeros(1,N_j); % placeholder, filled during for-loop below
Params.e.oldermother=zeros(1,N_j); % placeholder, filled during for-loop below
% childcare costs (if mother works) [F(a_t) in notation of ALSM2008]
Params.childcarecost.nochild=zeros(1,N_j);
Params.childcarecost.youngmother=zeros(1,N_j); % placeholder, filled during for-loop below
Params.childcarecost.oldermother=zeros(1,N_j); % placeholder, filled during for-loop below
Params.p_childcare=12791; % scales the childcare cost (Table 3 of ALSM2008)

% How to calculate e:
% e follows the McClements scale as explained in footnote 4 of ALSM2008: "a
% childless couple is equivalent to 1.67 adults. A couple with one child is
% equivalent to 1.9 adults if the child is under the age of 3, to 2 adults
% if the child is between 3 and 7, 2.07 adults if the child is between 8
% and 12, and 2.2 adults if the child is between 13 and 18. [We] assume
% that each couple has two children who arrive at a predetermined age and
% leave at age 18."
% Note: Is obvious from codes that it is a fixed addition for each child of
% that age. (mpiaer.f90; lines 329-354)

% How to calculate childcarecosts:
% ALSM2008 on pg 1528 say: F(a_t)=pG(a_t), and report a value for p in
% Table 2. They then say "We estimate the function G(a_t) from expenditure
% data of households with children of the relevant ages". Then on page 1531
% they say "We estimate the function G(a_t) directly from the data. In
% particular, for households where the mother is working, we regress total
% child care expenditure on the age of the youngest child, the age of the
% oldest child, the number of children, and a dummy that equals one if the
% age of the youngest child is zero. The shape of G(a_t) can then be
% derived from the coefficients of this regression function, considering
% that in our model eall women who have children have two of them, and at
% the same interval between children of two years. This implies that the
% child care cost  can be expressed as a function of the age of the oldest
% child."
% But the paper never actually goes on to provide the estimated regression
% coefficients (or alternatively the coefficients in G(a_t)).
% In the codes of ALSM2008: mpiaer.f90, line 321-325:
% infants(t) = 0.d0 
% if (t.ge.agekid.AND.t.le.agekid+7) then
%   infants(t) = kidalpha0*chcostfunc(t-agekid+1)/chcostfunc(1)
% endif
% And the chchostfunc is a vector of 8 numbers:
Params.childcostcoeffs=[649.38, 649.38, 713.18, 713.18, 667.67, 622.15, 576.64, 531.12];
% I just do the normalization now
Params.childcostcoeffs=Params.childcostcoeffs./Params.childcostcoeffs(1); 
% kidalpha0 in ALSM2008 codes is what their paper calls p.
% (Actually, in their codes kidalpha0=9000, while paper says p=12791.)
% (ALSM2008 codes also have a kidalpha1=-kidalpha0/18, which is kidalpha1=-450 in imported params. But this doesn't appear to be used for anything anyway.)

for jj=1:Params.J
    % Young mother
    nkids0to2=0;
    nkids3to7=0;
    nkids8to12=0;
    nkids13to18=0;
    agefirstchild=(jj+Params.agejshifter)-birthfirstchild_youngmother;
    agesecondchild=(jj+Params.agejshifter)-(birthfirstchild_youngmother+2); % second child born two year later
    if 0<=agefirstchild && agefirstchild<=2
        nkids0to2=nkids0to2+1;
    elseif 3<=agefirstchild && agefirstchild<=7
        nkids3to7=nkids3to7+1;
    elseif 8<=agefirstchild && agefirstchild<=12
        nkids8to12=nkids8to12+1;
    elseif 13<=agefirstchild && agefirstchild<=18
        nkids13to18=nkids13to18+1;
    end
    if 0<=agesecondchild && agesecondchild<=2
        nkids0to2=nkids0to2+1;
    elseif 3<=agesecondchild && agesecondchild<=7
        nkids3to7=nkids3to7+1;
    elseif 8<=agesecondchild && agesecondchild<=12
        nkids8to12=nkids8to12+1;
    elseif 13<=agesecondchild && agesecondchild<=18
        nkids13to18=nkids13to18+1;
    end
    Params.e.youngmother(jj)=1.67+0.233*nkids0to2+0.33*nkids3to7+0.4*nkids8to12+0.533*nkids13to18;
    if 0<=agefirstchild && agefirstchild<=7 % and childcarecosts
        Params.childcarecost.youngmother(jj)=Params.childcostcoeffs(agefirstchild+1);
    end
    % Older mother
    nkids0to2=0;
    nkids3to7=0;
    nkids8to12=0;
    nkids13to18=0;
    agefirstchild=(jj+Params.agejshifter)-birthfirstchild_oldermother;
    agesecondchild=(jj+Params.agejshifter)-(birthfirstchild_oldermother+2); % second child born two year later
    if 0<=agefirstchild && agefirstchild<=2
        nkids0to2=nkids0to2+1;
    elseif 3<=agefirstchild && agefirstchild<=7
        nkids3to7=nkids3to7+1;
    elseif 8<=agefirstchild && agefirstchild<=12
        nkids8to12=nkids8to12+1;
    elseif 13<=agefirstchild && agefirstchild<=18
        nkids13to18=nkids13to18+1;
    end
    if 0<=agesecondchild && agesecondchild<=2
        nkids0to2=nkids0to2+1;
    elseif 3<=agesecondchild && agesecondchild<=7
        nkids3to7=nkids3to7+1;
    elseif 8<=agesecondchild && agesecondchild<=12
        nkids8to12=nkids8to12+1;
    elseif 13<=agesecondchild && agesecondchild<=18
        nkids13to18=nkids13to18+1;
    end
    Params.e.oldermother(jj)=1.67+0.233*nkids0to2+0.33*nkids3to7+0.4*nkids8to12+0.533*nkids13to18;
    if 0<=agefirstchild && agefirstchild<=7 % and childcarecosts
        Params.childcarecost.oldermother(jj)=Params.childcostcoeffs(agefirstchild+1);
    end
end
% Done calculating the consumption equivalence scale

% Because these consumption equivalence scales are a key part of model, let's plot them and the childcare costs
if CreateFigures==1
    figure(1)
    subplot(1,2,1); plot(Params.agejshifter+(1:1:N_j),Params.e.nochild)
    hold on
    subplot(1,2,1); plot(Params.agejshifter+(1:1:N_j),Params.e.youngmother)
    subplot(1,2,1); plot(Params.agejshifter+(1:1:N_j),Params.e.oldermother)
    hold off
    title('Consumption Equivalence Scales')
    xlabel('Age')
    legend('Non-mothers','Young mothers','Old mothers')
    subplot(1,2,2); plot(Params.agejshifter+(1:1:N_j),Params.childcarecost.nochild)
    hold on
    subplot(1,2,2); plot(Params.agejshifter+(1:1:N_j),Params.childcarecost.youngmother)
    subplot(1,2,2); plot(Params.agejshifter+(1:1:N_j),Params.childcarecost.oldermother)
    hold off
    title('Childcare Costs')
    xlabel('Age')
    legend('Non-mothers','Young mothers','Old mothers')
end

%% Female labor force experience, h_f
Params.delta=0.074; % human capital depreciation rate
Params.eta0=0.0266; % human capital appreciation rate (eta0+eta1*age)
Params.eta1=-0.000038; % human capital appreciation rate
Params.h_f_accum=Params.eta0+Params.eta1*Params.age;

vfoptions.experienceasset=1; % Using an experience asset
% Note: by default, assumes it is the last d variable that controls the
% evolution of the experience asset (and that the last a variable is
% the experience asset).

% aprimeFn gives the value of h_fprime
vfoptions.aprimeFn=@(P,h_f,delta, h_f_accum) exp(log(h_f)+h_f_accum*P-delta*(1-P));

% We also need to tell simoptions about the experience asset
simoptions.experienceasset=1;
simoptions.aprimeFn=vfoptions.aprimeFn;
% simoptions.a_grid=a_grid; % This one has to be done later once a_grid exists

%% Grids
d_grid=[0;1]; % female participation decision

asset_grid=4*linspace(0,1,n_a(1))'.^3; % grid on assets
% Note: initially had max assets around 10, reran a few times setting lower
% max assets as was clear cdf of assets reached 1 way before this
h_f_grid=linspace(1,1+sum(Params.h_f_accum),n_a(2))'; % grid on human capital
% ALSM2008: "However, we impose a floor on how far human capital can decline so that human capital will not fall below its initial value."
% While the paper does not mention anything the codes seem to set the initial value of h_f to 1.
% Note that we know what the max of h_f is (start at one, P=1 every period)

a_grid=[asset_grid; h_f_grid];

% unclear what ALSM2008 did for the unit-roots.
% The file: transwagepoints.inp contains 521601 lines, which appear to be transition probabilities of some kind
% The file: wagepoints_f.inp contains 720 numbers, low of about 0.1 and max about 5.
% The file: wagepoints_m.inp contains 480 numbers, e-10^4, or e-10^3
% Since this is in effect just a VAR(1), and since the initial value for
% everyone is zero, I am just going to treat it as a life-cycle VAR(1) and
% discretize it as such.
tauchenoptions.initialj1mewz=0; % everyone starts with z=0 in period j=1
tauchenoptions.nSigmas=2; % Make shock grids range just +-2 standard deviation [means that exp(min(z_grid_J)) is around 0.1 and exp(max(z_grid_J)) is around 4.5, which seems something like what ALSM2008 codes looks like it uses]
[z_grid_J,pi_z_J]=discretizeLifeCycleVAR1_Tauchen(Params.Mew_z.*ones(1,N_j),Params.Rho_z.*ones(1,1,N_j),Params.SigmaSq_zeta.*ones(1,1,N_j),n_z,N_j,tauchenoptions);
% Note: return fn uses exponential of this grid

%% Return function
DiscountFactorParamNames={'beta'};

ReturnFn=@(P,aprime,a,h_f,z_m,z_f,childcarecost,e,h_m,gamma,phi1,phi2,y_f0,y_m0,R)...
    AttanasioLowSanchezMarcos2008_ReturnFn(P,aprime,a,h_f,z_m,z_f,childcarecost,e,h_m,gamma,phi1,phi2,y_f0,y_m0,R);

%% Solve value function
vfoptions.ptypestorecpu=1; % To keep free space on the GPU we move the solutions to the CPU memory (I turned this on because it solved the first two PTypes and then gave an out of memory error on the third)
vfoptions.lowmemory=1;
vfoptions.verbose=1;
vfoptions.verboseparams=1
tic;
[V, Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z, N_j,Names_i,d_grid, a_grid, z_grid_J, pi_z_J, ReturnFn, Params, DiscountFactorParamNames, vfoptions);
timevf=toc

% % % Just check some things
% % tempP=Policy.nochild(1,:,:,:,:,:);
% % tempa=Policy.nochild(2,:,:,:,:,:);
% % [min(tempP(:)),max(tempP(:))] % 1 or 2
% % [min(tempa(:)),max(tempa(:))] % 1 to n_a(1)
% % 
% % tempPj=reshape(Policy.nochild(1,:,:,:,:,:),[prod(n_a)*prod(n_z),N_j]);
% % [min(tempPj,[],1);max(tempPj,[],1)] % 1 or 2
% % 
% % tempaj=reshape(Policy.nochild(2,:,:,:,:,:),[prod(n_a)*prod(n_z),N_j]);
% % [min(tempaj,[],1);max(tempaj,[],1)] % 1 to n_a(1)
% % 
% % sum(sum(sum(sum(sum(tempa>(9/10*n_a(1)))))))/numel(tempa)


%% Take a look at policy

if CreateFigures==1
    % load('ALSM2008.mat','Policy','asset_grid','h_f_grid')
    figure(2)
    subplot(3,1,1); plot(asset_grid,asset_grid(Policy.nochild(2,:,1,6,6,1)))
    hold on
    subplot(3,1,1); plot(asset_grid,asset_grid(Policy.nochild(2,:,1,6,6,11)))
    subplot(3,1,1); plot(asset_grid,asset_grid(Policy.nochild(2,:,1,6,6,21)))
    subplot(3,1,1); plot(asset_grid,asset_grid(Policy.nochild(2,:,1,6,6,31)))
    subplot(3,1,1); plot(asset_grid,asset_grid(Policy.nochild(2,:,1,6,6,41)))
    hold off
    title('Policy, no child: assetsprime as fn of assets (given h_f=1)')
    ylabel('Next period assets')
    xlabel('Assets')
    legend('period 1','period 11','period 21','period 31','period 41')

    subplot(3,1,2); plot(asset_grid,asset_grid(Policy.nochild(2,:,floor(n_a(2)/3),6,6,1)))
    hold on
    subplot(3,1,2); plot(asset_grid,asset_grid(Policy.nochild(2,:,floor(n_a(2)/3),6,6,11)))
    subplot(3,1,2); plot(asset_grid,asset_grid(Policy.nochild(2,:,floor(n_a(2)/3),6,6,21)))
    subplot(3,1,2); plot(asset_grid,asset_grid(Policy.nochild(2,:,floor(n_a(2)/3),6,6,31)))
    subplot(3,1,2); plot(asset_grid,asset_grid(Policy.nochild(2,:,floor(n_a(2)/3),6,6,41)))
    hold off
    title('Policy, no child: assetsprime as fn of assets (given h_f=index 1/3 of max)')
    ylabel('Next period assets')
    xlabel('Assets')
    legend('period 1','period 11','period 21','period 31','period 41')

    subplot(3,1,3); plot(asset_grid,asset_grid(Policy.nochild(2,:,floor(2*n_a(2)/3),6,6,1)))
    hold on
    subplot(3,1,3); plot(asset_grid,asset_grid(Policy.nochild(2,:,floor(2*n_a(2)/3),6,6,11)))
    subplot(3,1,3); plot(asset_grid,asset_grid(Policy.nochild(2,:,floor(2*n_a(2)/3),6,6,21)))
    subplot(3,1,3); plot(asset_grid,asset_grid(Policy.nochild(2,:,floor(2*n_a(2)/3),6,6,31)))
    subplot(3,1,3); plot(asset_grid,asset_grid(Policy.nochild(2,:,floor(2*n_a(2)/3),6,6,41)))
    hold off
    title('Policy, no child: assetsprime as fn of assets (given h_f=index 2/3 of max)')
    ylabel('Next period assets')
    xlabel('Assets')
    legend('period 1','period 11','period 21','period 31','period 41')
    

    figure(3)
    subplot(3,1,1); plot(asset_grid,d_grid(Policy.nochild(1,:,1,6,6,1)))
    hold on
    subplot(3,1,1); plot(asset_grid,d_grid(Policy.nochild(1,:,1,6,6,11)))
    subplot(3,1,1); plot(asset_grid,d_grid(Policy.nochild(1,:,1,6,6,21)))
    subplot(3,1,1); plot(asset_grid,d_grid(Policy.nochild(1,:,1,6,6,31)))
    subplot(3,1,1); plot(asset_grid,d_grid(Policy.nochild(1,:,1,6,6,41)))
    hold off
    title('Policy, no child: P as fn of assets (given h_f=1)')
    ylabel('Participation')
    xlabel('Assets')
    legend('period 1','period 11','period 21','period 31','period 41')

    subplot(3,1,2); plot(asset_grid,d_grid(Policy.nochild(1,:,floor(n_a(2)/3),6,6,1)))
    hold on
    subplot(3,1,2); plot(asset_grid,d_grid(Policy.nochild(1,:,floor(n_a(2)/3),6,6,11)))
    subplot(3,1,2); plot(asset_grid,d_grid(Policy.nochild(1,:,floor(n_a(2)/3),6,6,21)))
    subplot(3,1,2); plot(asset_grid,d_grid(Policy.nochild(1,:,floor(n_a(2)/3),6,6,31)))
    subplot(3,1,2); plot(asset_grid,d_grid(Policy.nochild(1,:,floor(n_a(2)/3),6,6,41)))
    hold off
    title('Policy, no child: P as fn of assets (given h_f=index 1/3 of max)')
    ylabel('Participation')
    xlabel('Assets')
    legend('period 1','period 11','period 21','period 31','period 41')

    subplot(3,1,3); plot(asset_grid,d_grid(Policy.nochild(1,:,floor(2*n_a(2)/3),6,6,1)))
    hold on
    subplot(3,1,3); plot(asset_grid,d_grid(Policy.nochild(1,:,floor(2*n_a(2)/3),6,6,11)))
    subplot(3,1,3); plot(asset_grid,d_grid(Policy.nochild(1,:,floor(2*n_a(2)/3),6,6,21)))
    subplot(3,1,3); plot(asset_grid,d_grid(Policy.nochild(1,:,floor(2*n_a(2)/3),6,6,31)))
    subplot(3,1,3); plot(asset_grid,d_grid(Policy.nochild(1,:,floor(2*n_a(2)/3),6,6,41)))
    hold off
    title('Policy, no child: P as fn of assets (given h_f=index 2/3 of max)')
    ylabel('Participation')
    xlabel('Assets')
    legend('period 1','period 11','period 21','period 31','period 41')


end


%% Agent distribution
simoptions.ptypestorecpu=1; % Tell simoptions we used this for vfoptions


% ALSM2008, mid page 1530: "We assume that 12 percent of women never have a
% child, 41 percent have the first of their two children at age 24, and 47
% percent have their first child at 29. [...] All women in our model begin
% life at age 23 with zero assets."

% It seems like all females start with h_f=1
% Based on ALSM2008 codes, there is a parameter 'inithk' which is set to
% imported value of parameter 'inithkfem' which is 1.

% Zero assets, guessing median shocks for men too?
jequaloneDist=zeros([n_a,n_z],'gpuArray');
jequaloneDist(1,1,ceil(n_z(1)/2),ceil(n_z(2)/2))=1; % no assets, h_f=1, median shocks
% Model does not contain any survival/death risk, so
AgeWeightsParamNames={'ageweights'};
Params.ageweights=ones(N_j,1)/N_j;
% Need to give weights to each of the three permanent types
PTypeDistParamNames={'ptypeweights'};
Params.ptypeweights=[0.12,0.41,0.47];

% Note: ALSM2008 never seem to mention the age weights used. So I have just
% guessed they weight them all equally.

simoptions.d_grid=d_grid; % because we are using an experience asset
simoptions.a_grid=a_grid; % because we are using an experience asset


simoptions.verbose=1;
tic;
StationaryDist=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightsParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_z,N_j,Names_i,pi_z_J,Params,simoptions);
timedist=toc

% Look at agent dist to better understand model and see if grids seem appropriate
if CreateFigures==1
    % load('ALSM2008.mat','StationaryDist','asset_grid','h_f_grid')
    figure(4)
    subplot(2,1,1); plot(asset_grid,cumsum(sum(sum(sum(sum(StationaryDist.nochild,5),4),3),2)))
    hold on
    subplot(2,1,1); plot(asset_grid,cumsum(sum(sum(sum(sum(StationaryDist.youngmother,5),4),3),2)))
    subplot(2,1,1); plot(asset_grid,cumsum(sum(sum(sum(sum(StationaryDist.oldermother,5),4),3),2)))
    hold off
    title('Cdf over assets')
    legend('Non-mothers','Young mothers','Old mothers')
    subplot(2,1,2); plot(h_f_grid,cumsum(sum(sum(sum(sum(StationaryDist.nochild,5),4),3),1)))
    hold on
    subplot(2,1,2); plot(h_f_grid,cumsum(sum(sum(sum(sum(StationaryDist.youngmother,5),4),3),1)))
    subplot(2,1,2); plot(h_f_grid,cumsum(sum(sum(sum(sum(StationaryDist.oldermother,5),4),3),1)))
    hold off
    title('Cdf over h_f')
    legend('Non-mothers','Young mothers','Old mothers')

    figure(5)
    subplot(2,1,1); plot(asset_grid,cumsum(sum(sum(sum(StationaryDist.nochild(:,:,:,:,1),4),3),2)))
    hold on
    subplot(2,1,1); plot(asset_grid,cumsum(sum(sum(sum(StationaryDist.nochild(:,:,:,:,11),4),3),2)))
    subplot(2,1,1); plot(asset_grid,cumsum(sum(sum(sum(StationaryDist.nochild(:,:,:,:,21),4),3),2)))
    subplot(2,1,1); plot(asset_grid,cumsum(sum(sum(sum(StationaryDist.nochild(:,:,:,:,31),4),3),2)))
    subplot(2,1,1); plot(asset_grid,cumsum(sum(sum(sum(StationaryDist.nochild(:,:,:,:,41),4),3),2)))
    hold off
    title('Cdf over assets: no child')
    legend('period 1','11','21','31','41')
end

%% Just while working on this
save ALSM2008.mat -v7.3
% so I can skip rerunning all the above
% load ALSM2008.mat

%% Calculate some aggregates
FnsToEvaluate.participation=@(P,aprime,a,h_f,z_m,z_f) P;
FnsToEvaluate.femaleshadowwage=@(P,aprime,a,h_f,z_m,z_f,y_f0) y_f0*h_f*exp(z_f);
FnsToEvaluate.femalewage=@(P,aprime,a,h_f,z_m,z_f,y_f0) y_f0*h_f*exp(z_f)*P; % note, this is female wage on average across households, with zeros for households where female is not working; so typically we will be interested in this divided by participation
FnsToEvaluate.malewage=@(P,aprime,a,h_f,z_m,z_f,y_m0,h_m) y_m0*h_m*exp(z_m);
FnsToEvaluate.h_f=@(P,aprime,a,h_f,z_m,z_f) h_f;
FnsToEvaluate.h_m=@(P,aprime,a,h_f,z_m,z_f,h_m) h_m;
FnsToEvaluate.exp_z_m=@(P,aprime,a,h_f,z_m,z_f) exp(z_m);
FnsToEvaluate.exp_z_f=@(P,aprime,a,h_f,z_m,z_f) exp(z_f);
FnsToEvaluate.z_m=@(P,aprime,a,h_f,z_m,z_f) z_m;
FnsToEvaluate.z_f=@(P,aprime,a,h_f,z_m,z_f) z_f;
FnsToEvaluate.assets=@(P,aprime,a,h_f,z_m,z_f) a;


AllStats=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid_J,simoptions);

AgeConditionalStats=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist, Policy, FnsToEvaluate, Params,n_d,n_a,n_z,N_j,Names_i,d_grid, a_grid, z_grid_J, simoptions);

save ALSM2008_2.mat AllStats AgeConditionalStats Params N_j
% load ALSM2008_2.mat
% This save file contains everything needed by the codes below this point

% To better understand model, look at h_f



% Second half of Table 3 of ALSM2008
T3_Participation=AllStats.participation.Mean;
% Participation of mothers with children under age of 3 is calculated below as part of Table 4
T3_Participation_nonmothers=AllStats.participation.nochild.Mean;
% MedianAnnualizedWageLossForNonPartipants % Not obvious to me how this is
% calculated, I guess you compare the wage when non-participants return to
% work with their wage before, then annualize this?
ObservedGenderWageGap=(AllStats.femalewage.Mean/AllStats.participation.Mean)/AllStats.malewage.Mean;
% ALSM2008 explain in (vi) on page 1532 how they do wage growth from a
% regression on simulated panel data. I am just going to calculate it a
% differnt way which would give a similar answer.
FemaleWages=AgeConditionalStats.femalewage.Mean./AgeConditionalStats.participation.Mean; % female wage (for those who participate)
FemaleWageGrowth_youngerthan35=mean(((FemaleWages(2:13)-FemaleWages(1:12))./FemaleWages(1:12))); % note, 13 corresponds to age 35
FemaleWageGrowth_olderthan35=mean(((FemaleWages(14:end)-FemaleWages(13:end-1))./FemaleWages(13:end-1))); % note, 13 corresponds to age 35
% Note: not clear what to do with 35 year olds in above two lines, but unlikely to matter at all anyway


% Table 4 f ALSM2008
% The 'nochild' type of woman never has children
% The 'youngmother' type of woman has her first child at 24, with a second child arriving two years after the first
% The 'oldermother' type of woman has her first child at 29, with a second child arriving two years after the first
% Params.agejshifter=22; % So initial age is 23
Participation_youngmother_childunder3=sum(AgeConditionalStats.participation.youngmother.Mean(2:6).*Params.ageweights(2:6)'); % Note: ages 24-28; first child at 24, second at 26, so there is a child under 3 until she is 29 years old
Participation_oldermother_childunder3=sum(AgeConditionalStats.participation.youngmother.Mean(7:11).*Params.ageweights(7:11)'); % Note: ages 29-33; first child at 29, second at 31, so there is a child under 3 until she is 34 years old
Participation_mother_childunder3=Params.ptypeweights(2)*Participation_youngmother_childunder3+Params.ptypeweights(3)*Participation_oldermother_childunder3;
Participation_youngmother_child4to18=sum(AgeConditionalStats.participation.youngmother.Mean(6:22).*Params.ageweights(6:22)'); % Note: ages 24-28; first child at 24, second at 26, so there is a child 4-18 when she is 28-44 years old
Participation_oldermother_child4to18=sum(AgeConditionalStats.participation.youngmother.Mean(11:27).*Params.ageweights(11:27)'); % Note: ages 29-33; first child at 29, second at 31, so there is a child 4-18 when she is 33-49 years old
Participation_mother_child4to18=Params.ptypeweights(2)*Participation_youngmother_child4to18+Params.ptypeweights(3)*Participation_oldermother_child4to18;

% I skip the Median duration of return to work, easy to do but I would need to simulate panel data and I can't be bothered just now
% Plus they say 7 years, but doing a similar calculation their Table 1
% reports much smaller numbers, so not clear exactly what calculation
% underlies this.

fprintf('Some numbers relating to second half of Table 3 of ALSM2008 \n')
fprintf('Participation                                            %8.2f \n', T3_Participation)
fprintf('Participation mothers with children under age of 3       %8.2f \n', Participation_mother_childunder3)
fprintf('Participation nonmothers                                 %8.2f \n', T3_Participation_nonmothers)
fprintf(' - \n')
fprintf('Observed wage-gender gap                                 %8.2f \n', ObservedGenderWageGap)
fprintf('Female wage growth (if younger than 35)                  %8.3f \n', FemaleWageGrowth_youngerthan35)
fprintf('Female wage growth (if older than 35)                    %8.3f \n', FemaleWageGrowth_olderthan35)
% For comparison, in ALSM2008 paper these six numbers in same order are
% 0.72, 0.47, 0.83
% 0.69, 0.025, 0.006

fprintf(' \n')
fprintf('Some numbers relating to Table 4 of ALSM2008 \n')
fprintf('Participation mothers with children under age of 3           %8.2f \n', Participation_mother_childunder3)
fprintf('Participation young mothers with children under age of 3     %8.2f \n', Participation_youngmother_childunder3)
fprintf('Participation old mothers with children under age of 3       %8.2f \n', Participation_oldermother_childunder3)
fprintf('Participation mothers with children age 4-18                 %8.2f \n', Participation_mother_child4to18)
fprintf(' - \n')
% For comparison, in ALSM2008 paper these four numbers in same order are
% 0.47, 0.42, 0.53, 0.70



%% Plot the employment (participation) rates
% Figure 10 of ALSM2008

if CreateFigures==1
    figure(10)
    plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.participation.nochild.Mean)
    hold on
    plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.participation.youngmother.Mean)
    plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.participation.oldermother.Mean)
    plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.participation.Mean)
    hold off
    title('Participation Profiles by Age at Childbirth')
    ylabel('Employment Rate')
    xlabel('Age')
    legend('Non-mothers','Young mothers','Old mothers','All')
end



%% Some plots to better understand the model
if CreateFigures==1
    figure(6)
    subplot(3,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.malewage.nochild.Mean)
    hold on
    subplot(3,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.malewage.youngmother.Mean)
    subplot(3,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.malewage.oldermother.Mean)
    hold off
    title('Mean Male Earnings (conditional on age)')
    ylabel('Earnings')
    xlabel('Age')
    legend('Non-mothers','Young mothers','Old mothers')

    subplot(3,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.femalewage.nochild.Mean)
    hold on
    subplot(3,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.femalewage.youngmother.Mean)
    subplot(3,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.femalewage.oldermother.Mean)
    hold off
    title('Mean Female Earnings (conditional on age; includes zeros for non-working)')
    ylabel('Earnings')
    xlabel('Age')
    legend('Non-mothers','Young mothers','Old mothers')

    subplot(3,1,3); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.assets.nochild.Mean)
    hold on
    subplot(3,1,3); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.assets.youngmother.Mean)
    subplot(3,1,3); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.assets.oldermother.Mean)
    hold off
    title('Mean Assets (conditional on age)')
    ylabel('Earnings')
    xlabel('Age')
    legend('Non-mothers','Young mothers','Old mothers')
end

% The role of components of earnings
if CreateFigures==1
    figure(7)
    subplot(4,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.exp_z_m.nochild.Mean)
    hold on
    subplot(4,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.exp_z_m.youngmother.Mean)
    subplot(4,1,1); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.exp_z_m.oldermother.Mean)
    hold off
    title('Mean exp(z_m)')
    xlabel('Age')
    legend('Non-mothers','Young mothers','Old mothers')
    subplot(4,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.exp_z_f.nochild.Mean)
    hold on
    subplot(4,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.exp_z_f.youngmother.Mean)
    subplot(4,1,2); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.exp_z_f.oldermother.Mean)
    hold off
    title('Mean exp(z_f)')
    xlabel('Age')
    legend('Non-mothers','Young mothers','Old mothers')
    subplot(4,1,3); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.h_m.nochild.Mean)
    hold on
    subplot(4,1,3); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.h_m.youngmother.Mean)
    subplot(4,1,3); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.h_m.oldermother.Mean)
    hold off
    title('Mean h_m')
    xlabel('Age')
    legend('Non-mothers','Young mothers','Old mothers')
    subplot(4,1,3); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.h_f.nochild.Mean)
    hold on
    subplot(4,1,3); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.h_f.youngmother.Mean)
    subplot(4,1,3); plot(Params.agejshifter+(1:1:N_j),AgeConditionalStats.h_f.oldermother.Mean)
    hold off
    title('Mean h_f')
    xlabel('Age')
    legend('Non-mothers','Young mothers','Old mothers')
end








