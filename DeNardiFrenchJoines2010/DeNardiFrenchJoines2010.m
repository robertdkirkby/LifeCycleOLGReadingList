% De Nardi, French & Joines (2010) - Why do the Elderly Save? The Role of Medical Expenses

% Because VFI Toolkit does not presently have a way to allow survival
% probabilities to depend on the state of the model, we instead have to
% implement this via adding a 'dead' value to the health state, so it takes 
% three values instead of two. To be able to use this probability to
% determine the warm-glow of bequests, we have to further differentiate
% 'death' from 'dead', so the health state now takes four values, and
% 'death' indicates the first period of being dead in which a warm-glow is
% received.

% A few notational changes: I use j for age (DFJ2010 use t). I use q_c to
% index the income quintiles and q1,q2,...,q5 as names for them (DFJ2010 use I).

%% I CANNOT FIND THREE PARAMETER VALUES. JUST USING PLACEHOLDERS FOR THE MOMENT
% tau_e and estateexemption (tautilde and xtilde in DFJ2010 notation)
% The estate tax rate and the exemption amount for the estate tax.
% r, the interest rate.
% I have made up numbers: tau_e=0.05, estateexpemption=1000000, r=0.05;

%% Setup the model action and state space
n_d=0;
n_a=501; % assets
n_z=[4,15]; % h (health status, 0=good, 1=bad, 2=death, 3=dead) and zeta (markov for medical expense shocks)
n_e=7; % xi (iid for medical expense shocks)
N_j=33; % 31 is 70 to 100
% Permanent types: 2 genders, 5 income groups (quintiles of permanent income)
Names_i={'femaleq1','femaleq2','femaleq3','femaleq4','femaleq5',...
    'maleq1','maleq2','maleq3','maleq4','maleq5'}; % q1=first (bottom) quintile, q5=fifth (top) quintile

% Note: DFJ2010 used 8 points for zeta and 8 points for xi (lines 68,69 of their wealth20_Singles.c)

% Near the beginning of DFJ2010 wealth20_Singles.c code, there is a comment
%    "- We need to check if the difference in indexes exists also for singles (as it did for couples)
%    Note the different indexes for males and females: for a male age t=1, 
%    means age=70,  for a female t=1 means age=67. Both 
%    live up to age 100 (t=31 for males and 34 for females)."
% Reading through the code it seems to use same index for age in various
% matrices regardless of gender, so pretty sure this was resolved and model
% is just ages 70 to 100.

figure_c=0; % counter for figures

%% Parameters
% Six are just initial guesses as they will be GMM estimated later:
% upsilon, delta, theta, k, cfloor, 

% Demographics
Params.agejshifter=69;
Params.J=N_j;

% Discount factor (initial guess, will be estimated)
Params.beta=0.97;
% Utility fn (initial guesses, will be estimated)
Params.upsilon=3.66; % CES utility param
Params.delta=-0.36; % utility cost of bad health
% Warm-glow of bequests (initial guesses, will be estimated)
Params.theta=2419; % intensity of bequest motive
Params.k=215; % determines curvature of warm-glow in bequests, and hence the extent to which bequests are a luxury good

% Interest rate
Params.r=0.05;

% Goverment transfers (initial guess, will be estimated)
Params.cfloor=2653;

% Income taxes (from lines 189-194 of their wealth20_Singles.c)
% Based on tax brackets and marginal tax rates for each bracket
Params.taxbracket1=6250;
Params.taxbracket2=40200;
Params.taxbracket3=68400;
Params.taxbracket4=93950;
Params.taxbracket5=148250;
Params.taxbracket6=284700;
Params.margtaxrate0=0.0765; % tax rate 0 is for < bracket 1
Params.margtaxrate1=0.2616; % tax rate 1 is for brackets 1-to-2
Params.margtaxrate2=0.4119;
Params.margtaxrate3=0.3499;
Params.margtaxrate4=0.3834;
Params.margtaxrate5=0.4360;
Params.margtaxrate6=0.4761; 

% Estate taxes
Params.tau_e=0.05; % estate tax rate (DFJ2010 call this tautilde)
Params.estateexemption=1000000; % estate tax rate (DFJ2010 call this xtilde)
% These are not in paper, got them from DFJ2010 codes.
% line 1782-3:
% tauBeq = asstvecPtr[1]; /* Estate tax rate */
% exBeq  = asstvecPtr[2]; /* Estate tax exemption level */
% which is read from a file, line 394
% asstvecPtr    = gau5read(strcat(strcpy(fullpath,rootdir),"asstvec.fmt"));
% But "asstvec.fmt" is not in replication materials, for now I just use a
% placeholder.

% Medical expenses 1/2. (from Table 1):
% Two stochastic components: zeta which is AR(1) markov, and xi which is iid
% Note: these later get 'scaled up' by the sigma_coeff
Params.rho_m=0.922; % autocorrelation, persistent part
Params.sigma_mepsilon=sqrt(0.050); % std dev of innovations, persistent part
Params.sigma_mxi=sqrt(0.665); % std dev of innovations, transitory part
% Footnote 3: We assume that medical expenses are log-normally distributed, 
% so the predicted level of medical expenses is exp (m <= 1/2*sigma^2), where 
% m denotes predicted log medical expenses and sigma^2 denotes the variance of 
% the idiosyncratic shock psi. The ratio of the level of medical expenses 
% two standard deviations above the mean to average medical expenses is
%    exp(m+2*sigma)/exp(m+1/2*sigma^2) =exp(2*sigma-1/2*sigma^2)
% which =6.8 if sigma=sqrt(2.53)
% Added Note: I'm sure why they use 2.23 as being a relevant number?

% Medical expenses 2/2.
% Two deterministic coefficients, m(g,h,I,j) and sigma(g,h,I,j) from eqn (6) of DFJ2010
% Call them m_coeff and sigma_coeff, and due to the dependence on
% health status I create four parameters (m_coeff_healthgood, m_coeff_healthbad, 
% sigma_coeff_healthgood, sigma_coeff_healthbad) each of which is
% age-and-ptype dependent. Then in ReturnFn use health status in
% if-statement to decide which of _healthgood and _healthbad to use).
% sigma_coeff scales the standard deviation of the stochastic component
% This is setup below


%% Health transition probabilities
% copy-paste from DFJ2010 docs:
% "the parameters of the health transition matrix (healthprof.out)"
% "file describing how to map the two-year transition probabilities that we estimate into one-year transition probabilities (two_to_one.pdf)"

% Start with getting an interpreting data from healthprof.out (copy-paste from DFJ2010 docs):
% "The first column in healthprof.out is age1, which is just an index.
% The remaining columns are logit coefficients in the equation:
% Pr(health_t= bad | X_{t – 2 years}) = exp(X_t b)/(1+exp(X_t b)), 
% The coefficients are as follows: 
% ageshiftFINAL = age-specific intercept term
% healshift = age-specific coefficients for being in bad health 2 years ago (good health is the omitted category)
% maleshift = age-specific coefficients for being male (female is the omitted category)
% PIshift = age-specific coefficients on PI, which is measured as the individual s percentile rank (0 to 1) in the PI distribution 
% bPI2 = coefficients on PI percentile squared.  For now there is no variation across ages, but that can be added.
% Suppose you had an 80-year-old male who was in good health and at the 60th percentile of the PI distribution.
% The probability that this man will be in bad health at age 82 is: 
% Xb = ageshiftFINAL82+maleshift82+PI82*.6+bPI82*.36, 
% where the "82” subscript denotes the row with age index number 82."

% Following is copy-pasted from contents of healthprof.out
healthprof.age=70:1:102;
healthprof.ageshiftFinal=[-.7599814,-.7719013,-.7743046,-.7680932,-.7541689,-.7334337,-.7067895,-.675138,-.6393813, -.6004211, -.5591595, -.5164981, -.473339, -.4305839, -.3891349, -.3498937, -.3137623, -.2816426, -.2544363, -.2330454, -.2183719, -.2113174, -.2127841, -.2236736, -.2448879, -.277329, -.3218986, -.3794986, -.451031, -.5373976, -.6395003, -.7582409, -.8945215]; 
healthprof.healshift=[2.600759,2.55415,2.507542,2.460934,2.414325,2.367717,2.321109,2.274501,2.227893, 2.181284,2.134676, 2.088068, 2.04146, 1.994851, 1.948243, 1.901635, 1.855027, 1.808418, 1.76181,1.715202,1.668594, 1.621985, 1.575377, 1.528769, 1.482161, 1.435552, 1.388944, 1.342336, 1.295728, 1.249119, 1.202511, 1.155903, 1.109295];
healthprof.maleshift=[.3348944,.3210961,.3072977,.2934994,.2797011,.2659027, .2521044, .2383061, .2245077, .2107094, .1969111, .1831127, .1693144, .1555161, .1417177, .1279194, .1141211, .1003227, .0865244, .0727261, .0589277, .0451294, .0313311, .0175327, .0037344, -.0100639, -.0238623, -.0376606, -.0514589, -.0652573, -.0790556,-.0928539, -.1066523];
healthprof.PIshift=[-1.270953, -1.259267, -1.247581, -1.235895, -1.224209, -1.212523, -1.200836, -1.18915, -1.177464, -1.165778, -1.154092, -1.142406, -1.13072, -1.119033, -1.107347, -1.095661, -1.083975, -1.072289, -1.060603, -1.048917, -1.037231, -1.025545, -1.013858, -1.002172, -.9904861, -.9788, -.9671139, -.9554277, -.9437416, -.9320555, -.9203693, -.9086832, -.8969971];
healthprof.bPI2=.3519772*ones(1,33);

% So need to use these to evaluate Xb, then use this to evaluate
% Pr(health=bad|X, two years ago); note gives Prob for 'current age'
probbadhealth_2yr=zeros(2,102-70+1,2,5); % current h, age, gender, income q
for h_c=1:2
    for jj=70:102
        for g_c=1:2 % female, male
            for q_c=1:5 % income quintile
                age=jj;
                ageindex=jj-69;
                cpercentile=(20*(q_c-1)+10)/100; % should be 0.6 is 60th percentile; I identify each quintile with its midpoint, so bottom quintile is 10th percentile [not clear what DFJ2010 did]
                Xb=healthprof.ageshiftFinal(ageindex)+healthprof.healshift(ageindex)*(h_c-1)+healthprof.maleshift(ageindex)*(g_c-1)+healthprof.PIshift(ageindex)*cpercentile+healthprof.bPI2(ageindex)*(cpercentile);
                probbadhealth_2yr(h_c,ageindex,g_c,q_c)=exp(Xb)/(1+exp(Xb));
            end
        end
    end
end

% Convert two-year transitions into one-year transitions.
% Idea: interpret two-year transtion from good-to-good health (gg) and two-year 
% transition from bad-to-bad health (bb) as coming from annual (g) and (b). 
% Then a bit of algebra, spelled out in two-to-one.pdf (of
% DJF2010 replication materials) gives us this following.
% Let 
%   A=2-gg-bb
%   B=2*(gg-1)
%   C=1-bb-gg+bb^2
% and then we get that 
%   b=(-B+sqrt(B^2-4AC))/2A  [note, quadratic formula]
% and we can calculate g=(1-bb-b*(1-b))/(1-b)
% As double checks that solution is correct, we should have:
%   gg=g^2+(1-g)*(1-b)

probbadhealth_2yr=zeros(2,102-70+1,2,5); % current h, age, gender, income q
probhealth=zeros(2,2,102-70+1,2,5); % previous h, current h, age, gender, income q
for jj=70:102
    for g_c=1:2 % female, male
        for q_c=1:5 % income quintile
            age=jj;
            ageindex=jj-69;

            hlag_c=1; % h=0 good health
            gg=1-probbadhealth_2yr(hlag_c,ageindex,g_c,q_c);
            hlag_c=2; % h=1 bad health
            bb=probbadhealth_2yr(hlag_c,ageindex,g_c,q_c);

            A=2-gg-bb;
            B=2*(gg-1);
            C=1-bb-gg+bb^2;
            b=(-B+sqrt(B^2-4*A*C))/(2*A);
            % balt=(-B-sqrt(B^2-4AC))/2A % there is another root, but
            % the pdf 'two_to_one.pdf' suggests this is not used
            g=(1-bb-b*(1-b))/(1-b);

            hlag_c=1; % h=0 good health
            if h_c==1
                probhealth(hlag_c,h_c,ageindex,g_c,q_c)=g; % good-to-good
            elseif h_c==2
                probhealth(hlag_c,h_c,ageindex,g_c,q_c)=1-g; % good-to-bad
            end
            hlag_c=2; % h=0 good health
            if h_c==1
                probhealth(hlag_c,h_c,ageindex,g_c,q_c)=1-b; % bad-to-good
            elseif h_c==2
                probhealth(hlag_c,h_c,ageindex,g_c,q_c)=b; % bad-to-bad
            end

        end
    end
end
% Note: not entirely clear from DFJ2010 docs which 2yr period becomes which
% 1yr transition (conversion is clear, just not which period we allocate
% the result to). I assume that it is based on the age that we are looking
% at transitions to (not from).

% Done

%% Conditional survival probabilities
% parameters of the survival probability transitions are in deathprof.out

% Following is copy-paste from DFJ2010 docs:
% 2. The line at the bottom of death.do will help you understand deathprof.out
% it works just like healthprof.out: 
% order age1 ageshiftFINAL healshift maleshift PIshift bPI2; 
% The first column in deathprof.out is age1, which is an index.  
% The remaining columns are logit coefficients in the equation:
% Pr(alive at time t |Xt – 2 years, alive at time t – 2 years) = exp(Xtb)/(1+exp(Xtb)), 
% These coefficients are completely analogous to the coefficients of the health 
% transition equation.
% 3. Note that these are 2-year survival probabilities – take the square root of the 2-year survival probability to get the 1-year survival probability.

% Following is copy-pasted from contents of deathprof.out
deathprof.age=70:1:102;
deathprof.ageshiftFinal=[2.546009, 2.565764, 2.575288, 2.574989, 2.565275, 2.546554, 2.519237, 2.483731, 2.440444, 2.389787, 2.332167, 2.267992, 2.197673, 2.121617, 2.040233, 1.95393, 1.863116, 1.7682, 1.669591, 1.567698, 1.462929, 1.355692, 1.246397, 1.135452, 1.023266, .9102477, .7968054, .6833481, .5702842, .4580225, .3469718, .2375407, .1301378];
deathprof.healshift=[-1.026724, -1.012185, -.9976472, -.9831089, -.9685707, -.9540324, -.9394941, -.9249559, -.9104176, -.8958793, -.8813411, -.8668028, -.8522645, -.8377263, -.823188, -.8086497, -.7941114, -.7795732, -.7650349, -.7504966, -.7359584, -.7214201, -.7068818, -.6923436, -.6778053, -.663267, -.6487288, -.6341905, -.6196522, -.605114, -.5905757, -.5760374, -.5614991];
deathprof.maleshift=[-1.017154, -.9981147, -.9790752, -.9600356, -.9409961, -.9219566, -.902917, -.8838775, -.864838, -.8457984, -.8267589, -.8077194, -.7886798, -.7696403, -.7506008, -.7315612, -.7125217, -.6934822, -.6744426, -.6554031, -.6363636, -.6173241, -.5982845, -.579245, -.5602055, -.5411659, -.5221264, -.5030869, -.4840473, -.4650078, -.4459683, -.4269287, -.4078892];
deathprof.PIshift=[.879185, .838859, .7985329, .7582068, .7178808, .6775548, .6372288, .5969027, .5565767, .5162507, .4759246, .4355986, .3952726, .3549465, .3146205, .2742945, .2339684, .1936424, .1533164, .1129903, .0726643, .0323383, -.0079878, -.0483138, -.0886398, -.1289659, -.1692919, -.2096179, -.249944, -.29027, -.3305961, -.3709221, -.4112481];
deathprof.bPI2=.1008181*ones(1,33);

% So need to use these to evaluate Xb, then use this to evaluate
% Pr(alive at t|X, two years ago); note gives Prob for 'current age'
probalive_2yr=zeros(2,102-70+1,2,5); % current h, age, gender, income q
for h_c=1:2
    for jj=70:102
        for g_c=1:2 % female, male
            for q_c=1:5 % income quintile
                age=jj;
                ageindex=jj-69;
                cpercentile=(20*(q_c-1)+10)/100; % should be 0.6 is 60th percentile; I identify each quintile with its midpoint, so bottom quintile is 10th percentile [not clear what DFJ2010 did]
                Xb=deathprof.ageshiftFinal(ageindex)+deathprof.healshift(ageindex)*(h_c-1)+deathprof.maleshift(ageindex)*(g_c-1)+deathprof.PIshift(ageindex)*cpercentile+deathprof.bPI2(ageindex)*(cpercentile);
                probalive_2yr(h_c,ageindex,g_c,q_c)=exp(Xb)/(1+exp(Xb));
            end
        end
    end
end

% Convert to annual by taking the square root
probalive_1yr=sqrt(probalive_2yr); % current h, age, gender, income q


%% Combine health probabilities and conditional survival probabilities to
% construct the transtion matrix for pi_h_J for each ptype.

pi_h_J.femaleq1=zeros(n_z(1),n_z(1),N_j);
g_c=1; % female
q_c=1; % 1st income quintile
for hlag_c=1:2
    for jj=1:N_j
        % conditional survival probabilities
        sj=1-probalive_1yr(hlag_c,jj,g_c,q_c); % current h, age+1, gender, income q
        % health transitions
        pi_h_J.femaleq1(hlag_c,3,jj)=(1-sj); % probability of dying
        for h_c=1:2
            pi_h_J.femaleq1(hlag_c,h_c,jj)=sj*probhealth(hlag_c,h_c,jj,g_c,q_c); % prob of survival*prod of health status
        end
    end
end
% Now just add the death-to-dead, and make dead absorbing
pi_h_J.femaleq1(3:4,4,:)=1;

% Rest is just repeating this for the other ptypes
pi_h_J.femaleq2=zeros(n_z(1),n_z(1),N_j);
g_c=1; % female
q_c=2; % 2nd income quintile
for hlag_c=1:2
    for jj=1:N_j
        % conditional survival probabilities
        sj=1-probalive_1yr(hlag_c,jj,g_c,q_c); % current h, age+1, gender, income q
        % health transitions
        pi_h_J.femaleq2(hlag_c,3,jj)=(1-sj); % probability of dying
        for h_c=1:2
            pi_h_J.femaleq2(hlag_c,h_c,jj)=sj*probhealth(hlag_c,h_c,jj,g_c,q_c); % prob of survival*prod of health status
        end
    end
end
% Now just add the death-to-dead, and make dead absorbing
pi_h_J.femaleq2(3:4,4,:)=1;

pi_h_J.femaleq3=zeros(n_z(1),n_z(1),N_j);
g_c=1; % female
q_c=3; % 3rd income quintile
for hlag_c=1:2
    for jj=1:N_j
        % conditional survival probabilities
        sj=1-probalive_1yr(hlag_c,jj,g_c,q_c); % current h, age+1, gender, income q
        % health transitions
        pi_h_J.femaleq3(hlag_c,3,jj)=(1-sj); % probability of dying
        for h_c=1:2
            pi_h_J.femaleq3(hlag_c,h_c,jj)=sj*probhealth(hlag_c,h_c,jj,g_c,q_c); % prob of survival*prod of health status
        end
    end
end
% Now just add the death-to-dead, and make dead absorbing
pi_h_J.femaleq3(3:4,4,:)=1;

pi_h_J.femaleq4=zeros(n_z(1),n_z(1),N_j);
g_c=1; % female
q_c=4; % 4th income quintile
for hlag_c=1:2
    for jj=1:N_j
        % conditional survival probabilities
        sj=1-probalive_1yr(hlag_c,jj,g_c,q_c); % current h, age+1, gender, income q
        % health transitions
        pi_h_J.femaleq4(hlag_c,3,jj)=(1-sj); % probability of dying
        for h_c=1:2
            pi_h_J.femaleq4(hlag_c,h_c,jj)=sj*probhealth(hlag_c,h_c,jj,g_c,q_c); % prob of survival*prod of health status
        end
    end
end
% Now just add the death-to-dead, and make dead absorbing
pi_h_J.femaleq4(3:4,4,:)=1;

pi_h_J.femaleq5=zeros(n_z(1),n_z(1),N_j);
g_c=1; % female
q_c=5; % 5th income quintile
for hlag_c=1:2
    for jj=1:N_j
        % conditional survival probabilities
        sj=1-probalive_1yr(hlag_c,jj,g_c,q_c); % current h, age+1, gender, income q
        % health transitions
        pi_h_J.femaleq5(hlag_c,3,jj)=(1-sj); % probability of dying
        for h_c=1:2
            pi_h_J.femaleq5(hlag_c,h_c,jj)=sj*probhealth(hlag_c,h_c,jj,g_c,q_c); % prob of survival*prod of health status
        end
    end
end
% Now just add the death-to-dead, and make dead absorbing
pi_h_J.femaleq5(3:4,4,:)=1;


% Now do male
pi_h_J.maleq1=zeros(n_z(1),n_z(1),N_j);
g_c=2; % male
q_c=1; % 1st income quintile
for hlag_c=1:2
    for jj=1:N_j
        % conditional survival probabilities
        sj=1-probalive_1yr(hlag_c,jj,g_c,q_c); % current h, age+1, gender, income q
        % health transitions
        pi_h_J.maleq1(hlag_c,3,jj)=(1-sj); % probability of dying
        for h_c=1:2
            pi_h_J.maleq1(hlag_c,h_c,jj)=sj*probhealth(hlag_c,h_c,jj,g_c,q_c); % prob of survival*prod of health status
        end
    end
end
% Now just add the death-to-dead, and make dead absorbing
pi_h_J.maleq1(3:4,4,:)=1;

pi_h_J.maleq2=zeros(n_z(1),n_z(1),N_j);
g_c=2; % male
q_c=2; % 2nd income quintile
for hlag_c=1:2
    for jj=1:N_j
        % conditional survival probabilities
        sj=1-probalive_1yr(hlag_c,jj,g_c,q_c); % current h, age+1, gender, income q
        % health transitions
        pi_h_J.maleq2(hlag_c,3,jj)=(1-sj); % probability of dying
        for h_c=1:2
            pi_h_J.maleq2(hlag_c,h_c,jj)=sj*probhealth(hlag_c,h_c,jj,g_c,q_c); % prob of survival*prod of health status
        end
    end
end
% Now just add the death-to-dead, and make dead absorbing
pi_h_J.maleq2(3:4,4,:)=1;

pi_h_J.maleq3=zeros(n_z(1),n_z(1),N_j);
g_c=2; % male
q_c=3; % 3rd income quintile
for hlag_c=1:2
    for jj=1:N_j
        % conditional survival probabilities
        sj=1-probalive_1yr(hlag_c,jj,g_c,q_c); % current h, age+1, gender, income q
        % health transitions
        pi_h_J.maleq3(hlag_c,3,jj)=(1-sj); % probability of dying
        for h_c=1:2
            pi_h_J.maleq3(hlag_c,h_c,jj)=sj*probhealth(hlag_c,h_c,jj,g_c,q_c); % prob of survival*prod of health status
        end
    end
end
% Now just add the death-to-dead, and make dead absorbing
pi_h_J.maleq3(3:4,4,:)=1;

pi_h_J.maleq4=zeros(n_z(1),n_z(1),N_j);
g_c=2; % male
q_c=4; % 4th income quintile
for hlag_c=1:2
    for jj=1:N_j
        % conditional survival probabilities
        sj=1-probalive_1yr(hlag_c,jj,g_c,q_c); % current h, age+1, gender, income q
        % health transitions
        pi_h_J.maleq4(hlag_c,3,jj)=(1-sj); % probability of dying
        for h_c=1:2
            pi_h_J.maleq4(hlag_c,h_c,jj)=sj*probhealth(hlag_c,h_c,jj,g_c,q_c); % prob of survival*prod of health status
        end
    end
end
% Now just add the death-to-dead, and make dead absorbing
pi_h_J.maleq4(3:4,4,:)=1;

pi_h_J.maleq5=zeros(n_z(1),n_z(1),N_j);
g_c=2; % male
q_c=5; % 5th income quintile
for hlag_c=1:2
    for jj=1:N_j
        % conditional survival probabilities
        sj=1-probalive_1yr(hlag_c,jj,g_c,q_c); % current h, age+1, gender, income q
        % health transitions
        pi_h_J.maleq5(hlag_c,3,jj)=(1-sj); % probability of dying
        for h_c=1:2
            pi_h_J.maleq5(hlag_c,h_c,jj)=sj*probhealth(hlag_c,h_c,jj,g_c,q_c); % prob of survival*prod of health status
        end
    end
end
% Now just add the death-to-dead, and make dead absorbing
pi_h_J.maleq5(3:4,4,:)=1;


%% Income
% (i) the parameters of the income regression (incprof.out)
% The left-hand variable in the income regression is the log of income.  
% The right-hand-side variables are just like those for survival and health:
% age ageshiftFINAL healshift maleshift PIshift bPI2

% Following is copy-pasted from contents of incprof.out
incprof.age=70:1:102;
incprof.ageshiftFinal=[8.024719, 8.031113, 8.035824, 8.039002, 8.040796, 8.041354, 8.040826, 8.039358, 8.037103, 8.034205, 8.030817, 8.027085, 8.023159, 8.019187, 8.015318, 8.011701, 8.008484, 8.005816, 8.003847, 8.002726, 8.002598, 8.003615, 8.005926, 8.009678, 8.01502, 8.022102, 8.031073, 8.042079, 8.055271, 8.070798, 8.088807, 8.109447, 8.132869];
incprof.healshift=0*ones(1,33);
incprof.maleshift=.019004*ones(1,33);
incprof.PIshift=2.458653*ones(1,33);
incprof.bPI2=-.4810466*ones(1,33);

% Create Params from these for use in model
% healshift is 0, so omit from formula in the following (note, female=0, male=1 in the regression of DFJ2010)
cpercentile=10/100; % q1
Params.earnings.femaleq1=exp(incprof.ageshiftFinal+incprof.maleshift*0+incprof.PIshift*cpercentile+incprof.bPI2*cpercentile^2);
Params.earnings.maleq1=exp(incprof.ageshiftFinal+incprof.maleshift*1+incprof.PIshift*cpercentile+incprof.bPI2*cpercentile^2);
cpercentile=30/100; % q2
Params.earnings.femaleq2=exp(incprof.ageshiftFinal+incprof.maleshift*0+incprof.PIshift*cpercentile+incprof.bPI2*cpercentile^2);
Params.earnings.maleq2=exp(incprof.ageshiftFinal+incprof.maleshift*1+incprof.PIshift*cpercentile+incprof.bPI2*cpercentile^2);
cpercentile=50/100; % q3
Params.earnings.femaleq3=exp(incprof.ageshiftFinal+incprof.maleshift*0+incprof.PIshift*cpercentile+incprof.bPI2*cpercentile^2);
Params.earnings.maleq3=exp(incprof.ageshiftFinal+incprof.maleshift*1+incprof.PIshift*cpercentile+incprof.bPI2*cpercentile^2);
cpercentile=70/100; % q4
Params.earnings.femaleq4=exp(incprof.ageshiftFinal+incprof.maleshift*0+incprof.PIshift*cpercentile+incprof.bPI2*cpercentile^2);
Params.earnings.maleq4=exp(incprof.ageshiftFinal+incprof.maleshift*1+incprof.PIshift*cpercentile+incprof.bPI2*cpercentile^2);
cpercentile=90/100; % q5
Params.earnings.femaleq5=exp(incprof.ageshiftFinal+incprof.maleshift*0+incprof.PIshift*cpercentile+incprof.bPI2*cpercentile^2);
Params.earnings.maleq5=exp(incprof.ageshiftFinal+incprof.maleshift*1+incprof.PIshift*cpercentile+incprof.bPI2*cpercentile^2);


%% Medical expenses: m_coeff and sigma_coeff
% (i) the parameters of the medical expense regression where the left hand 
% side variable is average medical expenses over two years (medexprofX.out); 
% (ii) the parameters of the medical expense regression, where the left hand 
% side variable is annual medical expenses (medexprof_adj.out)

% Copy-paste from DFJ2010 docs
% The list of variables in medexprofX.out and medexprof_adj.out is:
% age ageshiftFINAL healshift maleshift PIshift bPI2 vage vheal vmale vPI vPI2
% This list contains the results from two regressions.  
% The first regression is structured identically to the income regression, 
% except that the left-hand-side variable is log medical expenses.  
% The results of this regression are in the variables
% ageshiftFINAL healshift maleshift PIshift bPI2 
% which are just as before.  
% The second regression is a regression of the variance of medical expenses 
% on the same covariates.  The results of this regression are in the variables
% vage vheal vmale vPI vPI2
% where: 'vage' is an age-specific intercept for the (conditional on X) variance 
% of log medical expenditures; "vheal" is the age-specific set of coefficients 
% for being in bad health; and so on.

% Comment: so first regression is m_coeff, and second regression is sigma_coeff
% Note: the m_coeff is not the mean medical expenses, because of the
% log-expenses formulation the mean expenses is E(m) = exp(mu + 1/2 *sigma^2).

% Following is copy-pasted from contents of medexprof_adj.out
% m_coeff (what DFJ2010 docs refer to as first regression)
medexprof_adj.m_coeff.age=70:1:102;
medexprof_adj.m_coeff.ageshiftFinal=[5.54198301, 5.582088892, 5.609991653, 5.627357155, 5.635732872, 5.636611686, 5.63136711, 5.62131496, 5.607654811, 5.591529413, 5.573962382, 5.555921485, 5.538257273, 5.521761574, 5.507109874, 5.494920062, 5.485692578, 5.479871309, 5.477782758, 5.479694619, 5.485758686, 5.49606959, 5.510602932, 5.529278296, 5.551898496, 5.578208324, 5.607833954, 5.640346602, 5.675199651, 5.711786678, 5.74938749, 5.787221871, 5.824393202];
medexprof_adj.m_coeff.healshift=[-0.126597891, -0.117190485, -0.107783, -0.098375494, -0.088967987, -0.079560502, -0.070152996, -0.06074549, -0.051337984, -0.041930499, -0.032522993, -0.023115486, -0.01370798, -0.004300474, 0.005107011, 0.014514517, 0.023922023, 0.033329508, 0.042737015, 0.052144521, 0.061552006, 0.070959512, 0.080367018, 0.089774503, 0.099182009, 0.108589516, 0.117997001, 0.127404428, 0.136812013, 0.146219435, 0.155626947, 0.165034347, 0.174442059];
medexprof_adj.m_coeff.maleshift=[-0.011614417, -0.021258834, -0.03090325, -0.040547567, -0.050191984, -0.0598363, -0.069480738, -0.079125155, -0.088769371, -0.098413888, -0.108058205, -0.117702642, -0.127346959, -0.136991376, -0.146635692, -0.156280009, -0.165924425, -0.175568763, -0.18521328, -0.194857596, -0.204502013, -0.21414633, -0.223790646, -0.233435084, -0.243079401, -0.252723817, -0.262368134, -0.27201265, -0.281656967, -0.291301405, -0.300945721, -0.310590038, -0.320234455];
medexprof_adj.m_coeff.PIshift=[1.707721202, 1.794494695, 1.881268188, 1.968040892, 2.054814385, 2.141587878, 2.228361583, 2.315135076, 2.401907781, 2.488681274, 2.575454766, 2.662227471, 2.749000964, 2.835774457, 2.922547162, 3.009320655, 3.096093147, 3.182866852, 3.269640345, 3.35641405, 3.443187543, 3.529960036, 3.61673374, 3.703507233, 3.790280726, 3.877053431, 3.963826924, 4.050599629, 4.137373122, 4.224146614, 4.310920319, 4.397692812, 4.484466305];
medexprof_adj.m_coeff.bPI2=-1.537435863*ones(1,33);
% sigma_coeff (what DFJ2010 docs refer to as second regression)
medexprof_adj.sigma_coeff.age=70:1:102;
medexprof_adj.sigma_coeff.ageshiftFinal=[0.859672381, 0.904219617, 0.934751994, 0.95437459, 0.965906155, 0.972021629, 0.97510948, 0.977415081, 0.980896477, 0.987369474, 0.998361836, 1.01525923, 1.039159854, 1.071018951, 1.111505152, 1.161143277, 1.220172244, 1.288686881, 1.366495683, 1.453264763, 1.548370627, 1.651049821, 1.760250136, 1.874778407, 1.993151008, 2.113743353, 2.234646092, 2.353801797, 2.468869698, 2.577366643, 2.67651902, 2.763412258, 2.834845595];
medexprof_adj.sigma_coeff.healshift=[0.641515082, 0.66864537, 0.6957758, 0.722906087, 0.750036375, 0.777166805, 0.804297092, 0.83142738, 0.858557668, 0.885688098, 0.912818385, 0.939948673, 0.96707896, 0.994209248, 1.021339678, 1.048469965, 1.075600253, 1.102730683, 1.129860971, 1.156991258, 1.184121688, 1.211251976, 1.238382263, 1.265512693, 1.292642981, 1.319773269, 1.346903699, 1.374033844, 1.401164274, 1.428295131, 1.455424707, 1.482555706, 1.509685282];
medexprof_adj.sigma_coeff.maleshift=[0.329451635, 0.307158368, 0.284865101, 0.262571834, 0.240278567, 0.217985301, 0.195692176, 0.173398909, 0.151105643, 0.128812376, 0.106519109, 0.084225985, 0.061932718, 0.039639451, 0.017346184, -0.004947083, -0.027240349, -0.049533474, -0.071826741, -0.094120007, -0.116413274, -0.138706541, -0.160999808, -0.183292932, -0.205586199, -0.227879466, -0.250172732, -0.272465999, -0.294759266, -0.31705239, -0.339345657, -0.361638924, -0.383932191];
medexprof_adj.sigma_coeff.PIshift=[1.460772596, 1.51177261, 1.562772625, 1.613771215, 1.664771229, 1.715771244, 1.766769834, 1.817769848, 1.868768439, 1.919768453, 1.970768467, 2.021767058, 2.072767072, 2.123767086, 2.174765677, 2.225765691, 2.276765705, 2.327764295, 2.37876431, 2.4297629, 2.480762914, 2.531762929, 2.582761519, 2.633761533, 2.684761548, 2.735760138, 2.786760152, 2.837758743, 2.888758757, 2.939758771, 2.990757362, 3.041757376, 3.09275739];
medexprof_adj.sigma_coeff.bPI2=-1.384958574*ones(1,33);

% Create Params from these for use in model
% Note, we need two separate Params for health, and can use if-statement in
% ReturnFn to figure out which one applies for current health state
% note: female=0 and male=1 for maleshift (following DFJ2010 regressions)
% note: good health=0 and bad health =1 for healshift (following DFJ2010 regressions)
cpercentile=10/100; % q1
Params.m_coeff_healthgood.femaleq1=medexprof_adj.m_coeff.ageshiftFinal+medexprof_adj.m_coeff.healshift*0+medexprof_adj.m_coeff.maleshift*0+medexprof_adj.m_coeff.PIshift*cpercentile+medexprof_adj.m_coeff.bPI2*cpercentile^2;
Params.m_coeff_healthbad.femaleq1=medexprof_adj.m_coeff.ageshiftFinal+medexprof_adj.m_coeff.healshift*1+medexprof_adj.m_coeff.maleshift*0+medexprof_adj.m_coeff.PIshift*cpercentile+medexprof_adj.m_coeff.bPI2*cpercentile^2;
Params.m_coeff_healthgood.maleq1=medexprof_adj.m_coeff.ageshiftFinal+medexprof_adj.m_coeff.healshift*0+medexprof_adj.m_coeff.maleshift*1+medexprof_adj.m_coeff.PIshift*cpercentile+medexprof_adj.m_coeff.bPI2*cpercentile^2;
Params.m_coeff_healthbad.maleq1=medexprof_adj.m_coeff.ageshiftFinal+medexprof_adj.m_coeff.healshift*1+medexprof_adj.m_coeff.maleshift*1+medexprof_adj.m_coeff.PIshift*cpercentile+medexprof_adj.m_coeff.bPI2*cpercentile^2;
cpercentile=30/100; % q2
Params.m_coeff_healthgood.femaleq2=medexprof_adj.m_coeff.ageshiftFinal+medexprof_adj.m_coeff.healshift*0+medexprof_adj.m_coeff.maleshift*0+medexprof_adj.m_coeff.PIshift*cpercentile+medexprof_adj.m_coeff.bPI2*cpercentile^2;
Params.m_coeff_healthbad.femaleq2=medexprof_adj.m_coeff.ageshiftFinal+medexprof_adj.m_coeff.healshift*1+medexprof_adj.m_coeff.maleshift*0+medexprof_adj.m_coeff.PIshift*cpercentile+medexprof_adj.m_coeff.bPI2*cpercentile^2;
Params.m_coeff_healthgood.maleq2=medexprof_adj.m_coeff.ageshiftFinal+medexprof_adj.m_coeff.healshift*0+medexprof_adj.m_coeff.maleshift*1+medexprof_adj.m_coeff.PIshift*cpercentile+medexprof_adj.m_coeff.bPI2*cpercentile^2;
Params.m_coeff_healthbad.maleq2=medexprof_adj.m_coeff.ageshiftFinal+medexprof_adj.m_coeff.healshift*1+medexprof_adj.m_coeff.maleshift*1+medexprof_adj.m_coeff.PIshift*cpercentile+medexprof_adj.m_coeff.bPI2*cpercentile^2;
cpercentile=50/100; % q3
Params.m_coeff_healthgood.femaleq3=medexprof_adj.m_coeff.ageshiftFinal+medexprof_adj.m_coeff.healshift*0+medexprof_adj.m_coeff.maleshift*0+medexprof_adj.m_coeff.PIshift*cpercentile+medexprof_adj.m_coeff.bPI2*cpercentile^2;
Params.m_coeff_healthbad.femaleq3=medexprof_adj.m_coeff.ageshiftFinal+medexprof_adj.m_coeff.healshift*1+medexprof_adj.m_coeff.maleshift*0+medexprof_adj.m_coeff.PIshift*cpercentile+medexprof_adj.m_coeff.bPI2*cpercentile^2;
Params.m_coeff_healthgood.maleq3=medexprof_adj.m_coeff.ageshiftFinal+medexprof_adj.m_coeff.healshift*0+medexprof_adj.m_coeff.maleshift*1+medexprof_adj.m_coeff.PIshift*cpercentile+medexprof_adj.m_coeff.bPI2*cpercentile^2;
Params.m_coeff_healthbad.maleq3=medexprof_adj.m_coeff.ageshiftFinal+medexprof_adj.m_coeff.healshift*1+medexprof_adj.m_coeff.maleshift*1+medexprof_adj.m_coeff.PIshift*cpercentile+medexprof_adj.m_coeff.bPI2*cpercentile^2;
cpercentile=70/100; % q4
Params.m_coeff_healthgood.femaleq4=medexprof_adj.m_coeff.ageshiftFinal+medexprof_adj.m_coeff.healshift*0+medexprof_adj.m_coeff.maleshift*0+medexprof_adj.m_coeff.PIshift*cpercentile+medexprof_adj.m_coeff.bPI2*cpercentile^2;
Params.m_coeff_healthbad.femaleq4=medexprof_adj.m_coeff.ageshiftFinal+medexprof_adj.m_coeff.healshift*1+medexprof_adj.m_coeff.maleshift*0+medexprof_adj.m_coeff.PIshift*cpercentile+medexprof_adj.m_coeff.bPI2*cpercentile^2;
Params.m_coeff_healthgood.maleq4=medexprof_adj.m_coeff.ageshiftFinal+medexprof_adj.m_coeff.healshift*0+medexprof_adj.m_coeff.maleshift*1+medexprof_adj.m_coeff.PIshift*cpercentile+medexprof_adj.m_coeff.bPI2*cpercentile^2;
Params.m_coeff_healthbad.maleq4=medexprof_adj.m_coeff.ageshiftFinal+medexprof_adj.m_coeff.healshift*1+medexprof_adj.m_coeff.maleshift*1+medexprof_adj.m_coeff.PIshift*cpercentile+medexprof_adj.m_coeff.bPI2*cpercentile^2;
cpercentile=90/100; % q5
Params.m_coeff_healthgood.femaleq5=medexprof_adj.m_coeff.ageshiftFinal+medexprof_adj.m_coeff.healshift*0+medexprof_adj.m_coeff.maleshift*0+medexprof_adj.m_coeff.PIshift*cpercentile+medexprof_adj.m_coeff.bPI2*cpercentile^2;
Params.m_coeff_healthbad.femaleq5=medexprof_adj.m_coeff.ageshiftFinal+medexprof_adj.m_coeff.healshift*1+medexprof_adj.m_coeff.maleshift*0+medexprof_adj.m_coeff.PIshift*cpercentile+medexprof_adj.m_coeff.bPI2*cpercentile^2;
Params.m_coeff_healthgood.maleq5=medexprof_adj.m_coeff.ageshiftFinal+medexprof_adj.m_coeff.healshift*0+medexprof_adj.m_coeff.maleshift*1+medexprof_adj.m_coeff.PIshift*cpercentile+medexprof_adj.m_coeff.bPI2*cpercentile^2;
Params.m_coeff_healthbad.maleq5=medexprof_adj.m_coeff.ageshiftFinal+medexprof_adj.m_coeff.healshift*1+medexprof_adj.m_coeff.maleshift*1+medexprof_adj.m_coeff.PIshift*cpercentile+medexprof_adj.m_coeff.bPI2*cpercentile^2;



cpercentile=10/100; % q1
Params.sigma_coeff_healthgood.femaleq1=medexprof_adj.sigma_coeff.ageshiftFinal+medexprof_adj.sigma_coeff.healshift*0+medexprof_adj.sigma_coeff.maleshift*0+medexprof_adj.sigma_coeff.PIshift*cpercentile+medexprof_adj.sigma_coeff.bPI2*cpercentile^2;
Params.sigma_coeff_healthbad.femaleq1=medexprof_adj.sigma_coeff.ageshiftFinal+medexprof_adj.sigma_coeff.healshift*1+medexprof_adj.sigma_coeff.maleshift*0+medexprof_adj.sigma_coeff.PIshift*cpercentile+medexprof_adj.sigma_coeff.bPI2*cpercentile^2;
Params.sigma_coeff_healthgood.maleq1=medexprof_adj.sigma_coeff.ageshiftFinal+medexprof_adj.sigma_coeff.healshift*0+medexprof_adj.sigma_coeff.maleshift*1+medexprof_adj.sigma_coeff.PIshift*cpercentile+medexprof_adj.sigma_coeff.bPI2*cpercentile^2;
Params.sigma_coeff_healthbad.maleq1=medexprof_adj.sigma_coeff.ageshiftFinal+medexprof_adj.sigma_coeff.healshift*1+medexprof_adj.sigma_coeff.maleshift*1+medexprof_adj.sigma_coeff.PIshift*cpercentile+medexprof_adj.sigma_coeff.bPI2*cpercentile^2;
cpercentile=30/100; % q2
Params.sigma_coeff_healthgood.femaleq2=medexprof_adj.sigma_coeff.ageshiftFinal+medexprof_adj.sigma_coeff.healshift*0+medexprof_adj.sigma_coeff.maleshift*0+medexprof_adj.sigma_coeff.PIshift*cpercentile+medexprof_adj.sigma_coeff.bPI2*cpercentile^2;
Params.sigma_coeff_healthbad.femaleq2=medexprof_adj.sigma_coeff.ageshiftFinal+medexprof_adj.sigma_coeff.healshift*1+medexprof_adj.sigma_coeff.maleshift*0+medexprof_adj.sigma_coeff.PIshift*cpercentile+medexprof_adj.sigma_coeff.bPI2*cpercentile^2;
Params.sigma_coeff_healthgood.maleq2=medexprof_adj.sigma_coeff.ageshiftFinal+medexprof_adj.sigma_coeff.healshift*0+medexprof_adj.sigma_coeff.maleshift*1+medexprof_adj.sigma_coeff.PIshift*cpercentile+medexprof_adj.sigma_coeff.bPI2*cpercentile^2;
Params.sigma_coeff_healthbad.maleq2=medexprof_adj.sigma_coeff.ageshiftFinal+medexprof_adj.sigma_coeff.healshift*1+medexprof_adj.sigma_coeff.maleshift*1+medexprof_adj.sigma_coeff.PIshift*cpercentile+medexprof_adj.sigma_coeff.bPI2*cpercentile^2;
cpercentile=50/100; % q3
Params.sigma_coeff_healthgood.femaleq3=medexprof_adj.sigma_coeff.ageshiftFinal+medexprof_adj.sigma_coeff.healshift*0+medexprof_adj.sigma_coeff.maleshift*0+medexprof_adj.sigma_coeff.PIshift*cpercentile+medexprof_adj.sigma_coeff.bPI2*cpercentile^2;
Params.sigma_coeff_healthbad.femaleq3=medexprof_adj.sigma_coeff.ageshiftFinal+medexprof_adj.sigma_coeff.healshift*1+medexprof_adj.sigma_coeff.maleshift*0+medexprof_adj.sigma_coeff.PIshift*cpercentile+medexprof_adj.sigma_coeff.bPI2*cpercentile^2;
Params.sigma_coeff_healthgood.maleq3=medexprof_adj.sigma_coeff.ageshiftFinal+medexprof_adj.sigma_coeff.healshift*0+medexprof_adj.sigma_coeff.maleshift*1+medexprof_adj.sigma_coeff.PIshift*cpercentile+medexprof_adj.sigma_coeff.bPI2*cpercentile^2;
Params.sigma_coeff_healthbad.maleq3=medexprof_adj.sigma_coeff.ageshiftFinal+medexprof_adj.sigma_coeff.healshift*1+medexprof_adj.sigma_coeff.maleshift*1+medexprof_adj.sigma_coeff.PIshift*cpercentile+medexprof_adj.sigma_coeff.bPI2*cpercentile^2;
cpercentile=70/100; % q4
Params.sigma_coeff_healthgood.femaleq4=medexprof_adj.sigma_coeff.ageshiftFinal+medexprof_adj.sigma_coeff.healshift*0+medexprof_adj.sigma_coeff.maleshift*0+medexprof_adj.sigma_coeff.PIshift*cpercentile+medexprof_adj.sigma_coeff.bPI2*cpercentile^2;
Params.sigma_coeff_healthbad.femaleq4=medexprof_adj.sigma_coeff.ageshiftFinal+medexprof_adj.sigma_coeff.healshift*1+medexprof_adj.sigma_coeff.maleshift*0+medexprof_adj.sigma_coeff.PIshift*cpercentile+medexprof_adj.sigma_coeff.bPI2*cpercentile^2;
Params.sigma_coeff_healthgood.maleq4=medexprof_adj.sigma_coeff.ageshiftFinal+medexprof_adj.sigma_coeff.healshift*0+medexprof_adj.sigma_coeff.maleshift*1+medexprof_adj.sigma_coeff.PIshift*cpercentile+medexprof_adj.sigma_coeff.bPI2*cpercentile^2;
Params.sigma_coeff_healthbad.maleq4=medexprof_adj.sigma_coeff.ageshiftFinal+medexprof_adj.sigma_coeff.healshift*1+medexprof_adj.sigma_coeff.maleshift*1+medexprof_adj.sigma_coeff.PIshift*cpercentile+medexprof_adj.sigma_coeff.bPI2*cpercentile^2;
cpercentile=90/100; % q5
Params.sigma_coeff_healthgood.femaleq5=medexprof_adj.sigma_coeff.ageshiftFinal+medexprof_adj.sigma_coeff.healshift*0+medexprof_adj.sigma_coeff.maleshift*0+medexprof_adj.sigma_coeff.PIshift*cpercentile+medexprof_adj.sigma_coeff.bPI2*cpercentile^2;
Params.sigma_coeff_healthbad.femaleq5=medexprof_adj.sigma_coeff.ageshiftFinal+medexprof_adj.sigma_coeff.healshift*1+medexprof_adj.sigma_coeff.maleshift*0+medexprof_adj.sigma_coeff.PIshift*cpercentile+medexprof_adj.sigma_coeff.bPI2*cpercentile^2;
Params.sigma_coeff_healthgood.maleq5=medexprof_adj.sigma_coeff.ageshiftFinal+medexprof_adj.sigma_coeff.healshift*0+medexprof_adj.sigma_coeff.maleshift*1+medexprof_adj.sigma_coeff.PIshift*cpercentile+medexprof_adj.sigma_coeff.bPI2*cpercentile^2;
Params.sigma_coeff_healthbad.maleq5=medexprof_adj.sigma_coeff.ageshiftFinal+medexprof_adj.sigma_coeff.healshift*1+medexprof_adj.sigma_coeff.maleshift*1+medexprof_adj.sigma_coeff.PIshift*cpercentile+medexprof_adj.sigma_coeff.bPI2*cpercentile^2;


%% Grids
a_grid=2*(10^6)*linspace(0,1,n_a)'.^3; % asset grid, put more points nearer zero where curvature of V is higher
% note: borrowing constraint aprime>=0 is implicitly imposed by grid
% Note: DFJ2010 codes have max 'cash-on-hand' (assets +income) of
% 1.6million, 1million 200,000 depending on type. (lines 75-77 of their wealth20_Singles.c)

h_grid=[0;1;2;3]; % health status, 0=good, 1=bad, 2=death, 3=dead

% AR(1) markov medical expense shock, zeta
[zeta_grid,pi_zeta]=discretizeAR1_FarmerToda(0,Params.rho_m, Params.sigma_mepsilon, n_z(2));

% i.i.d. medical expense shock, xi
[xi_grid,pi_xi]=discretizeAR1_FarmerToda(0,0, Params.sigma_mxi, n_e);
pi_xi=pi_xi(1,:)';

z_grid=[h_grid; zeta_grid];
e_grid=pi_xi;
pi_e=pi_xi;

d_grid=[];

% put e into vfoptions and simoptions
vfoptions.n_e=n_e;
vfoptions.e_grid=e_grid;
vfoptions.pi_e=pi_e;
simoptions.n_e=vfoptions.n_e;
simoptions.e_grid=vfoptions.e_grid;
simoptions.pi_e=vfoptions.pi_e;

% pi_z_J: combine pi_h_J and pi_zeta
% female
pi_z_J.femaleq1=zeros(prod(n_z),prod(n_z),N_j);
for jj=1:N_j
    pi_z_J.femaleq1(:,:,jj)=kron(pi_zeta,pi_h_J.femaleq1(:,:,jj)); % kron() in reverse order
end
pi_z_J.femaleq2=zeros(prod(n_z),prod(n_z),N_j);
for jj=1:N_j
    pi_z_J.femaleq2(:,:,jj)=kron(pi_zeta,pi_h_J.femaleq2(:,:,jj)); % kron() in reverse order
end
pi_z_J.femaleq3=zeros(prod(n_z),prod(n_z),N_j);
for jj=1:N_j
    pi_z_J.femaleq3(:,:,jj)=kron(pi_zeta,pi_h_J.femaleq3(:,:,jj)); % kron() in reverse order
end
pi_z_J.femaleq4=zeros(prod(n_z),prod(n_z),N_j);
for jj=1:N_j
    pi_z_J.femaleq4(:,:,jj)=kron(pi_zeta,pi_h_J.femaleq4(:,:,jj)); % kron() in reverse order
end
pi_z_J.femaleq5=zeros(prod(n_z),prod(n_z),N_j);
for jj=1:N_j
    pi_z_J.femaleq5(:,:,jj)=kron(pi_zeta,pi_h_J.femaleq5(:,:,jj)); % kron() in reverse order
end
% male
pi_z_J.maleq1=zeros(prod(n_z),prod(n_z),N_j);
for jj=1:N_j
    pi_z_J.maleq1(:,:,jj)=kron(pi_zeta,pi_h_J.maleq1(:,:,jj)); % kron() in reverse order
end
pi_z_J.maleq2=zeros(prod(n_z),prod(n_z),N_j);
for jj=1:N_j
    pi_z_J.maleq2(:,:,jj)=kron(pi_zeta,pi_h_J.maleq2(:,:,jj)); % kron() in reverse order
end
pi_z_J.maleq3=zeros(prod(n_z),prod(n_z),N_j);
for jj=1:N_j
    pi_z_J.maleq3(:,:,jj)=kron(pi_zeta,pi_h_J.maleq3(:,:,jj)); % kron() in reverse order
end
pi_z_J.maleq4=zeros(prod(n_z),prod(n_z),N_j);
for jj=1:N_j
    pi_z_J.maleq4(:,:,jj)=kron(pi_zeta,pi_h_J.maleq4(:,:,jj)); % kron() in reverse order
end
pi_z_J.maleq5=zeros(prod(n_z),prod(n_z),N_j);
for jj=1:N_j
    pi_z_J.maleq5(:,:,jj)=kron(pi_zeta,pi_h_J.maleq5(:,:,jj)); % kron() in reverse order
end

%% ReturnFn
DiscountFactorParamNames={'beta'}; 
% Note: conditional survival probabilities are handled by pi_h_J
% Specifically the transtions from h=0 & h=1 (good and bad health) to h=2 (death).

ReturnFn=@(aprime,a,h,zeta,xi,r,upsilon,delta,theta,k,earnings,m_coeff_healthbad,m_coeff_healthgood,sigma_coeff_healthbad,sigma_coeff_healthgood, cfloor, tau_e, estateexemption, taxbracket1, taxbracket2, taxbracket3, taxbracket4, taxbracket5, taxbracket6, margtaxrate0, margtaxrate1, margtaxrate2, margtaxrate3, margtaxrate4, margtaxrate5, margtaxrate6) ...
    DeNardiFrenchJoines2010_ReturnFn(aprime,a,h,zeta,xi,r,upsilon,delta,theta,k,earnings,m_coeff_healthbad,m_coeff_healthgood,sigma_coeff_healthbad,sigma_coeff_healthgood, cfloor, tau_e, estateexemption, taxbracket1, taxbracket2, taxbracket3, taxbracket4, taxbracket5, taxbracket6, margtaxrate0, margtaxrate1, margtaxrate2, margtaxrate3, margtaxrate4, margtaxrate5, margtaxrate6);


%% Solve value fn and policy fn
vfoptions.verbose=1;
vfoptions.verboseparams=1; % show params being sent to each ptype
vfoptions.divideandconquer=1;
% vfoptions.level1n=9; % you could reduce runtimes through playing with this a bit
vfoptions.gridinterplayer=1;
vfoptions.ngridinterp=10;
tic;
[V,Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
vftime=toc
% Note: the level of accuracy being used on assets is very high (501 points, plus 10
% interpolation points between each two consecutive grid points). 
% DFJ2010 appear to have used 110-to-235 points on cash-on-hand (very similar to using 
% them on assets), although their interpolation would not be limited to
% grid.

% When using grid interpolation layer, you have to tell simoptions
simoptions.gridinterplayer=vfoptions.gridinterplayer;
simoptions.ngridinterp=vfoptions.ngridinterp;


%% Quick look a value fn in final period
figure_c=figure_c+1;
figure(figure_c);
subplot(2,1,1); surf(squeeze(V.maleq1(:,:,1,1,end)))
% The warm-glow of bequests is so immense it is the only thing you can see
subplot(2,1,2); surf(squeeze(V.maleq1(:,1:2,1,1,end)))
% Can see that in the first two 'health' states, utility is positive, and then becomes zero when dead

% NOTE: According to eqns in DFJ2010 paper, the warm-glow function phi(e) does not have a '-1' in numerator, while 
% utility fn does. Seems odd given their choice to use same upsilon for both. This left the warm-glow as being very 
% close to zero, so I have put the -1 in here so that there is actually a warm-glow worth keeping assets for.
% [This should be in their code somewhere, but I have not yet tried to find it.]

%% Plot the policy fn, take a look at it
PolicyVals=PolicyInd2Val_Case1_FHorz_PType(Policy,n_d,n_a,n_z,N_j,d_grid,a_grid,simoptions);
% The only policy is aprime
% First, for the middle (third) quintile, look at male and female at different ages
% row: male/female
% column: good health (1st index), bad health (2nd index)
figure_c=figure_c+1;
figure(figure_c);
subplot(2,2,1); plot(a_grid,PolicyVals.maleq3(1,:,1,8,4,1),a_grid,PolicyVals.maleq3(1,:,1,8,4,10),a_grid,PolicyVals.maleq3(1,:,1,8,4,20),a_grid,PolicyVals.maleq3(1,:,1,8,4,30))
legend('age 70', 'age 79', 'age 89', 'age 99')
title('Next period assets: male, good health; mid values for med shocks')
subplot(2,2,2); plot(a_grid,PolicyVals.maleq3(1,:,2,8,4,1),a_grid,PolicyVals.maleq3(1,:,2,8,4,10),a_grid,PolicyVals.maleq3(1,:,2,8,4,20),a_grid,PolicyVals.maleq3(1,:,2,8,4,30))
legend('age 70', 'age 79', 'age 89', 'age 99')
title('Next period assets: male, bad health; mid values for med shocks')
subplot(2,2,3); plot(a_grid,PolicyVals.femaleq3(1,:,1,8,4,1),a_grid,PolicyVals.femaleq3(1,:,1,8,4,10),a_grid,PolicyVals.femaleq3(1,:,1,8,4,20),a_grid,PolicyVals.femaleq3(1,:,1,8,4,30))
legend('age 70', 'age 79', 'age 89', 'age 99')
title('Next period assets: female, good health; mid values for med shocks')
subplot(2,2,4); plot(a_grid,PolicyVals.femaleq3(1,:,2,8,4,1),a_grid,PolicyVals.femaleq3(1,:,2,8,4,10),a_grid,PolicyVals.femaleq3(1,:,2,8,4,20),a_grid,PolicyVals.femaleq3(1,:,2,8,4,30))
legend('age 70', 'age 79', 'age 89', 'age 99')
title('Next period assets: female, bad health; mid values for med shocks')

% DFJ2010 don't appear to report any policy functions, so not sure if these
% are reasonable or not.
% The fact there is a 'cliff' at low assets looks odd, but maybe this is
% the impact of the consumption floor? (I've not solved a consumption floor
% model before; you have to have aprime=0 to qualify for the consumption floor 
% transfers.) I doubt it, cliff is around $1million, and consumption floor is
% just 2653

%% Setup initial agent distribution
% DFJ2010 say this is based on the 1996 data, but we don't have that so
% just making something up for now.
jequaloneDist=struct();

% Based on Figures in DFJ2010, looks like different quintiles start with
% anything from 10000 to 200,000 dollars (median for quintile). So I will
% start them all with roughly this amount based on figure 11.
q1initassets=2000;   [~,q1initassetsind]=min(abs(a_grid-q1initassets)); % value, and then get corresponding index in a_grid (asset grid)
q2initassets=20000;  [~,q2initassetsind]=min(abs(a_grid-q2initassets));
q3initassets=65000;  [~,q3initassetsind]=min(abs(a_grid-q3initassets));
q4initassets=100000; [~,q4initassetsind]=min(abs(a_grid-q4initassets));
q5initassets=800000; [~,q5initassetsind]=min(abs(a_grid-q5initassets));

% For simplicity, start everyone healthy and min medical expenses (mid zeta, first e)
% So put them in (q1initassets,1,8,1)
jequaloneDist.maleq1=zeros([n_a,n_z,n_e],'gpuArray');
jequaloneDist.maleq1(q1initassetsind,1,8,1)=1;
jequaloneDist.maleq2=zeros([n_a,n_z,n_e],'gpuArray');
jequaloneDist.maleq2(q2initassetsind,1,8,1)=1;
jequaloneDist.maleq3=zeros([n_a,n_z,n_e],'gpuArray');
jequaloneDist.maleq3(q3initassetsind,1,8,1)=1;
jequaloneDist.maleq4=zeros([n_a,n_z,n_e],'gpuArray');
jequaloneDist.maleq4(q4initassetsind,1,8,1)=1;
jequaloneDist.maleq5=zeros([n_a,n_z,n_e],'gpuArray');
jequaloneDist.maleq5(q5initassetsind,1,8,1)=1;

jequaloneDist.femaleq1=zeros([n_a,n_z,n_e],'gpuArray');
jequaloneDist.femaleq1(q1initassetsind,1,8,1)=1;
jequaloneDist.femaleq2=zeros([n_a,n_z,n_e],'gpuArray');
jequaloneDist.femaleq2(q2initassetsind,1,8,1)=1;
jequaloneDist.femaleq3=zeros([n_a,n_z,n_e],'gpuArray');
jequaloneDist.femaleq3(q3initassetsind,1,8,1)=1;
jequaloneDist.femaleq4=zeros([n_a,n_z,n_e],'gpuArray');
jequaloneDist.femaleq4(q4initassetsind,1,8,1)=1;
jequaloneDist.femaleq5=zeros([n_a,n_z,n_e],'gpuArray');
jequaloneDist.femaleq5(q5initassetsind,1,8,1)=1;

% Presumably the age weights is based on the 1996 data, I will just put
% equal weights. Note that 'death' is happening within the model.
AgeWeightsParamNames={'mewj'};
Params.mewj=ones(1,N_j)/N_j;

% The weight of each quintile-gender permanent type should also follow the
% 1996 data that I don't have.
PTypeDistParamNames={'ptypemasses'};
Params.ptypemasses=ones(1,10)/10; % 1/5 for each quintile, and 1/2 for each gender


%% Agent distribution
tic;
StationaryDist=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightsParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_z,N_j,Names_i,pi_z_J,Params,simoptions);
disttime=toc

%% Calculate the model statistics that are relevant to the GMM estimation targets
tic;

statstime=toc;

%% Runtimes summary
fprintf(' \n')
fprintf('When doing GMM, you have solve the model a lot of times, as a rule of thumb this has to take at most a few minutes each time \n')
fprintf('So lets look at runtimes: \n')
fprintf('  Value fn took:    %4.1f seconds \n',vftime)
fprintf('  Agent dist took:  %4.1f seconds \n',disttime)
fprintf('  Model stats took: %4.1f seconds \n',statstime)
fprintf('So total runtime was: %4.1f seconds \n', vftime+disttime+statstime)
fprintf(' \n')

%% GMM estimation (NOT YET IMPLEMENTED AS DFJ2010 MATERIALS DO NOT INCLUDE THE TARGET MOMENTS AND THEIR COVAR MATRIX)

EstimParamNames={'upsilon','beta','delta','theta','k','cfloor'};
% Reminder:
% upsilon=CES utility fn param
% beta=discount factor
% delta=disutility of bad health
% theta=bequest intensity
% k=bequest curvature (luxury good)
% cfloor=floor on consumption (from gov transfers)

% "we match median assets by cohort, age, and permanent income"


















