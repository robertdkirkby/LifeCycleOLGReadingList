% Calculate the MPC in Carroll (1997)
% MPC=marginal propensity to consume
% Note that you have to define what this is the MPC out of?
% We will calculate the MPC out of permanent income
%   And the MPC out of transitory income
%   And the MPC out of wealth (so assets)
% We will calculate (each of these) MPC at the median income, but conditional on each different level of assets (and on age)
% We graph this in Figure 1
% We will then calculate the average across assets (based on agent distribution), conditional on age
% We graph this in Figure 2


% Rough definition of MPC=Delta C/Delta X.
%    Delta=change in
%    C=consumption
%    X=what we calculate MPC out of (so permanent income, wealth, etc.)

% So what we will need to calculate the MPC is to calculate the consumption
% for some given X, then calculate consumption again for a different X,
% then take the difference (giving us Delta C) and finally divide by Delta X.

% First, just run a code to solve the model of Carroll (1997):
Carroll1997
% [from: https://github.com/robertdkirkby/LifeCycleOLGReadingList/tree/main/Carroll1997]


%% Let's start with the MPC out of wealth. This will be nice as we already have an asset grid which is essentially wealth.
% So we just need to look at how much consumptions changes between different asset (wealth) levels.
% One of the things the Carroll1997 code creates is ValuesOnGrid.consumption, which gives us the values of consumption at all the grid points.
% So we can just compute the MPC out of wealth as
% DeltaC=ValuesOnGrid.consumption(2:end,:,:,:)-ValuesOnGrid.consumption(1:end-1,:,:,:); % I'm just going to ignore the top grid point
% DeltaW=a_grid(2:end)-a_grid(1:end-1);
% Decided to use +4 grid points (as otherwise gets a bit messing because increase in wealth from moving up one grid point is so tiny it iteracts with the small errors coming from forcing choices to be on grid)
DeltaC=ValuesOnGrid.consumption(5:end,:,:,:)-ValuesOnGrid.consumption(1:end-4,:,:,:); % I'm just going to ignore the top grid point
DeltaW=a_grid(5:end)-a_grid(1:end-4);


MPCwealth=DeltaC./DeltaW;
% Note the currently it is a different MPC at every point on the grid
size(MPCwealth)
% We can plot this conditional on age, at the median income (shocks)
figure(1)
subplot(3,1,1); plot(a_grid(1:end-4),MPCwealth(:,ceil(n_z/2),7,1))
hold on
subplot(3,1,1); plot(a_grid(1:end-4),MPCwealth(:,ceil(n_z/2),7,11))
subplot(3,1,1); plot(a_grid(1:end-4),MPCwealth(:,ceil(n_z/2),7,21))
hold off
title('MPC out of wealth (as fn of assets, at median income/shocks)')
legend('j=1','j=11','j=21')


% Alternatively, we might want some average MPC, so we can take MPC over agent distribution (keeping just the age dimension)
AvgMPCwealth_age=MPCwealth.*StationaryDist(1:end-4,:,:,:); % just hoping top four grid point of assets is anyway zero mass (you should check this)
AvgMPCwealth_age=squeeze(sum(sum(sum(AvgMPCwealth_age,1),2),3));
figure(2)
subplot(3,1,1); plot(1:1:N_j,AvgMPCwealth_age)
title('Average MPC out of wealth')

%% Next, MPC out of permanent income, but only at median permanent income shock
DeltaC=ValuesOnGrid.consumption(:,ceil(n_z/2)+1,:,:)-ValuesOnGrid.consumption(:,ceil(n_z/2),:,:);
DeltaPermanentIncome=reshape(P_grid_J(ceil(n_z/2)+1,:)-P_grid_J(ceil(n_z/2),:),[1,1,1,N_j]);
MPCpermanentincome=DeltaC./DeltaPermanentIncome;

size(MPCpermanentincome) % note the singular second dimension as just at median permanent shock

% We can plot this conditional on age, at the median income (shocks)
figure(1)
subplot(3,1,2); plot(a_grid,MPCpermanentincome(:,1,7,1))
hold on
subplot(3,1,2); plot(a_grid,MPCpermanentincome(:,1,7,11))
subplot(3,1,2); plot(a_grid,MPCpermanentincome(:,1,7,21))
hold off
title('MPC out of permanent income (as fn of assets, at median income/shocks)')
legend('j=1','j=11','j=21')


% Alternatively, we might want some average MPC, so we can take MPC over agent distribution (keeping just the age dimension)
AvgMPCpermanentincome_age=MPCpermanentincome.*(StationaryDist(:,ceil(n_z/2),:,:)/sum(sum(sum(sum(StationaryDist(:,ceil(n_z/2),:,:))))));
AvgMPCpermanentincome_age=squeeze(sum(sum(sum(AvgMPCpermanentincome_age,1),2),3));
figure(2)
subplot(3,1,2); plot(1:1:N_j,AvgMPCpermanentincome_age)
title('Average MPC out of permanent income)')


%% Lastly, MPC out of transitory income, but only at the median transitory income shock
DeltaC=ValuesOnGrid.consumption(:,:,8,:)-ValuesOnGrid.consumption(:,:,7,:);
DeltaTransitoryIncome=vfoptions.e_grid(8)-vfoptions.e_grid(7);
MPCtransitoryincome=DeltaC./DeltaTransitoryIncome;

size(MPCtransitoryincome) % note the singular second dimension as just at median permanent shock

% We can plot this conditional on age, at the median income (shocks)
figure(1)
subplot(3,1,3); plot(a_grid,MPCtransitoryincome(:,ceil(n_z/2),1,1))
hold on
subplot(3,1,3); plot(a_grid,MPCtransitoryincome(:,ceil(n_z/2),1,11))
subplot(3,1,3); plot(a_grid,MPCtransitoryincome(:,ceil(n_z/2),1,21))
hold off
title('MPC out of transitory income (as fn of assets, at median income/shocks)')
legend('j=1','j=11','j=21')


% Alternatively, we might want some average MPC, so we can take MPC over agent distribution (keeping just the age dimension)
AvgMPCtransitoryincome_age=MPCtransitoryincome.*(StationaryDist(:,:,7,:)/sum(sum(sum(sum(StationaryDist(:,:,7,:))))));
AvgMPCtransitoryincome_age=squeeze(sum(sum(sum(AvgMPCtransitoryincome_age,1),2),3));
figure(2)
subplot(3,1,3); plot(1:1:N_j,AvgMPCtransitoryincome_age)
title('Average MPC out of transitory income)')



%% Remarks
% The average MPC out of wealth and transitory income are 'lower' than out
% of permanent income, because the former two are just a one off increase
% so you want to spread consumption of them across your lifetime

% The MPC out of wealth/transitory income in the first few periods is high
% as you are borrowing constrained and this extra income helps loosen the
% borrowing constraint.

% The MPc out of wealth/transitory income goes up in last few periods as
% you are smoothing the consumption of that extra income over just the few
% remaining periods (there is no warm glow bequests in this model)


%% Final comments
% The smaller the DeltaX used to compute MPC=DeltaC/DeltaX, the better (as
% in principle in the limit as DeltaX goes to zero we are computing dC/dX
% which is the more theoretically based definition of MPC)

% I compute all of this with an increase in X (so DeltaX is positive). We
% might be interested in, e.g., is the MPC out of housing wealth
% asymmetric? In this case we want to compute MPC for both positive and negative
% DeltaX and compare them.

% Two things we might be interested in: (i) getting a theoretically clean MPC
% out of the model, (ii) getting an MPC out of the model that is calcuated
% in the same way we calculate the MPC in the data. Here I do something
% closer to (i), if you something more like (ii) you might want to
% simulate panel data and compute your MPC from that in the same way you do
% from the real-world data.

% MPCs can require you to use finer grids than would be necessary for many other model statistics.
% You can see the grids are largely irrelevant for the 'average MPC' but
% mess up some of those for each grid point (by a small magnitude each time)


