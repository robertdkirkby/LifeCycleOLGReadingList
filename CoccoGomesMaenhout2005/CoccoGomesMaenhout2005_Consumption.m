function c=CoccoGomesMaenhout2005_Consumption(riskyshare,savings,a,z,e,kappa_ij)
% d1 is share of savings in risky asset
% d2 is the amount of savings

income=exp(kappa_ij+z+e); % labor income
% kappa_ij is the deterministic component
% v is the persistent stochastic component
% epsilon is the transitory stochastic component

c=a+income-savings;

end