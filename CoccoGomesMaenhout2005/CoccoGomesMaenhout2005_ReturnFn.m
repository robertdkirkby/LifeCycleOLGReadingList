function F=CoccoGomesMaenhout2005_ReturnFn(riskyshare,savings,a,z,e,kappa_ij)
% There are the two decision variables
% riskyshare: is share of savings in risky asset
% savings: is the amount of savings

F=-Inf;

income=exp(kappa_ij+z+e); % labor income
% kappa_ij is the deterministic component
% v is the persistent stochastic component
% epsilon is the transitory stochastic component

c=a+income-savings;

if c>0
    F_c=c; % When using Epstein-Zin preferences we just set the return to c, the role of gamma and psi is handled by the codes
    
    F=F_c; % Bequests, but only with probability of death (1-prob of survival)
end


% Note that both the borrowing constraint and the short-selling contraint
% are actually being enforced by the grids in any case.
% (If savings is forced to be positive, then obviously both the risky and
% riskless assets must be positive as they are just fractions of this)

end