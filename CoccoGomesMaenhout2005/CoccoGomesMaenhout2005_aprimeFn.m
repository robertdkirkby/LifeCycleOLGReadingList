function aprime=CoccoGomesMaenhout2005_aprimeFn(riskyshare,savings,u, Rf)
% Note: because of how Case3 works we need to input (d1,d2,u) as the first arguements.
% That is, the first inputs must be the decision variables (d variables),
% followed by the shocks that are iid and occur between periods (u variables)

Rr=Rf+u; % Note Rf is the gross interest rate (1.02 in calibration)

aprime=Rf*(1-riskyshare)*savings+Rr*riskyshare*savings;

end