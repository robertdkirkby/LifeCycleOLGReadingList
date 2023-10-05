function laborincometax=ConesaKitaoKrueger2009_LaborIncomeTaxRevenueFn(l,kprime,k,eta,r,tau_ss,ybar,epsilon_j,alpha_i,I_j,A,alpha,delta, tau_kappa0, tau_kappa1, tau_kappa2) %i

w=(1-alpha)*(A^(1/(1-alpha)))*((r+delta)/alpha)^(alpha/(alpha-1)); % eqn 10 of CKK2009, combined with eqn 9
% Yhat=(rhat+delta)/(1-alpha);

y_pretaxlabor=w*epsilon_j*eta*alpha_i*l*I_j;
socialsecuritytax=0.5*tau_ss*min(y_pretaxlabor,ybar);
y_taxablelabor=y_pretaxlabor-socialsecuritytax;
laborincometax=0;
if y_taxablelabor>0
    % Incorporates the extension of Fehr & Kinderman (2015) to allow tau_kappa1=0 or infinity
    if tau_kappa1==0
        laborincometax=tau_kappa0*y_taxablelabor;
    elseif tau_kappa1==Inf
        laborincometax=tau_kappa0*max(y_taxablelabor-tau_kappa2,0);
    else
        laborincometax=tau_kappa0*(y_taxablelabor - (y_taxablelabor^(-tau_kappa1) + tau_kappa2)^(-1/tau_kappa1));
    end
end
% y_capitalincome=r*(k+Tr_beq);
% capitalincometax=tau_k*y_capitalincome;
% 
% socialsecuritybenefit=SS*(1-I_j);
% 
% c_pretax=(1+r)*(k+Tr_beq)-capitalincometax+y_pretaxlabor-socialsecuritytax-laborincometax+socialsecuritybenefit-kprime;
% c=c_pretax/(1-tau_c);
% 
% F=-Inf;
% 
% if c>0 && l<1
%     F=(((c^gamma_c) *((1-l)^(1-gamma_c))).^(1-gamma))/(1-gamma);
% end
% 
% if Model_AltUtilityFn==1
%     if c>0 && l<1
%         F=(c^(1-gamma_c))/(1-gamma_c)+chi*((1-l)^(1-gamma_l))/(1-gamma_l);
%     end
% end

    
end