function F=BorellaDeNardiYang2018_ReturnFn(n_f,n_m,aprime,a,z_f,z_m,htype,r,omega,gamma,tau_ss,participationcost_men,participationcost_women_single,participationcost_women_married,kappa_j_women,kappa_j_men,potentialhours)
% This is the Economy 4 model. 
% htype==1 is married couple, htype==2 is single female, htype==3 is single male, and hype==4 is deceased

% Note: all the equations that make up the return function are contained on page 11 of BDNY2018
% Note: the requirement that retirees work 0 hours is enforced via kappa_j
% (the age-dependent determinisitic labor productivity; denoted e^i_t by
% BDNY2018) taking value of zero in retirement (so since working has a
% utility cost, and will not generate earnings, everyone chooses zero labour supply)

F=-Inf;

if htype==1 % married couple

    l_f=1-n_f-participationcost_women_married*(n_f>0);
    l_m=1-n_m-participationcost_men*(n_m>0);
    earnings_f=kappa_j_women*z_f*potentialhours*n_f; % BDNY2018 refer to these as e, epsilon, n
    earnings_m=kappa_j_men*z_m*potentialhours*n_m; % BDNY2018 refer to these as e, epsilon, n
    c=(1+r)*a+(1-tau_ss)*(earnings_f+earnings_m)-aprime;
    % Note: aprime>=0 and n_f>=0 are both implicit in the grids
    if c>0 && l_f>0 && l_m>0
        F=((((c/2)^omega)*(l_f^(1-omega)))^(1-gamma)-1)/(1-gamma)+((((c/2)^omega)*(l_m^(1-omega)))^(1-gamma)-1)/(1-gamma);
    end

elseif htype==2 % single female


    l_f=1-n_f-participationcost_women_single*(n_f>0);
    earnings_f=kappa_j_women*z_f*potentialhours*n_f; % BDNY2018 refer to these as e, epsilon, n
    c=(1+r)*a+(1-tau_ss)*earnings_f-aprime;
    if c>0 && l_f>0
        F=(((c^omega)*(l_f^(1-omega)))^(1-gamma)-1)/(1-gamma);
    end     

    if n_m>0 % not important to anything, just felt like enforcing this
        F=-Inf;
    end
elseif htype==3 % single male

    l_m=1-n_m-participationcost_men*(n_m>0);
    earnings_m=kappa_j_men*z_m*potentialhours*n_m; % BDNY2018 refer to these as e, epsilon, n
    c=(1+r)*a+(1-tau_ss)*earnings_m-aprime;
    if c>0 && l_m>0
        F=(((c^omega)*(l_m^(1-omega)))^(1-gamma)-1)/(1-gamma);
    end     

    if n_f>0 % not important to anything, just felt like enforcing this
        F=-Inf;
    end
elseif htype==4 % deceased

    F=0;
    
    if n_f>0 || n_m>0 % not important to anything, just felt like enforcing this
        F=-Inf;
    end
end





end