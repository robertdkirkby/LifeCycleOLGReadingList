function F=AttanasioLowSanchezMarcos2008_ReturnFn(P,aprime,a,h_f,z_m,z_f,childcarecost,e,h_m,gamma,phi1,phi2,y_f0,y_m0,R)

F=-Inf;

y_f=y_f0*h_f*exp(z_f); % female earnings (if participates)
y_m=y_m0*h_m*exp(z_m); % male earnings

c=a+(y_f-childcarecost)*P+y_m-aprime/R; % Note: the timing of R is just used by original because they use the cash-on-hand setup, here it is pointless

if c>0
    U_c=(((c/e)^(1-gamma))/(1-gamma));
    F=U_c*exp(phi1*P)-phi2*P;
end

% To test, I added the following (and added agej as input). Gives results
% you would expect (assets drift ever up with time, and participation
% zig-zag up and down)
% if (aprime-a)<0.01
%     F=-Inf;
% end
% 
% if rem(agej,2)==P
%     F=-Inf;
% end

end