function F=Carroll1997_ReturnFn(aprime,a,P,V,r,gamma,agej,Jr,pension)

F=-Inf;

if agej<Jr
    y=V*P;
else
    y=pension*P; %note, because of pi_P_J, P will be constant in retirement.
end
c=a+y-aprime/(1+r); % Normally I would have a*(1+r), but want to follow precise calibration of Carroll (1997)

if c>0
    F=(c^(1-gamma))/(1-gamma);
end

end