function absW_cstarlstar_minus_W_cstarCEV_Ll0=chooseCEV_L_function(CEV_L, Policy_optimal, Policy_0, Policy_coptimal, Policy_loptimal, Policy_c0, Policy_l0,n_d,n_a,n_z,N_j,pi_z, Params,PTypeDistParamNames,DiscountFactorParamNames)
% See footnote 17 on pg 36 of CKK2009 for explanation of functions being used

N_d=prod(n_d);
N_a=prod(n_a);
N_z=prod(n_z);

% Note that none of the following parameters is age dependent
Model_AltUtilityFn=Params.Model_AltUtilityFn;
if Model_AltUtilityFn==0
    gamma_c=Params.gamma_c;
    gamma=Params.gamma;
else
    gamma_c=Params.gamma_c;
    gamma_l=Params.gamma_l;
    chi=Params.chi;
end

%% W_cstarlstar (Note: V1 should be equal to V_optimal, this provides a double-check that you can do on this code)
c=Policy_coptimal.ptype001; l=Policy_loptimal.ptype001; % Note that l is hours worked (not leisure)
F1.ptype001=(((c.^gamma_c).*((1-l).^(1-gamma_c))).^(1-gamma))/(1-gamma);
if Model_AltUtilityFn==1
    F1.ptype001=(c.^(1-gamma_c))/(1-gamma_c)+chi*((1-l).^(1-gamma_l))/(1-gamma_l);
end
c=Policy_coptimal.ptype002; l=Policy_loptimal.ptype002;
F1.ptype002=(((c.^gamma_c).*((1-l).^(1-gamma_c))).^(1-gamma))/(1-gamma);
if Model_AltUtilityFn==1
    F1.ptype002=(c.^(1-gamma_c))/(1-gamma_c)+chi*((1-l).^(1-gamma_l))/(1-gamma_l);
end

% Calculate the Value Fn by backward iteration
FofPolicy.ptype001=F1.ptype001;
FofPolicy.ptype002=F1.ptype002;

% % The following will also be needed to calculate the expectation of next period value fn, evaluated based on the policy.
% PolicyIndexesKron.ptype001=KronPolicyIndexes_FHorz_Case1(Policy.ptype001, n_d, n_a, n_z,N_j); % Actually with just one d, one a and one z varialbe this is redundant.
% PolicyIndexesKron.ptype002=KronPolicyIndexes_FHorz_Case1(Policy.ptype002, n_d, n_a, n_z,N_j); % Actually with just one d, one a and one z varialbe this is redundant.

V.ptype001=zeros(N_a,N_z,N_j,'gpuArray');
V.ptype002=zeros(N_a,N_z,N_j,'gpuArray');

V.ptype001(:,:,N_j)=FofPolicy.ptype001(:,:,N_j); % Initial guess
V.ptype002(:,:,N_j)=FofPolicy.ptype002(:,:,N_j); % Initial guess

EVnext_zj.ptype001=zeros(N_a,1,'gpuArray');
EVnext_zj.ptype002=zeros(N_a,1,'gpuArray');

for reverse_j=1:N_j-1
    jj=N_j-reverse_j; % current period, counts backwards from J-1
        
    beta=prod(gpuArray(CreateVectorFromParams(Params,DiscountFactorParamNames,jj))); % Both agent types have same beta in CKK2009
    for z_c=1:N_z
        EVnext_zj.ptype001=sum(V.ptype001(:,:,jj+1).*pi_z(z_c,:),2); % size N_z-by-1
        EVnext_zj.ptype002=sum(V.ptype002(:,:,jj+1).*pi_z(z_c,:),2); % size N_z-by-1

        if N_d==0 %length(n_d)==1 && n_d(1)==0
            optaprime.ptype001=Policy_optimal.ptype001(:,z_c,jj);
            optaprime.ptype002=Policy_optimal.ptype002(:,z_c,jj);
        else
            optaprime.ptype001=shiftdim(Policy_optimal.ptype001(2,:,z_c,jj),1);
            optaprime.ptype002=shiftdim(Policy_optimal.ptype002(2,:,z_c,jj),1);
        end
        EVnextOfPolicy_zj.ptype001=EVnext_zj.ptype001(optaprime.ptype001);
        EVnextOfPolicy_zj.ptype002=EVnext_zj.ptype002(optaprime.ptype002);
        
        V.ptype001(:,z_c,jj)=FofPolicy.ptype001(:,z_c,jj)+beta*EVnextOfPolicy_zj.ptype001;
        V.ptype002(:,z_c,jj)=FofPolicy.ptype002(:,z_c,jj)+beta*EVnextOfPolicy_zj.ptype002;
    end
end

%Transforming Value Fn out of Kronecker Form
V1.ptype001=reshape(V.ptype001,[n_a,n_z,N_j]);
V1.ptype002=reshape(V.ptype002,[n_a,n_z,N_j]);

W_cstarlstar=sum(Params.(PTypeDistParamNames{1}).*[V1.ptype001(1,1+floor(n_z/2),1); V1.ptype002(1,1+floor(n_z/2),1)]); % (1,4,1,:) is no assets (1), average labor productivity (1+floor(n_z/2)), age one (1), both permanent types (:).


%% W_cstarCEV_Ll0
c=(1+CEV_L)*Policy_coptimal.ptype001; l=Policy_l0.ptype001; % Note that l is hours worked (not leisure)
F2.ptype001=(((c.^gamma_c).*((1-l).^(1-gamma_c))).^(1-gamma))/(1-gamma);
if Model_AltUtilityFn==1
    F2.ptype001=(c.^(1-gamma_c))/(1-gamma_c)+chi*((1-l).^(1-gamma_l))/(1-gamma_l);
end
c=(1+CEV_L)*Policy_coptimal.ptype002; l=Policy_l0.ptype002;
F2.ptype002=(((c.^gamma_c).*((1-l).^(1-gamma_c))).^(1-gamma))/(1-gamma);
if Model_AltUtilityFn==1
    F2.ptype002=(c.^(1-gamma_c))/(1-gamma_c)+chi*((1-l).^(1-gamma_l))/(1-gamma_l);
end

% Calculate the Value Fn by backward iteration
FofPolicy.ptype001=F2.ptype001;
FofPolicy.ptype002=F2.ptype002;

% % The following will also be needed to calculate the expectation of next period value fn, evaluated based on the policy.
% PolicyIndexesKron.ptype001=KronPolicyIndexes_FHorz_Case1(Policy.ptype001, n_d, n_a, n_z,N_j); % Actually with just one d, one a and one z varialbe this is redundant.
% PolicyIndexesKron.ptype002=KronPolicyIndexes_FHorz_Case1(Policy.ptype002, n_d, n_a, n_z,N_j); % Actually with just one d, one a and one z varialbe this is redundant.

V.ptype001=zeros(N_a,N_z,N_j,'gpuArray');
V.ptype002=zeros(N_a,N_z,N_j,'gpuArray');


V.ptype001(:,:,N_j)=FofPolicy.ptype001(:,:,N_j); % Initial guess
V.ptype002(:,:,N_j)=FofPolicy.ptype002(:,:,N_j); % Initial guess


EVnext_zj.ptype001=zeros(N_a,1,'gpuArray');
EVnext_zj.ptype002=zeros(N_a,1,'gpuArray');

for reverse_j=1:N_j-1
    jj=N_j-reverse_j; % current period, counts backwards from J-1

    beta=prod(gpuArray(CreateVectorFromParams(Params,DiscountFactorParamNames,jj))); % Both agent types have same beta in CKK2009
    for z_c=1:N_z
        EVnext_zj.ptype001=sum(V.ptype001(:,:,jj+1).*pi_z(z_c,:),2); % size N_z-by-1
        EVnext_zj.ptype002=sum(V.ptype002(:,:,jj+1).*pi_z(z_c,:),2); % size N_z-by-1
        
        if N_d==0 %length(n_d)==1 && n_d(1)==0
            optaprime.ptype001=Policy_0.ptype001(:,z_c,jj); % I HAVE USED Policy_0 HERE FOR THE EXPECTATIONS (next period assets). NOT CLEAR THIS IS CORRECT.
            optaprime.ptype002=Policy_0.ptype002(:,z_c,jj);
        else
            optaprime.ptype001=shiftdim(Policy_0.ptype001(2,:,z_c,jj),1);
            optaprime.ptype002=shiftdim(Policy_0.ptype002(2,:,z_c,jj),1);
        end
        EVnextOfPolicy_zj.ptype001=EVnext_zj.ptype001(optaprime.ptype001);
        EVnextOfPolicy_zj.ptype002=EVnext_zj.ptype002(optaprime.ptype002);
        
        V.ptype001(:,z_c,jj)=FofPolicy.ptype001(:,z_c,jj)+beta*EVnextOfPolicy_zj.ptype001;
        V.ptype002(:,z_c,jj)=FofPolicy.ptype002(:,z_c,jj)+beta*EVnextOfPolicy_zj.ptype002;
    end
end

%Transforming Value Fn out of Kronecker Form
V2.ptype001=reshape(V.ptype001,[n_a,n_z,N_j]);
V2.ptype002=reshape(V.ptype002,[n_a,n_z,N_j]);

W_copitmalCEV_Ll0=sum(Params.(PTypeDistParamNames{1}).*[V2.ptype001(1,1+floor(n_z/2),1); V2.ptype002(1,1+floor(n_z/2),1)]); % (1,4,1,:) is no assets (1), average labor productivity (1+floor(n_z/2)), age one (1), both permanent types (:).


%% 
% absW_cstarlstar_minus_W_cstarCEV_Ll0=(W_cstarlstar-W_copitmalCEV_Ll0)^2; % This has the problem that everything is already of the order of around 10^(-4) so we have reached to optimization tolerance before we even begin
absW_cstarlstar_minus_W_cstarCEV_Ll0=((W_cstarlstar-W_copitmalCEV_Ll0)^2)/abs(W_cstarlstar); % Divide by absolute value (of the one that doesn't depend on CEV_L) so that it is about relative magnitude

absW_cstarlstar_minus_W_cstarCEV_Ll0=absW_cstarlstar_minus_W_cstarCEV_Ll0*10^4; % this is just to make it larger (was added as otherwise is anyway zero to four decimal places)

fprintf('Current iteration of choose CEV_L: CEV_L=%8.4f, W_cstarlstar=%8.8f, W_copitmalCEV_Ll0=%8.8f, target to zero=%8.8f \n', CEV_L, W_cstarlstar, W_copitmalCEV_Ll0, absW_cstarlstar_minus_W_cstarCEV_Ll0)

end