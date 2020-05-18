% Matlab script to compare different models when fitted to data for all
% serotypes.
%
% This implements the PSIS-LOO method described in Vehtari et al.
% Stat Comput 27, 1413-1432 (2017)

% Specify the models to be considered
mList=1:10;

%==========================================================================
% DATA
% Load the challenge data
M=load('BuffaloTransmissionData.txt');

% Extract the challange data. Vectors are:
% sType - serotype (1-3 for SAT1-SAT3)
% iType - infection type (1-inoculated, 2-contact)
% tVI - times of last negative VI result, first positive VI result, last
%       positive VI result, first negative VI result
sType=M(:,1);
iType=M(:,2);
tVI=M(:,3:6);
clear('M');

% Set the time at which contact buffalo are introduced to the inoculated
% ones
tIntro=2;

% Compute the number of animals
nA=length(iType);
%==========================================================================

%==========================================================================
% COMPUTE THE ELPDs AND WAICs
% Create vectors to store the expected log predictive densities and the
% effective number of parameters for the models
elpdi=zeros(length(mList),nA);
elpd=zeros(length(mList),1);
pLOO=zeros(length(mList),1);

% Create arrays to store WAIC and the effective number of parameters for
% the models
WAIC=zeros(length(mList),1);
pWAIC=zeros(length(mList),1);

% For each model ...
for m=1:length(mList)
    disp(['model ' num2str(mList(m))])

% Load the MCMC samples
    varload=load(['Model' num2str(mList(m)) '_MCMCSamples']);
    ParSamp=varload.ParSamp;

% Extract the parameters (including unobserved infection times, etc.) and
% the log likelihood
    PS=[];
    logL=[];
    for chain=1:length(ParSamp)
        PS=[PS; ParSamp{chain}(:,1:end-2)];
        logL=[logL; ParSamp{chain}(:,end)];
    end
    
% Determine the number of samples
    S=size(PS,1);
    M=0.2*S;
    
% Compute the log pointwise predictive density
    disp('   computing lppd')
    lppd=zeros(size(PS,1),nA);
    for j=1:S
        lppd(j,:)=PointwiseLhood_AcuteTransmission(PS(j,:)',mList(m),...
                                                   tVI,iType,sType,tIntro)';
        if abs(logL(j)-sum(lppd(j,:)))>1e-12
            disp(' warning: something wrong with likelihood')
        end
    end

% Compute the raw importance weights
    disp('   computing importance weights')
    w=1./exp(lppd);

% For each animal ...
    for a=1:nA
        
% Fit a generalised Pareto distribution to the largest 20% of the weights
        [wSort,ind]=sort(w(:,a));
        pHat=fminsearch(@(par)(gplike(par,wSort(S-M+1:S))),[0 1]);
        if pHat(1)>0.5
            display(['     tail parameter > 0.5 for animal ' num2str(a)])
        end

% Stabilize the importance weights by replacing the largest 20% with the
% Pareto-smoothed values
        z=((1:M)-0.5)/M;
        wNew=gpinv(z,pHat(1),pHat(2));
        w(ind(S-M+1:S),a)=wNew;

% Truncate the smoothed weights
        wMax=(S^0.75)*mean(w(:,a));
        w(w(:,a)>=wMax,a)=wMax;

    end

% Compute the PSIS-LOO estimate of the ELPD
    disp('   computing elpd')
    elpdi(m,:)=-2*log(sum(w.*exp(lppd),1)./sum(w,1));
    elpd(m)=sum(elpdi(m,:));

% Compute the effective number of parameters (based on LOO)
    pLOO(m)=sum(log(mean(exp(lppd),1)))-(-0.5*sum(elpdi(m,:)));

end
%==========================================================================

% Save the results
save('ModelComparisonToo_PSIS-LOO','elpdi','elpd','pLOO','mList')

% Tidy up
clear
