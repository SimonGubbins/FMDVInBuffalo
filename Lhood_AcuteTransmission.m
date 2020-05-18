function [logL, prior]=Lhood_AcuteTransmission(par,mFlag,tVI,iType,...
                                               sType,tIntro)
%
% [logL, prior]=Lhood_AcuteTransmission(par,mFlag,tVI,iType,...
%                                       sType,tIntro)
%
% Matlab function for computing the log likelihood and prior for
% transmission, latent and infectious period parameters for FMDV in buffalo
% based on challenge experiments
%
% This version includes a range of models (see below) to allow for
% differences amongst serotypes. It also imputes the infectious periods.
%
% Inputs:
% par - array containing the (transformed) model parameter
% mFlag - flag indicating model to use:
%         1-SIR, common parameters for all serotypes
%         2-SIR, common infectious period, different transmission rate
%         3-SIR, different infectious period, common transmission rate
%         4-SIR, different parameters for all serotypes
%         5-SEIR, common parameters for all serotypes
%         6-SEIR, common periods, different transmission rate
%         7-SEIR, different latent period, common infectious period and
%           transmission rate
%         8-SEIR, different infectious period, common latent period and
%           transmission rate
%         9-SEIR, different periods, common transmission rate
%         10-SEIR, different parameters for all serotypes
% tVI - array containing times of last negative VI result, first positive
%       VI result, last positive VI result, first negative VI result
% iType - vector identifying infection type (1-inoculated, 2-contact)
% sType - vector identifying serotype (1-3 for SAT1-3)
% tIntro - time at which contact animals introduced to inoculated ones
%          (assumed to be the same for all contacts)
%
% Outputs:
% logL - log likelihood for parameters
% prior - log prior probability for parameters

%==========================================================================
% PREPARE THE INPUTS
% Set the number of animals
nA=length(iType);

% Identify and compute the number of contact animals
c=(iType==2);
nC=length(find(c));

% Compute the number of serotypes
nS=max(sType);

% SIR, common parameters
if mFlag==1
    kI=par(nA+nC+1).*ones(nS,1);
    muI=par(nA+nC+2).*ones(nS,1);
    b=par(nA+nC+3).*ones(nS,1);

% SIR, common periods, different transmission rates
elseif mFlag==2
    kI=par(nA+nC+1).*ones(nS,1);
    muI=par(nA+nC+2).*ones(nS,1);
    b=par(nA+nC+3:nA+nC+nS+2);
    mb=par(nA+nC+nS+3);

% SIR, different periods, common transmission rates
elseif mFlag==3
    kI=par(nA+nC+1:nA+nC+nS);
    muI=par(nA+nC+nS+1:nA+nC+2*nS);
    b=par(nA+nC+2*nS+1).*ones(nS,1);
    mkI=par(nA+nC+2*nS+2);
    mmuI=par(nA+nC+2*nS+3);

% SIR, different periods and transmission rates
elseif mFlag==4
    kI=par(nA+nC+1:nA+nC+nS);
    muI=par(nA+nC+nS+1:nA+nC+2*nS);
    b=par(nA+nC+2*nS+1:nA+nC+3*nS);
    mkI=par(nA+nC+3*nS+1);
    mmuI=par(nA+nC+3*nS+2);
    mb=par(nA+nC+3*nS+3);

% SEIR, common parameters
elseif mFlag==5
    kE=par(nA+2*nC+1).*ones(nS,1);
    muE=par(nA+2*nC+2).*ones(nS,1);
    kI=par(nA+2*nC+3).*ones(nS,1);
    muI=par(nA+2*nC+4).*ones(nS,1);
    b=par(nA+2*nC+5).*ones(nS,1);

% SEIR, common periods, different transmission rates
elseif mFlag==6
    kE=par(nA+2*nC+1).*ones(nS,1);
    muE=par(nA+2*nC+2).*ones(nS,1);
    kI=par(nA+2*nC+3).*ones(nS,1);
    muI=par(nA+2*nC+4).*ones(nS,1);
    b=par(nA+2*nC+5:nA+2*nC+nS+4);
    mb=par(nA+2*nC+nS+5);

% SEIR, different latent periods, common infectious periods and
% transmission rates
elseif mFlag==7
    kE=par(nA+2*nC+1:nA+2*nC+nS);
    muE=par(nA+2*nC+nS+1:nA+2*nC+2*nS);
    kI=par(nA+2*nC+2*nS+1).*ones(nS,1);
    muI=par(nA+2*nC+2*nS+2).*ones(nS,1);
    b=par(nA+2*nC+2*nS+3).*ones(nS,1);
    mkE=par(nA+2*nC+2*nS+4);
    mmuE=par(nA+2*nC+2*nS+5);

% SEIR, different infectious periods, common latent periods and
% transmission rates
elseif mFlag==8
    kE=par(nA+2*nC+1).*ones(nS,1);
    muE=par(nA+2*nC+2).*ones(nS,1);
    kI=par(nA+2*nC+3:nA+2*nC+nS+2);
    muI=par(nA+2*nC+nS+3:nA+2*nC+2*nS+2);
    b=par(nA+2*nC+2*nS+3).*ones(nS,1);
    mkI=par(nA+2*nC+2*nS+4);
    mmuI=par(nA+2*nC+2*nS+5);

% SEIR, different periods, common transmission rates
elseif mFlag==9
    kE=par(nA+2*nC+1:nA+2*nC+nS);
    muE=par(nA+2*nC+nS+1:nA+2*nC+2*nS);
    kI=par(nA+2*nC+2*nS+1:nA+2*nC+3*nS);
    muI=par(nA+2*nC+3*nS+1:nA+2*nC+4*nS);
    b=par(nA+2*nC+4*nS+1).*ones(nS,1);
    mkE=par(nA+2*nC+4*nS+2);
    mmuE=par(nA+2*nC+4*nS+3);
    mkI=par(nA+2*nC+4*nS+4);
    mmuI=par(nA+2*nC+4*nS+5);

% SEIR, different periods and transmission rates
elseif mFlag==10
    kE=par(nA+2*nC+1:nA+2*nC+nS);
    muE=par(nA+2*nC+nS+1:nA+2*nC+2*nS);
    kI=par(nA+2*nC+2*nS+1:nA+2*nC+3*nS);
    muI=par(nA+2*nC+3*nS+1:nA+2*nC+4*nS);
    b=par(nA+2*nC+4*nS+1:nA+2*nC+5*nS);
    mkE=par(nA+2*nC+5*nS+1);
    mmuE=par(nA+2*nC+5*nS+2);
    mkI=par(nA+2*nC+5*nS+3);
    mmuI=par(nA+2*nC+5*nS+4);
    mb=par(nA+2*nC+5*nS+5);

end

% Set the infection times for the animals (inoculated animals are assumed
% to be infected at t=0)
tI=zeros(size(iType));
tI(c)=par(1:nC);

% Set the latent and infectious periods for the animals
E=zeros(size(iType));
if mFlag<=4
    I=par(nC+1:nC+nA);
elseif mFlag>=5
    E(c)=par(nC+1:2*nC);
    I=par(2*nC+1:2*nC+nA);
end
%==========================================================================

%==========================================================================
% COMPUTE THE PRIOR
% Calculate the log(prior) for the parameters
if mFlag==1
    prior=log(gampdf(kI(1),1,1))+...
          log(gampdf(muI(1),2,6/2))+...
          log(gampdf(b(1),1.3942,1.3304));
elseif mFlag==2
    prior=log(gampdf(kI(1),1,1))+...
          log(gampdf(muI(1),2,6/2))+...
          sum(log(exppdf(b,mb)))+...
          log(gampdf(mb,1.3942,1.3304));
elseif mFlag==3
    prior=sum(log(exppdf(kI,mkI)))+...
          sum(log(exppdf(muI,mmuI)))+...
          log(gampdf(b(1),1.3942,1.3304))+...
          log(exppdf(mkI,1))+...
          log(gampdf(mmuI,2,6/2));
elseif mFlag==4
    prior=sum(log(exppdf(kI,mkI)))+...
          sum(log(exppdf(muI,mmuI)))+...
          sum(log(exppdf(b,mb)))+...
          log(exppdf(mkI,1))+...
          log(gampdf(mmuI,2,6/2))+...
          log(gampdf(mb,1.3942,1.3304));
elseif mFlag==5
    prior=log(gampdf(kE(1),1,1))+...
          log(gampdf(muE(1),2,2/2))+...
          log(gampdf(kI(1),1,1))+...
          log(gampdf(muI(1),2,6/2))+...
          log(gampdf(b(1),1.3942,1.3304));
elseif mFlag==6
    prior=log(gampdf(kE(1),1,1))+...
          log(gampdf(muE(1),2,2/2))+...
          log(gampdf(kI(1),1,1))+...
          log(gampdf(muI(1),2,6/2))+...
          sum(log(exppdf(b,mb)))+...
          log(gampdf(mb,1.3942,1.3304));
elseif mFlag==7
    prior=sum(log(exppdf(kE,mkE)))+...
          sum(log(exppdf(muE,mmuE)))+...
          log(exppdf(kI(1),1))+...
          log(gampdf(muI(1),2,6/2))+...
          log(gampdf(b(1),1.3942,1.3304))+...
          log(exppdf(mkE,1))+...
          log(gampdf(mmuE,2,2/2));
elseif mFlag==8
    prior=log(exppdf(kE(1),1))+...
          log(gampdf(muE(1),2,2/2))+...
          sum(log(exppdf(kI,mkI)))+...
          sum(log(exppdf(muI,mmuI)))+...
          log(gampdf(b(1),1.3942,1.3304))+...
          log(exppdf(mkI,1))+...
          log(gampdf(mmuI,2,6/2));
elseif mFlag==9
    prior=sum(log(exppdf(kE,mkE)))+...
          sum(log(exppdf(muE,mmuE)))+...
          sum(log(exppdf(kI,mkI)))+...
          sum(log(exppdf(muI,mmuI)))+...
          log(gampdf(b(1),1.3942,1.3304))+...
          log(exppdf(mkE,1))+...
          log(gampdf(mmuE,2,2/2))+...
          log(exppdf(mkI,1))+...
          log(gampdf(mmuI,2,6/2));
elseif mFlag==10
    prior=sum(log(exppdf(kE,mkE)))+...
          sum(log(exppdf(muE,mmuE)))+...
          sum(log(exppdf(kI,mkI)))+...
          sum(log(exppdf(muI,mmuI)))+...
          sum(log(exppdf(b,mb)))+...
          log(exppdf(mkE,1))+...
          log(gampdf(mmuE,2,2/2))+...
          log(exppdf(mkI,1))+...
          log(gampdf(mmuI,2,6/2))+...
          log(gampdf(mb,1.3942,1.3304));
end

% Constrain the infection times to lie in appropriate ranges
if mFlag<=4
    prior=prior+sum(log(unifpdf(tI(c),tVI(c,1),tVI(c,2))));
elseif mFlag>=5
    prior=prior+sum(log(unifpdf(tI(c),tIntro,tVI(c,2))));
end

% Constrain the latent periods to lie in appropriate ranges
if mFlag>=5
    prior=prior+sum(log(unifpdf(tI(c)+E(c),tVI(c,1),tVI(c,2))));
end

% Constrain the infectious periods to lie in appropriate ranges
cens=isnan(tVI(:,4));
prior=prior+sum(log(unifpdf(tI(~cens)+E(~cens)+I(~cens),...
                            tVI(~cens,3),tVI(~cens,4))));
prior=prior+sum(log(unifpdf(tI(cens)+E(cens)+I(cens),...
                            tVI(cens,3),30)));
%==========================================================================

%==========================================================================
% COMPUTE THE LOG LIKELIHOOD
% Initialise log-likelihood
logL=0;

% Set the timestep for computing the force of infection
dt=0.01;

% For each serotype ...
for s=1:nS

% Extract the infection types, latent period, infectious period, infection
% times and virus isolation times for the serotype
    iTypeS=iType(sType==s);
    ES=E(sType==s);
    IS=I(sType==s);
    tIS=tI(sType==s);
    tVIS=tVI(sType==s,:);
    cS=c(sType==s);

% Set the time points for computing the probability of infection
    t=tIntro:dt:ceil(max(tIS));

% Compute the number of animals within each pen
    n=length(iTypeS);

% Compute the number of infectious individuals at each time point
    nI=zeros(size(t));
    for a=1:length(tIS)
        x=(tIS(a)+ES(a)<t & tIS(a)+ES(a)+IS(a)>t);
        if isempty(x)==0
            nI(x)=nI(x)+1;
        end
    end

% Compute the probability of infection for each buffalo
    pI=ones(size(tIS));
    for a=1:length(pI)
        
% If the buffalo is a contact ...
        if iTypeS(a)>1

% Determine the time of infection
            x=find(tIS(a)>=t,1,'last');

% Compute the probability of infection at that time
            if isempty(x)==0
                pI(a)=(b(s).*nI(x)./n).*exp(-(b(s).*sum(dt*nI(1:x-1))/n));
            else
                pI(a)=0;
            end

        end
    end

% Add the contribution of the probability of infection to the likelihood
    logL=logL+sum(log(pI(cS)));

% Add the contribution for the latent periods to the likelihood (assuming
% there is one in the model)
    if mFlag>=5
        logL=logL+sum(log(gampdf(ES(cS),kE(s),muE(s)/kE(s))));
    end

% Add the contribution for of infectious periods for the in-contact animals
% to the likelihood, right-censoring observations as necessary
    cens=isnan(tVIS(:,4));
    logL=logL+sum(log(gampdf(IS(cS & ~cens),kI(s),muI(s)/kI(s))))...
             +sum(log(1-gamcdf(tVIS(cS & cens,3)-tIS(cS & cens)-ES(cS & cens),...
                               kI(s),muI(s)/kI(s))));

% Add the contribution of the infectious periods for the inoculated animals
% to the likelihood
    if mFlag<=4

% For the model with no latent period, I is the infectious period ...
        logL=logL+sum(log(gampdf(IS(~cS),kI(s),muI(s)/kI(s))));

    elseif mFlag>=5

% For the model with a latent period, I is the time at which the infectious
% period ends and it could have begun any time between t=0 and t=2 ...
        logL=logL+sum(log(gamcdf(IS(~cS),kI(s),muI(s)/kI(s))-...
                          gamcdf(IS(~cS)-2,kI(s),muI(s)/kI(s))));

    end
    
end
%==========================================================================
