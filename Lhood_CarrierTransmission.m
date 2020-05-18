function [logL, prior]=Lhood_CarrierTransmission(par,mFlag)
%
% [logL, prior]=Lhood_CarrierTransmission(par,mFlag)
%
% Matlab function to compute the log likelihood and prior for the
% transmission rate for FMDV from carrier buffalo
%
% Inputs:
% par - vector of model parameters
% mFlag - flag indicating the model to use:
%         1 - transmission rate common to all serotypes
%         2 - transmission rates differ amongst serotypes
%
% Outputs:
% logL - log likelihood for the transmission rate(s)
% prior - log prior for the transmission rate(s)

%==========================================================================
% EXTRACT THE PARAMETERS
% Transmission rates
if mFlag==1
    b=par*ones(3,1);
elseif mFlag==2
    b=par(1:3);
end

% Hierarchical parameters
if mFlag==2
    mub=par(4);
end
%==========================================================================

%==========================================================================
% COMPUTE THE PRIOR
% If the transmission rate is common to all serotypes ...
if mFlag==1
    prior=log(exppdf(b(1),100));

% If the transmission rate varies amongst serotypes ...
elseif mFlag==2
    prior=sum(log(exppdf(b,mub)))+...
          log(exppdf(mub,1));

end
%==========================================================================

%==========================================================================
% COMPUTE THE LOG LIKELIHOOD
% In the outcome for each serotype:
% t0 - time of last sampling at which all in-contact buffalo were 
%      PCR-negative (if transmission occurred)
% t1 - time of first sampling at which in-contact buffalo were PCR-positive
%      (if transmission occurred) or time of last sampling at which both
%      carriers were positive (if transmission did not occur)
% nU - number of uninfected contact buffalo
% nI - number of infected contact buffalo
% tC - time of last postive and time of first negative sample for each
%      carrier buffalo (where transmission did not occur)

% Initialise the log likelihood
logL=0;

% Add the contribution for SAT-1
% Group 1: t0=0, t1=14, nU=3, nI=7
% Group 2: t0=14, t1=44, nU=3, nI=6
logL=logL+log((exp(-b(1).*(2./12).*14.*3).*(1-(1-(1-exp(-b(1).*(2./12).*14))).^7)).*...
              (exp(-b(1).*(2./11).*44.*3).*(1-(1-(exp(-b(1).*(2./11).*14)-...
                                                  exp(-b(1).*(2./11).*44))).^6)));

% Add the contribution for SAT-2
% Group 1: no transmission, nU=10, t1=14
% Group 2: no transmission, nU=9, t1=104
logL=logL+log(exp(-b(2).*(2/12).*14.*10).*...
              exp(-b(2).*(2/11).*104.*9));

% Add the contribution for SAT-3
% Group 1: t0=0, t1=14, nU=7, nI=3
% Group 2: no transmission, nU=9, t1=44
logL=logL+log((exp(-b(3).*(2/12).*14.*7).*(1-(1-(1-exp(-b(3).*(2/12).*14))).^3)).*...
              (exp(-b(3).*(2/11).*44.*9)));
%==========================================================================
