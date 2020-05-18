function [logL, prior]=Lhood_CarrierDuration(par,mFlag,durC,sType)
%
% [logL, prior]=Lhood_CarrierDuration(par,mFlag,durC,sType)
%
% Matlab function to compute the log likelihood and prior for a model of
% the duration of carrier state for FMDV in buffalo
%
% Inputs:
% par - vector of model parameters
% mFlag - model to fit:
%         1-parameters common to all serotypes
%         2-means differ amongst serotypes, shape common
%         3-shapes differ amongst serotypes, mean common
%         4-means and shapes differ amongst serotypes
% durC - times of last positive titre and first negative tonsil swab,
%        respectively
% sType - vector identifying serotype for each buffalo
%
% Outputs:
% logL - log likelihood
% prior - log prior probability

% Extract the parameters
if mFlag==1
    s=par(1)*ones(3,1);
    mu=par(2)*ones(3,1);
elseif mFlag==2
    s=par(1)*ones(3,1);
    mu=par(2:4);
    kM=par(5);
    muM=par(6);
elseif mFlag==3
    s=par(1:3);
    mu=par(4)*ones(3,1);
    kS=par(5);
    muS=par(6);
elseif mFlag==4
    s=par(1:3);
    mu=par(4:6);
    kS=par(7);
    muS=par(8);
    kM=par(9);
    muM=par(10);
end

% Compute the prior
if mFlag==1
    prior=log(exppdf(s(1),100))+log(exppdf(mu(1),307));
elseif mFlag==2
    prior=log(exppdf(s(1),100))+...
          sum(log(gampdf(mu,kM,muM./kM)))+...
          log(exppdf(kM,100))+...
          log(exppdf(muM,307));
elseif mFlag==3
    prior=sum(log(gampdf(s,kS,muS./kS)))+...
          log(exppdf(mu(1),307))+...
          log(exppdf(kS,100))+...
          log(exppdf(muS,100));
elseif mFlag==4
    prior=sum(log(gampdf(s,kS,muS./kS)))+...
          sum(log(gampdf(mu,kM,muM./kM)))+...
          log(exppdf(kS,100))+...
          log(exppdf(muS,100))+...
          log(exppdf(kM,100))+...
          log(exppdf(muM,307));
end

% Compute the log likelihood
c=isnan(durC(:,2));
logL=sum(log(gamcdf(durC(~c,2),s(sType(~c)),mu(sType(~c))./s(sType(~c)))-...
             gamcdf(durC(~c,1),s(sType(~c)),mu(sType(~c))./s(sType(~c)))))+...
     sum(log(1-gamcdf(durC(c,1),s(sType(c)),mu(sType(c))./s(sType(c)))));
