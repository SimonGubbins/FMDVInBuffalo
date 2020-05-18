function [logL, prior]=Lhood_MaternalAbDuration(par,mFlag,age0,age1)
%
% [logL, prior]=Lhood_MaternalAbDuration(par,mFlag,age0,age1)
%
% Matlab function to compute the log likelihood and prior for a model of
% the duration of maternal antibodies for FMDV in buffalo
%
% Inputs:
% par - vector of model parameters
% mFlag - flag indicating model to use:
%         1-parameters common to serotypes
%         2-parameters differ amongst to serotypes
% age0, age1 - arrays of ages of last protective titre and first
%              non-protective titre, respectively (cols for each serotype)
%
% Outputs:
% logL - log likelihood
% prior - log prior probability

% Extract the parameters
if mFlag==1
    s=par(1);
    mu=par(2);
elseif mFlag==2
    s=par(1:3);
    mu=par(4:6);
    p_s=par(7:8);
    p_mu=par(9:10);
end

% Compute the prior
if mFlag==1
    prior=log(gampdf(s,1,100))+log(gampdf(mu,1,0.5));
elseif mFlag==2
    prior=sum(log(gampdf(s,p_s(1),p_s(2)./p_s(1))))+...
          sum(log(gampdf(mu,p_mu(1),p_mu(2)./p_mu(1))))+...
          sum(log(gampdf(p_s,1,100)))+...
          sum(log(gampdf(p_mu,1,0.5)));
end

% Compute the log likelihood
if mFlag==1
    logL=sum(log(gamcdf(age1(~isnan(age1)),s,mu./s)-...
                 gamcdf(age0(~isnan(age0)),s,mu./s)));
elseif mFlag==2
    s=repmat(s',size(age0,1),1);
    mu=repmat(mu',size(age0,1),1);
    x=~isnan(age0);
    logL=sum(sum(log(gamcdf(age1(x),s(x),mu(x)./s(x))-...
                     gamcdf(age0(x),s(x),mu(x)./s(x)))));
end
