function ParEst_AcuteTransmission(mFlag,seeds,nchains,nsamp,nburnin,nthin)
%
% ParEst_AcuteTransmission(mFlag,seeds,nchains,nsamp,nburnin,nthin)
%
% Matlab function to find implement a Bayesian MCMC scheme to estimate
% transmission parameters for FMDV in buffalo based on the outcome of
% group challenge experiments.
%
% This version includes a range of models (see below) to allow for
% differences amongst serotypes. It also imputes the infectious periods.
%
% Inputs:
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
% seeds - vector of seeds (must be of length nchains) to use for the random
%         number generator for each chain
% nchains - number of chains to run
% nsamp - number of samples to use when estimating parameters
% nburnin - number of samples to discard before estimating parameters
% nthin - number of samples by which to thin each chain
%
% Outputs: none (N.B. All chains are saved rather than provided as output
% arguments)

%==========================================================================
% PREPARE THE DATA
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

% Set the time at which contact buffalo are introduced to the inoculated
% ones
tIntro=2;

% Compute the number of animals
nA=length(iType);

% Compute the number of contact animals in the study
nC=length(find(iType==2));

% Compute the number of serotypes
nS=max(sType);
%==========================================================================

% Set the number of parameters (i.e. model parameters and unobserved times
% of infection latent periods for each animal)
nmp=[3 nS+2+1 2*nS+1+2 3*nS+3 5 nS+4+1 2*nS+3+2 2*nS+3+2 4*nS+1+4 5*nS+5];
if mFlag<=4
    npar=nA+nC+nmp(mFlag);
elseif mFlag>=5
    npar=nA+2*nC+nmp(mFlag);
end

% Create the arrays storing the output for the chain
ParSamp=cell(1,nchains);

% For each chain ...
parfor chain=1:nchains
    
% Initialise the random number generator
    rng(seeds(chain),'twister');

%==========================================================================
% INITIALISE THE CHAIN
% Set the initial scaling factor for the proposal distribution
    sf=(2.38.^2)/npar;
    SIG=eye(npar);
    
% Set the counter for the number of accepted samples
    n_accept=0;

% Create the arrays storing the output for the chain
    ParSampC=zeros(nsamp/nthin,npar+2);
    iter=1;

% Generate the initial parameters for the chain, ensuring they generate a
% finite log likelihood and prior
    disp('Initialising chain')
    CurrL=NaN;
    prior=NaN;
    while ~isfinite(CurrL+prior)

% Set the initial infection times and latent and infectious periods
        if mFlag<=4
            tI=unifrnd(tVI(iType==2,1),tVI(iType==2,2));
            tIToo=zeros(size(iType));
            tIToo(iType==2)=tI;
            I=unifrnd(tVI(:,3),tVI(:,4))-tIToo;
        elseif mFlag>=5
            tI=unifrnd(0,tVI(iType==2,2));
            E=unifrnd(tI,tVI(iType==2,2))-tI;
            tIToo=zeros(size(iType));
            tIToo(iType==2)=tI;
            EToo=zeros(size(iType));
            EToo(iType==2)=E;
            I=unifrnd(tVI(:,3),tVI(:,4))-tIToo-EToo;
        end
        I(isnan(I))=unifrnd(3,9,length(find(isnan(I))));

% Sample an initial set of parameters
        if mFlag==1
            kI=unifrnd(1,10);
            muI=unifrnd(1,6);
            b=unifrnd(0,5);
            par0=[kI; muI; b];
        elseif mFlag==2
            kI=unifrnd(1,10);
            muI=unifrnd(1,6);
            b=unifrnd(0,5,nS,1);
            mb=expfit(b);
            par0=[kI; muI; b; mb];
        elseif mFlag==3
            kI=unifrnd(1,10,nS,1);
            muI=unifrnd(1,6,nS,1);
            b=unifrnd(0,5);
            mkI=expfit(kI);
            mmuI=expfit(muI);
            par0=[kI; muI; b; mkI; mmuI];
        elseif mFlag==4
            kI=unifrnd(1,10,nS,1);
            muI=unifrnd(1,6,nS,1);
            b=unifrnd(0,5,nS,1);
            mkI=expfit(kI);
            mmuI=expfit(muI);
            mb=expfit(b);
            par0=[kI; muI; b; mkI; mmuI; mb];
        elseif mFlag==5
            kE=unifrnd(1,10);
            muE=unifrnd(1,6);
            kI=unifrnd(1,10);
            muI=unifrnd(1,6);
            b=unifrnd(0,5);
            par0=[kE; muE; kI; muI; b];
        elseif mFlag==6
            kE=unifrnd(1,10);
            muE=unifrnd(1,6);
            kI=unifrnd(1,10);
            muI=unifrnd(1,6);
            b=unifrnd(0,5,nS,1);
            mb=expfit(b);
            par0=[kE; muE; kI; muI; b; mb];
        elseif mFlag==7
            kE=unifrnd(1,10,nS,1);
            muE=unifrnd(1,6,nS,1);
            kI=unifrnd(1,10);
            muI=unifrnd(1,6);
            b=unifrnd(0,5);
            mkE=expfit(kE);
            mmuE=expfit(muE);
            par0=[kE; muE; kI; muI; b; mkE; mmuE];
        elseif mFlag==8
            kE=unifrnd(1,10);
            muE=unifrnd(1,6);
            kI=unifrnd(1,10,nS,1);
            muI=unifrnd(1,6,nS,1);
            b=unifrnd(0,5);
            mkI=expfit(kI);
            mmuI=expfit(muI);
            par0=[kE; muE; kI; muI; b; mkI; mmuI];
        elseif mFlag==9
            kE=unifrnd(1,10,nS,1);
            muE=unifrnd(1,6,nS,1);
            kI=unifrnd(1,10,nS,1);
            muI=unifrnd(1,6,nS,1);
            b=unifrnd(0,5);
            mkE=expfit(kE);
            mmuE=expfit(muE);
            mkI=expfit(kI);
            mmuI=expfit(muI);
            par0=[kE; muE; kI; muI; b; mkE; mmuE; mkI; mmuI];
        elseif mFlag==10
            kE=unifrnd(1,10,nS,1);
            muE=unifrnd(1,6,nS,1);
            kI=unifrnd(1,10,nS,1);
            muI=unifrnd(1,6,nS,1);
            b=unifrnd(0,5,nS,1);
            mkE=expfit(kE);
            mmuE=expfit(muE);
            mkI=expfit(kI);
            mmuI=expfit(muI);
            mb=expfit(b);
            par0=[kE; muE; kI; muI; b; mkE; mmuE; mkI; mmuI; mb];
        end

% Combine the times of infection, latent periods and model parameters
        if mFlag<=4
            par=[tI; I; par0];
        elseif mFlag>=5
            par=[tI; E; I; par0];
        end

% Compute the log-likelihood
        [CurrL, prior]=Lhood_AcuteTransmission(par,mFlag,tVI,iType,...
                                               sType,tIntro);

    end
%==========================================================================

%==========================================================================
% ESTIMATE THE PARAMETERS
% Sample parameter space
    disp('Sampling parameter space')
    for samp=1:nsamp+nburnin

% Indicate what's going on
        disp(['Chain: ' num2str(chain) ', Sample: ' num2str(samp) ';'...
              ' Acceptance: ' num2str(100*n_accept/samp,3) '%'])

% Update the variance-covariance matrix for the proposal distribution
        if samp<=nburnin && (samp<=2*npar || n_accept==0)
            SIGp=0.001*eye(npar);
        else
            SIGp=sf.*(SIG+0.001*eye(npar));
        end

% Generate the new set of probabilities
        par_new=par+mvnrnd(zeros(1,length(par)),SIGp)';

% Compute the log likelihood and prior for the new parameter set
        [NewL, prior_new]=Lhood_AcuteTransmission(par_new,mFlag,tVI,...
                                                  iType,sType,tIntro);

% Test whether to accept the new parameter set
        u=unifrnd(0,1);
        if isfinite(NewL+prior_new) && ...
           u<min(1,exp((NewL+prior_new)-(CurrL+prior)))

% Update the counter
            n_accept=n_accept+1;

% Update the covariance matrix for the proposal distribution
            if n_accept==1
                pbar=mean([par par_new],2);
                SIG=cov([par'; par_new']);
            elseif samp<=nburnin && n_accept>1
                pbar_prev=pbar;
                pbar=(n_accept./(n_accept+1)).*pbar_prev+...
                     (1./(n_accept+1)).*par_new;
                SIG=((n_accept-1)./n_accept).*SIG+...
                    (1./n_accept).*(n_accept.*(pbar_prev*pbar_prev')-...
                                    (n_accept+1).*(pbar*pbar')+...
                                    (par_new*par_new'));
            end

% Update the chain
            CurrL=NewL;
            prior=prior_new;
            par=par_new;

        end

% Every one hundred samples during burn-in, tune the scaling factor
% for the proposal distribution to ensure an acceptance rate of 20-40%
        if samp<=nburnin && mod(samp+1,100)==1 && n_accept/samp<0.2
            sf=sf/2;
        elseif samp<=nburnin && mod(samp+1,100)==1 && n_accept/samp>0.4
            sf=2*sf;
        end
%==========================================================================

%==========================================================================
% STORE THE OUTPUT
% After burn in, save iterations of the chain, thinning as specified
        if nthin==1
            ParSampC(samp,:)=[par' prior CurrL];
        elseif samp>nburnin && mod(samp,nthin)==1
            ParSampC(iter,:)=[par' prior CurrL];
            iter=iter+1;
        end
%==========================================================================

    end

% Store the chain
    ParSamp{chain}=ParSampC;

end

%==========================================================================
% SAVE THE OUTPUTS
% Save to a file
save(['Model' num2str(mFlag) '_MCMCSamples'],...
     'ParSamp','nburnin','nsamp','seeds')
%==========================================================================
