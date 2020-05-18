function ParEst_CarrierDuration(mFlag,nchains,nsamp,nburnin,nthin)
%
% ParEst_CarrierDuration(mFlag,nchains,nsamp,nburnin,nthin)
%
% Matlab function to implement a Bayesian MCMC scheme to estimate the
% duration of carrier state for FMDV in buffalo. This version assumes
% hierarchical structure for the parameters
%
% Inputs:
% mFlag - model to fit:
%         1-parameters common to all serotypes
%         2-means differ amongst serotypes, shape common
%         3-shapes differ amongst serotypes, mean common
%         4-both means and shape differ amongst serotypes
% nchains - number of chains to run
% nsamp - number of samples to use when estimating parameters
% nburnin - number of samples to discard before estimating parameters
% nthin - number of samples by which to thin each chain
%
% Outputs: none (N.B. All chains are saved rather than provided as output
% arguments)

%==========================================================================
% DATA ON CARRIAGE IN BUFFALO
% Duration for each serotype. In each array, the first column is the last
% positive and the second column is the first negative result for virus
% isolation. A NaN indicates no last negative result.
durC=[316 336;  % SAT-1; Maree et al. 2016, carriers only
      185 NaN;
      400 NaN;
      200 214;
      185 NaN;
      168 185;
      214 231;
	  136 155;  % SAT-2; Maree et al. 201, carriers only
      185 NaN;
      35  80;
      35  80;
      126 136;  % SAT-3; Maree et al. 2016, carriers only
      185 NaN;
      95  109;
      109 126;
      35  80];
sType=[ones(7,1); 2*ones(4,1); 3*ones(5,1)];
%==========================================================================

% Set the number of parameters
if mFlag==1
    npar=2;
elseif mFlag==2 || mFlag==3
    npar=6;
elseif mFlag==4
    npar=10;
end

% Create the arrays storing the output for the chain
ParSamp=cell(1,nchains);

% For each chain ...
parfor chain=1:nchains

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

% Sample an initial set of parameters
        if mFlag==1
            s=unifrnd(1,2);
            mu=unifrnd(150,200);
            par=[s; mu];
        elseif mFlag==2
            s=unifrnd(1,2);
            mu=unifrnd(150,200,3,1);
            pm=gamfit(mu);
            par=[s; mu; pm(1); pm(1)*pm(2)];
        elseif mFlag==3
            s=unifrnd(1,2,3,1);
            ps=gamfit(s);
            mu=unifrnd(150,200);
            par=[s; mu; ps(1); ps(1)*ps(2)];
        elseif mFlag==4
            s=unifrnd(1,2,3,1);
            ps=gamfit(s);
            mu=unifrnd(150,200,3,1);
            pm=gamfit(mu);
            par=[s; mu; ps(1); ps(1)*ps(2); pm(1); pm(1)*pm(2)];
        end

% Compute the log-likelihood
        [CurrL, prior]=Lhood_CarrierDuration(par,mFlag,durC,sType);

    end
%==========================================================================

%==========================================================================
% UPDATE THE PARAMETERS
% Sample parameter space
    disp('Sampling parameter space')
    for samp=1:nsamp+nburnin

% Indicate what's going on
        disp(['Chain: ' num2str(chain) ', Sample: ' num2str(samp) ';'...
              ' Accept: ' num2str(100*n_accept/samp,3) '%'])

% Update the variance-covariance matrix for the proposal distribution
        if samp<=nburnin && (samp<=2*npar || n_accept==0)
            SIGp=0.01*eye(npar);
        else
            SIGp=sf.*(SIG+0.01*eye(npar));
        end

% Generate the new set of probabilities
        par_new=par+mvnrnd(zeros(1,length(par)),SIGp)';

% Compute the log likelihood and prior for the new parameter set
        [NewL, prior_new]=Lhood_CarrierDuration(par_new,mFlag,durC,sType);

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
% COMPUTE DIC AND pD
% Compute the deviance for each sample
Dev=[];
PS=[];
for chain=1:nchains
    Dev=[Dev; -2*ParSamp{chain}(:,end)];
    PS=[PS; ParSamp{chain}(:,1:end-2)];
end

% Compute the mean deviance
Dbar=mean(Dev);

% Compute the deviance at the posterior mean for the parameters
Dhat=-2*Lhood_CarrierDuration(mean(PS,1)',mFlag,durC,sType);

% Compute the DIC
DIC=2*Dbar-Dhat;

% Compute the effective number of parameters
pD=Dbar-Dhat;
%==========================================================================

%==========================================================================
% SAVE THE CHAINS
% Save the outputs
save(['CarrierDuration_MCMCSamples_Model' num2str(mFlag)],...
     'ParSamp','nburnin','nsamp','DIC','pD')
%==========================================================================
