function LLH = compute_LLH_cosSAswap(pars, data, gvar)

Jm = pars(1); % Mean unit resource chunk
K = pars(2); % Capacity
Jf = pars(3); % Fluctuation of unit resource
muf = pars(4); % Fluctuation of bias
s = pars(5); % Swap rate

N_vec = unique(data.N);
K = floor(K); % K must be a discrete quantity

%compute model predictions and log likelihood
LLH = 0;

J1=exp(Jm+Jf*cos(4*gvar.sample)); % Stimulus-specific precision fluctuation
bias=abs(muf*cos(4*gvar.sample-pi/2)); % Bias fluctuation (Given that it is absolute error in our case, we take the absolute value of bias)

for ii=1:length(N_vec)
    N = N_vec(ii); % memory load
    if N <= K % memory load within capacity with probability K/N
        
        S1 = floor(K/N);
        S2 = floor(K/N)+1;
        p1 = 1-mod(K,N)/N;
        p2 = mod(K,N)/N;
        
        % convert to kappa parameter
        kappa1 = interp1(gvar.J_map, gvar.kappa_map, S1*J1, 'pchip');
        kappa2 = interp1(gvar.J_map, gvar.kappa_map, S2*J1, 'pchip');
        
        % calculate kc
        p_error = zeros(length(gvar.error_range),length(gvar.sample));
        for jj=1:length(gvar.error_range)
            mu=gvar.error_range(jj)-bias;
            p_error(jj,:) = p2 * 1/(2*pi)/besseli0_fast(kappa2,1).*exp(kappa2.*cos(mu)) + ...
                p1 * 1/(2*pi)/besseli0_fast(kappa1,1).*exp(kappa1.*cos(mu));
        end
        for i=1:length(gvar.sample)
            p_error(:,i) = p_error(:,i) / sum(p_error(:,i))/2;
        end
        % compute probabilities of reponses, take log, and sum
        p_T=zeros(1,length(error_nt));
        error=data.error_idx{ii}; sample=data.sample_idx{ii};
        for t=1:length(data.error_idx{ii})
            p_T(t) = (1-s)*p_error(error(t),sample(t));
        end
        % Swap
        error_nt=data.error_nt_idx{ii}; sample_nt=data.sample_nt_idx{ii};
        p_NT=zeros(1,length(error_nt));
        for t=1:length(data.error_nt_idx{ii})
            if N_vec(ii)>1
                for j=1:N_vec(ii)-1
                    p_NT(t)=p_NT(t)+s/(N_vec(ii)-1)*p_error(error_nt(t,j), sample_nt(t,j));
                end
            end
        end
        p_resp=p_T+p_NT;
        
    else % memory exceeds capacity with probability 1-K/Ns
        
        kappa1 = interp1(gvar.J_map, gvar.kappa_map, J1, 'pchip');
        
        p_error = zeros(length(gvar.error_range),length(gvar.sample));
        for jj=1:length(gvar.error_range) % loop error
            mu=gvar.error_range(jj)-bias;
            p_error(jj,:) = (K/N)*1/(2*pi) / besseli0_fast(kappa1, 1).*exp(kappa1*cos(mu)) + ...
                (1-K/N)*(0.5/90);
        end
        
        for i=1:length(gvar.sample)
            p_error(:,i) = p_error(:,i) / sum(p_error(:,i))/2;
        end
        % compute probabilities of reponses, take log, and sum
        p_T=zeros(1,length(error_nt));
        error=data.error_idx{ii}; sample=data.sample_idx{ii};
        for t=1:length(data.error_idx{ii})
            p_T(t) = (1-s)*p_error(error(t),sample(t));
        end
        % Swap
        error_nt=data.error_nt_idx{ii}; sample_nt=data.sample_nt_idx{ii};
        p_NT=zeros(1,length(error_nt));
        for t=1:length(data.error_nt_idx{ii})
            if N_vec(ii)>1
                for j=1:N_vec(ii)-1
                    p_NT(t)=p_NT(t)+s/(N_vec(ii)-1)*p_error(error_nt(t,j), sample_nt(t,j));
                end
            end
        end
        p_resp=p_T+p_NT;
    end
    
    LLH = LLH + sum(log(p_resp));  % note that this LLH is a negative value
end

% We output postive likelihood, maximizing negative likelihood is equivalent to
% miniziming postive likelihood.
LLH = -LLH;

% This should never happen
if isnan(LLH) || isinf(LLH)
    LLH = exp(666);
end
end