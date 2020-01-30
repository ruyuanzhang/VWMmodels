function LLH = compute_LLH_MIXswap(pars, data, gvar)
kappa_r = pars(1:end-2); % motor noise vector
K = pars(end-1); % capacity
s = pars(end); % swap rate

N_vec = unique(data.N);
setSize = length(N_vec);
K = floor(K); % K should be a discrete value;

% compute model predictions and log likelihood
LLH = 0;
for ii=1:setSize
    N = N_vec(ii); % memory load
    kappa_r_tmp = kappa_r(ii);
    p_error = zeros(length(gvar.error_range),1);
    for jj=1:length(gvar.error_range)
        if N <= K % load smaller than capacity
            p_error(jj) = 1/(2*pi*besseli0_fast(kappa_r_tmp,1)) * exp(kappa_r_tmp*cos(gvar.error_range(jj)));
        else % load exceed capacity
            p_error(jj) = (K/N)*1/(2*pi*besseli0_fast(kappa_r_tmp,1)*exp(kappa_r_tmp)) * exp(kappa_r_tmp*cos(gvar.error_range(jj))) + (1-K/N)*(0.5/90);
        end
    end
    
    % because we only consider [0,90],full distribution should be 0.5
    p_error = p_error/sum(p_error) / 2;
    %compute probabilities of reponses, take log, and sum
    p_T = (1-s)*p_error(data.error_idx{ii});
    
    error_nt=data.error_nt_idx{ii};
    p_NT=zeros(1,length(error_nt));
    for t=1:length(data.error_nt_idx{ii})
        if N_vec(ii)>1
            for j=1:N_vec(ii)-1
                p_NT(t)=p_NT(t)+s/(N_vec(ii)-1)*p_error(error_nt(t,j));
            end
        end
    end
    p_resp=p_T+p_NT;
    
    % calc loglikeli
    LLH = LLH + sum(log(p_resp));  % note that this LLH is a negative value
    
end

% We output postive likelihood, maximizing negative likelihood is equivalent to
% miniziming postive likelihood.
LLH = -LLH;

% This should never happen
if isnan(LLH)
    LLH = Inf;
end
end