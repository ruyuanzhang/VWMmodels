function LLH = VWM_MIX_nll(pars, data)
% This is the likelihood function of mixture model
% <pars>: parameters
% <data>: data.error_idx, data.N

kappa_r = pars(1:end-1); % motor noise vector
K = pars(end); % capacity
K = floor(K); % K should be a discrete value

setSize = length(data.unique_N);
% calculate log likelihood
LLH = 0;
for ii=1:setSize
    N = data.unique_N(ii); % memory load
    kappa_r_tmp = kappa_r(ii); % for every set size, we use a kappa
    
    if N <= K % load smaller than capacity
        p_error = 1/(2*pi*besseli0_fast(kappa_r_tmp,1)) * exp(kappa_r_tmp*cos(data.gvar.error_range));
    else % load exceed capacity
        p_error = (K/N)*1/(2*pi*besseli0_fast(kappa_r_tmp, 1)*exp(kappa_r_tmp)) * exp(kappa_r_tmp*cos(data.gvar.error_range)) + (1-K/N)*(1/180);
    end
    
    % Normalize to 1, 
    % IF use PROBABILITY DENSITY
    % p_error = p_error/sum(p_error * diff(gvar.error_range(1:2)));
    % IF use PROBABILITY
    p_error = p_error/sum(p_error);
    p_error(p_error==0)= eps; % replace 0 to a small number
    
    %compute probabilities of reponses, take log, and sum
    p_resp = p_error(data.error_idx{ii});
    
    % calc loglikeli
    LLH = LLH - sum(log(p_resp));  % note that this LLH is a negative value
    
end


% This should never happen
if isnan(LLH)
    LLH = Inf;
end
end