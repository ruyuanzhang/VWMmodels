function LLH = VWM_MIXnt_nll(pars, data)
% This is the likelihood function of mixture model
% <pars>: parameters
% <data>: data.error_idx, data.N
%   N: a vector of set size levels
%   error_idx: a cell, each element is indices of error [1, 180].
%   distrError_idx: a length(N) x nTrial cell, each element is indices of
%       error[1, 180].
%   gvar: some auxilinary variables, to improve computational efficiency

kappa_r = pars(1:end-1); % motor noise vector
K = pars(end-1); % capacity
pNT=pars(end); % non-target swap probability

setSize = length(data.unique_N);
K = floor(K); % K should be a discrete value

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
    
    % get the probability of true response distribution    
    % no swap probability
    p_resp_tmp1 = p_error(data.error_idx{ii});
    % swap probability
    if N>1
        p_resp = nan(1, length(p_resp_tmp1));
        tmp = data.distrError_idx{ii};
        for iTrial=1:length(p_resp_tmp1)
            p_resp_tmp2 = p_error(tmp{iTrial});
            p_resp(iTrial) = (1-pNT)*p_resp_tmp1(iTrial) + sum(pNT/(N-1) * p_resp_tmp2);
        end
    else
        p_resp=p_resp_tmp1;
    end   
    % calc loglikeli
    LLH = LLH - sum(log(p_resp));  % note that this LLH is a negative value
    
end


% This should never happen
if isnan(LLH)
    LLH = Inf;
end
end