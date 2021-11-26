function LLH = VWM_ILnt_nll(pars, data)
% This is the likelihood function of item-limited model  
% <pars>: parameters
% <data>: data.error_idx, data.N
%   N: a vector of set size levels
%   error_idx: a cell, each element is indices of error [1, 180].
%   distrError_idx: a length(N) x nTrial cell, each element is indices of
%       error[1, 180].
%   gvar: some auxilinary variables, to improve computational efficiency

kappa_r = pars(1); % motor noise
K = pars(2); % capacity
pNT = pars(3); % non-target swap probablity

data.unique_N = unique(data.N);
K = floor(K); % K should be a discrete value;

% Computational modeling predictions and log likelihood
LLH = 0;
for ii=1:length(data.unique_N) % loop setSize level
    N = data.unique_N(ii); % memory load
    
    % We discretize errors [-pi/2, pi/2) into 180 equal space range, then
    % calculate their probability
    if N <= K % load smaller than capacity
        p_error = 1/(2*pi*besseli0_fast(kappa_r, 1)) * exp(kappa_r*cos(data.gvar.error_range));
    else % load exceed capacity
        p_error = (K/N)*1/(2*pi*besseli0_fast(kappa_r,1)*exp(kappa_r)) * exp(kappa_r*cos(data.gvar.error_range)) + (1-K/N)*(1/180);
    end
    
    % if use PROBABILITY DENSITY
    %p_error = p_error / sum(p_error * diff(gvar.error_range(1:2))); % pi/180 is the binWidth
    % if use PROBABILITY
    p_error = p_error / sum(p_error); % pi/180 is the binWidth
    p_error(p_error==0) = eps; % to avoid 0
    
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
    
    % calc negative loglikeli
    LLH = LLH - sum(log(p_resp));  %
end

% This should never happen
if isnan(LLH)
    LLH = Inf;
end
end