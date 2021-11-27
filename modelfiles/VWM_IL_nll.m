function LLH = VWM_IL_nll(pars, data)
% This is the negative likelihood function of item-limited model  
% <pars>: parameters
% <data>: data.error_idx, data.N
%   error_idx: index of error [1, 180]
%   N: a vector of set size levels
%   gvar: some auxilinary variables, to improve computational efficiency

kappa_r = pars(1); % motor noise
K = pars(2); % capacity

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
    p_error(p_error==0) = eps;
    
    % get the probability of true response distribution
    p_resp = p_error(data.error_idx{ii});
    
    % calc negative loglikeli
    LLH = LLH - sum(log(p_resp));  %
end

% This should never happen
if isnan(LLH)
    LLH = Inf;
end
end