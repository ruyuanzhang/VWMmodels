function LLH = VWM_SAnt_nll(pars, data)
% <pars>: params 1x4 vector, [J1, K, Kappa_r, pNT]
%   J1: unit resources
%   K: capacity
%   kappa_r: motor noise
%   pNT: non-target swap probability
% <data>:
%   N: a vector containing set size levels, discrete
%   error: [-90,89],response error between response and error
%   distrError: a cell, each element contains Error with respect to distractors [-90, 89]
%   gvar: some auxilinary variables

J1 = pars(1); % unit resource chunk
K = pars(2); % capacity
kappa_r = pars(3); % motor noise
pNT = pars(4); % non-target swap probability


K = floor(K); % K must be a discrete quantity

% compute model predictions and log likelihood
LLH = 0;
for ii=1:length(data.unique_N)
    N = data.unique_N(ii); % memory load
    if N <= K % memory load within capacity with probability K/N
        S1 = floor(K/N); % low
        S2 = floor(K/N)+1; % high
        p1 = 1-mod(K,N)/N;
        p2 = mod(K,N)/N;
        
        % convert to kappa parameter
        kappa1 = interp1(data.gvar.J_map, data.gvar.kappa_map, S1*J1, 'pchip');
        kappa2 = interp1(data.gvar.J_map, data.gvar.kappa_map, S2*J1, 'pchip');
        
        k_c1 = sqrt(kappa1^2+kappa_r^2+2*kappa1*kappa_r*cos(data.gvar.error_range));
        k_c2 = sqrt(kappa2^2+kappa_r^2+2*kappa2*kappa_r*cos(data.gvar.error_range));
        p_error = p2 * 1/(2*pi)*besseli0_fast(k_c2,1)/(besseli0_fast(kappa2,1)*besseli0_fast(kappa_r,1)).*exp(k_c2-kappa2-kappa_r) + ...
            p1 * 1/(2*pi)*besseli0_fast(k_c1, 1)/(besseli0_fast(kappa1, 1)*besseli0_fast(kappa_r,1)).*exp(k_c1-kappa1-kappa_r);                       

    else % Memory exceeds capacity with probability 1-K/N        
        kappa1 = interp1(data.gvar.J_map, data.gvar.kappa_map, J1, 'pchip');
        k_c1 = sqrt(kappa1^2+kappa_r^2+2*kappa1*kappa_r*cos(data.gvar.error_range));
        p_error = (K/N)*1/(2*pi) * besseli0_fast(k_c1,1) / (besseli0_fast(kappa1, 1)*besseli0_fast(kappa_r,1)).*exp(k_c1-kappa1-kappa_r) + ...
                (1-K/N)*(1/180);                
    end
    
    % compute probabilities of responses, normalize to 1
    % If use PROBABILITY DENSITY
    %p_error = p_error / sum(p_error * diff(gvar.error_range(1:2)));
    % Or If use PROBABILITY
    p_error = p_error / sum(p_error);
    p_error(p_error==0) = eps;
    
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
    
    LLH = LLH - sum(log(p_resp));  % note that this LLH is a negative value
end

% This should never happen
if isnan(LLH)
    LLH = Inf;
end
end