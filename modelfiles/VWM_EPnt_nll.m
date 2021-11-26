function LLH = VWM_EPnt_nll(pars, data)
% <pars>: 1x3 params, [J1bar power kappa_r, pNT]
% <data>: data.error_idx, data.N
%   "N": a vector of set size levels
%   "error_idx": a cell, each element is indices of error [1, 180].
%   "distrError_idx": a length(N) x nTrial cell, each element is indices of
%       error[1, 180].
% <gvar>: some auxilinary variables, to improve computational efficiency

J1bar = pars(1);
power = pars(2);
kappa_r = pars(3);
pNT = pars(4); % non-target swap error probability

if numel(data.unique_N)==1 % only one set size level
    power = 0;   % such that Jbar = J1bar
end

%compute model predictions and log likelihood
LLH = 0;
for ii=1:length(data.unique_N)
    
    J = J1bar * data.unique_N(ii)^(-power);  % note we use negative power here    
    J = min(J,max(data.gvar.J_map));
    kappa = interp1(data.gvar.J_map, data.gvar.kappa_map, J, 'pchip'); % convert J to kappa
    
    % calculate k_c
    k_c = sqrt(kappa_r^2 + kappa.^2 + 2*kappa_r*kappa*cos(data.gvar.error_range));
    p_error = bsxfun(@rdivide, besseli0_fast(k_c,1), 2*pi*besseli0_fast(kappa,1)*besseli0_fast(kappa_r,1)).*exp(bsxfun(@minus, k_c, kappa+kappa_r));
    
    % Normalize to 1
    % if use PROBABILITY DENSITY
    % p_error = p_error/sum(p_error * diff(gvar.error_range(1:2)));
    % if use PROBABILITY
    p_error = p_error/sum(p_error);
    p_error(p_error==0) = eps; % avoid nan in log
    
    % get the probability of true response distribution    
    % no swap probability
    p_resp_tmp1 = p_error(data.error_idx{ii});
    % swap probability
    N=data.unique_N(ii);
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
    
    LLH = LLH - sum(log(p_resp));  % note that we output negative loglikelihood
end

% This should never happen
if isnan(LLH)
    LLH = Inf;
end
end