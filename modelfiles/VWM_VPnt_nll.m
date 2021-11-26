function LLH = VWM_VPnt_nll(pars, data)
%-%-%-%-%-%-%-%-%-%-%
% Likelihood function 
% note that this function return the positive value of likelihood, so should be
% minimized
%
% <pars>: params, 1x5 vector [J1bar power tau kappa_r, pNT]
% <data>: 
%   N: set size levels
%   error_idx: error index with respect to target
%   distrError_idx: error index with respect to distractors
%   gvar: some auxilinary variables

%-%-%-%-%-%-%-%-%-%-%
J1bar = pars(1);
power = pars(2);
tau = pars(3);
kappa_r = pars(4);
pNT = pars(5); % non-target swapping error


if numel(data.unique_N)==1
    power = 0;   % such that Jbar = J1bar
end

%compute model predictions and log likelihood
LLH = 0;
for ii=1:length(data.unique_N)
    %draw values for J and compute corresponding kappas
    
    Jbar = ones(1, data.gvar.nMCSamples)*J1bar*data.unique_N(ii)^(-power);  % note we use negative power here
    
    J = gamrnd(Jbar/tau, tau); % we draw many samples of J
    J = min(J,max(data.gvar.J_map));
    kappa = interp1(data.gvar.J_map, data.gvar.kappa_map,J,'pchip'); % convert J to kappa
    
    %compute response distribution (add motor noise, marginalize over J)
    k_c = zeros(length(data.gvar.error_range),numel(J));
    for jj=1:length(data.gvar.error_range)
        k_c(jj,:) = sqrt(kappa_r^2 + kappa.^2 + 2*kappa_r*kappa*cos(data.gvar.error_range(jj)));
    end
    % calculate p_error
    p_error = bsxfun(@rdivide,besseli0_fast(k_c, 1),2*pi*besseli0_fast(kappa, 1)*besseli0_fast(kappa_r,1)).*exp(bsxfun(@minus, k_c, kappa+kappa_r));
    % we marginalize across J, to obtain the probability distribution
    p_error = mean(p_error, 2); 
    
    % normalize to 1
    % if use PROBABILITY DENSITY
    % p_error = p_error/sum(p_error * diff(gvar.error_range(1:2)));
    % if use PROBABILITY
    p_error = p_error/sum(p_error);
    p_error(p_error==0) = eps;
    
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
    
    LLH = LLH - sum(log(p_resp));  % note that this LLH is a negative value
end

% This should never happen
if isnan(LLH)
    LLH = Inf;
end
end