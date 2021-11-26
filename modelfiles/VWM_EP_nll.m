function LLH = VWM_EP_nll(pars, data)
% <pars>: 1x3 params, [J1bar power kappa_r]
% <data>: 
%   error_idx:
%   data.gvar: some auxilinary variables

J1bar = pars(1);
power = pars(2);
kappa_r = pars(3);

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
    % p_error = p_error/sum(p_error * diff(data.gvar.error_range(1:2)));
    % if use PROBABILITY
    p_error = p_error/sum(p_error);
    p_error(p_error==0) = eps; % avoid nan in log
    
    % Compute probabilities of reponses, take log, and sum
    p_resp = p_error(data.error_idx{ii});
    
    LLH = LLH - sum(log(p_resp));  % note that we output negative loglikelihood
end

% This should never happen
if isnan(LLH)
    LLH = Inf;
end
end