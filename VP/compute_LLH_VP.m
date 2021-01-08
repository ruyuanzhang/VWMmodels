function LLH = compute_LLH_VP(pars, data, gvar)
%-%-%-%-%-%-%-%-%-%-%
% Likelihood function 
% note that this function return the positive value of likelihood, so should be
% minimized
%
% <pars>: params, 1x4 vector [J1bar power tau kappa_r]
% <data>: error_idx
% <gvar>: some auxilinary variables

%-%-%-%-%-%-%-%-%-%-%
J1bar = pars(1);
power = pars(2);
tau = pars(3);
kappa_r = pars(4);

N_vec = unique(data.N);
if numel(N_vec)==1
    power = 0;   % such that Jbar = J1bar
end

%compute model predictions and log likelihood
LLH = 0;
for ii=1:length(N_vec)
    %draw values for J and compute corresponding kappas
    
    Jbar = ones(1, gvar.nMCSamples)*J1bar*N_vec(ii)^(-power);  % note we use negative power here
    
    J = gamrnd(Jbar/tau, tau); % we draw many samples of J
    J = min(J,max(gvar.J_map));
    kappa = interp1(gvar.J_map,gvar.kappa_map,J,'pchip'); % convert J to kappa
    
    %compute response distribution (add motor noise, marginalize over J)
    k_c = zeros(length(gvar.error_range),numel(J));
    for jj=1:length(gvar.error_range)
        k_c(jj,:) = sqrt(kappa_r^2 + kappa.^2 + 2*kappa_r*kappa*cos(gvar.error_range(jj)));
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
    
    %compute probabilities of reponses, take log, and sum
    p_resp = p_error(data.error_idx{ii});
    
    LLH = LLH - sum(log(p_resp));  % note that this LLH is a negative value
end

% This should never happen
if isnan(LLH)
    LLH = Inf;
end
end