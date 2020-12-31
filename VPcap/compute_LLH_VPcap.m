function LLH = compute_LLH_VPcap(pars, data, gvar)
J1bar = pars(1);
power = pars(2);
tau = pars(3);
kappa_r = pars(4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K = pars(5);  % the capacity parameters
K = floor(K); % K must be an integer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


N_vec = unique(data.N);
if numel(N_vec)==1
    power = 0;   % such that Jbar = J1bar
end

%compute model predictions and log likelihood
LLH = 0;
for ii=1:length(N_vec)
    %draw values for J and compute corresponding kappas
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % judge whether N exceeds item capacity K, 
    % If yes, then can only remember K items
    % If not, can remember N item
    n_to_rememb = N_vec(ii)*(K>=N_vec(ii)) + K*(K<N_vec(ii));  % how many item to remember
 
    Jbar = ones(1,gvar.nMCSamples)*J1bar*n_to_rememb^(-power);  % note that we use negative power here.
    J = gamrnd(Jbar/tau, tau);
    J = min(J,max(gvar.J_map));
    kappa = interp1(gvar.J_map,gvar.kappa_map,J,'pchip');
    
    %compute response distribution (add motor noise, marginalize over J)
    kappa_sq = kappa.^2;
    k_c = zeros(length(gvar.error_range),numel(J));
    for jj=1:length(gvar.error_range)
        k_c(jj,:) = sqrt(kappa_r^2 + kappa_sq + 2*kappa_r*kappa*cos(gvar.error_range(jj)));
    end
    p_error = mean(bsxfun(@rdivide,besseli0_fast(k_c,1),2*pi*besseli0_fast(kappa,1)*besseli0_fast(kappa_r,1)).*exp(bsxfun(@minus, k_c, kappa+kappa_r)),2);
    
    % Normalize to 1
    % If use PROBABILITY DENSITY
    % p_error = p_error/sum(p_error) * 1/diff(gvar.error_range(1:2));
    % If use PROBABILITY
    p_error = p_error/sum(p_error);
    
    %compute probabilities of reponses, take log, and sum
    p_resp = p_error(data.error_idx{ii});
    
    if N_vec(ii) > K % N iterm exceeds capacity limits 
        p_resp = K/N_vec(ii) * p_resp + (1-K/N_vec(ii)) * 0.5/90;  
    end
    
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