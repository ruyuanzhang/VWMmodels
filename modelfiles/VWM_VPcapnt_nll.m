function LLH = VWM_VPcapnt_nll(pars, data)
J1bar = pars(1);
power = pars(2);
tau = pars(3);
kappa_r = pars(4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K = pars(5);  % the capacity parameters
K = floor(K); % K must be an integer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pNT = pars(6);

if numel(data.unique_N)==1
    power = 0;   % such that Jbar = J1bar
end

%compute model predictions and log likelihood
LLH = 0;
for ii=1:length(data.unique_N)
    %draw values for J and compute corresponding kappas
    N = data.unique_N(ii);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % judge whether N exceeds item capacity K, 
    % If yes, then can only remember K items
    % If not, can remember N item
    n_to_rememb = data.unique_N(ii)*(K>=data.unique_N(ii)) + K*(K<data.unique_N(ii));  % how many item to remember
 
    Jbar = ones(1,data.gvar.nMCSamples)*J1bar*n_to_rememb^(-power);  % note that we use negative power here.
    J = gamrnd(Jbar/tau, tau);
    J = min(J,max(data.gvar.J_map));
    kappa = interp1(data.gvar.J_map,data.gvar.kappa_map,J,'pchip');
    
    %compute response distribution (add motor noise, marginalize over J)
    kappa_sq = kappa.^2;
    k_c = zeros(length(data.gvar.error_range),numel(J));
    for jj=1:length(data.gvar.error_range)
        k_c(jj,:) = sqrt(kappa_r^2 + kappa_sq + 2*kappa_r*kappa*cos(data.gvar.error_range(jj)));
    end
    p_error = mean(bsxfun(@rdivide,besseli0_fast(k_c,1),2*pi*besseli0_fast(kappa,1)*besseli0_fast(kappa_r,1)).*exp(bsxfun(@minus, k_c, kappa+kappa_r)),2);
    
    % Normalize to 1
    % If use PROBABILITY DENSITY
    % p_error = p_error/sum(p_error) * 1/diff(gvar.error_range(1:2));
    % If use PROBABILITY
    p_error = p_error/sum(p_error);
    p_error(p_error==0) = eps;
    
    %compute probabilities of reponses, take log, and sum
    p_resp = p_error(data.error_idx{ii});
    
    %
    if data.unique_N(ii) > K % N iterm exceeds capacity limits 
        p_resp = K/data.unique_N(ii) * p_resp + (1-K/data.unique_N(ii)) * 1/180;  
    end
    
    % consider non-target swap error
    p_resp_tmp1 = p_resp;
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