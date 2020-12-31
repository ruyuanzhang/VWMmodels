function LLH = compute_LLH_SA(pars, data, gvar)
% <pars>: params 1x2 vector, [J1, K]
%   J1: unit resources
%   K: capacity
% <data>:
%   error_idx: the index of error in [-90 89] range
% <gvar>: some auxilinary variables

J1 = pars(1); % unit resource chunk
K = pars(2); % capacity
kappa_r = pars(3); % motor noise

N_vec = unique(data.N);
K = floor(K); % K must be a discrete quantity

% compute model predictions and log likelihood
LLH = 0;
for ii=1:length(N_vec)
    N = N_vec(ii); % memory load
    if N <= K % memory load within capacity with probability K/N
        S1 = floor(K/N); % low
        S2 = floor(K/N)+1; % high
        p1 = 1-mod(K,N)/N;
        p2 = mod(K,N)/N;
        
        % convert to kappa parameter
        kappa1 = interp1(gvar.J_map, gvar.kappa_map, S1*J1, 'pchip');
        kappa2 = interp1(gvar.J_map, gvar.kappa_map, S2*J1, 'pchip');
        
        k_c1 = sqrt(kappa1^2+kappa_r^2+2*kappa1*kappa_r*cos(gvar.error_range));
        k_c2 = sqrt(kappa2^2+kappa_r^2+2*kappa2*kappa_r*cos(gvar.error_range));
        p_error = p2 * 1/(2*pi)*besseli0_fast(k_c2,1)/(besseli0_fast(kappa2,1)*besseli0_fast(kappa_r,1)).*exp(k_c2-kappa2-kappa_r) + ...
            p1 * 1/(2*pi)*besseli0_fast(k_c1, 1)/(besseli0_fast(kappa1, 1)*besseli0_fast(kappa_r,1)).*exp(k_c1-kappa1-kappa_r);                       

    else % Memory exceeds capacity with probability 1-K/N        
        kappa1 = interp1(gvar.J_map, gvar.kappa_map, J1, 'pchip');
        k_c1 = sqrt(kappa1^2+kappa_r^2+2*kappa1*kappa_r*cos(gvar.error_range));
        p_error = (K/N)*1/(2*pi) * besseli0_fast(k_c1,1) / (besseli0_fast(kappa1, 1)*besseli0_fast(kappa_r,1)).*exp(k_c1-kappa1-kappa_r) + ...
                (1-K/N)*(1/180);                
    end
    
    % compute probabilities of responses, normalize to 1
    % If use PROBABILITY DENSITY
    %p_error = p_error / sum(p_error * diff(gvar.error_range(1:2)));
    % Or If use PROBABILITY
    p_error = p_error / sum(p_error);
    p_error(p_error==0) = eps;
    % compute probabilities of reponses
    p_resp = p_error(data.error_idx{ii});
    
    LLH = LLH - sum(log(p_resp));  % note that this LLH is a negative value
end

% This should never happen
if isnan(LLH)
    LLH = Inf;
end
end