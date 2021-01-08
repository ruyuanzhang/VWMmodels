function LLH = compute_LLH_cosSA(pars, data, gvar)
% <pars>:1x4 vector, [Jm K Jf muf]
% <data>: error_idx
% <gvar>: some auxlinary variables

Jm = pars(1); % Mean unit resource chunk
K = pars(2); % Capacity
Jf = pars(3); % Fluctuation of unit resource  
muf = pars(4); % Fluctuation of bias

N_vec = unique(data.N);
K = floor(K); % K must be a discrete quantity

%compute model predictions and log likelihood
LLH = 0;

J1=exp(Jm+Jf*cos(4*gvar.sample)); % Stimulus-specific precision fluctuation
bias=muf*sin(4*gvar.sample); % Bias fluctuation

for ii=1:length(N_vec)
    N = N_vec(ii); % memory load
    if N <= K % memory load within capacity with probability K/N
        
        S1 = floor(K/N);
        S2 = floor(K/N)+1;
        p1 = 1-mod(K,N)/N;
        p2 = mod(K,N)/N;
        
        % convert to kappa parameter
        kappa1 = interp1(gvar.J_map, gvar.kappa_map, S1*J1, 'pchip');
        kappa2 = interp1(gvar.J_map, gvar.kappa_map, S2*J1, 'pchip');
        
        % calculate p_error table
        p_error = zeros(length(gvar.error_range),length(gvar.sample));
        for jj=1:length(gvar.sample)
            mu=gvar.error_range-bias(jj);
            p_error(:,jj) = p2 * 1/(2*pi)./besseli0_fast(kappa2,1).*exp(kappa2.*cos(mu)) + ...
                    p1 * 1/(2*pi)./besseli0_fast(kappa1,1).*exp(kappa1.*cos(mu));
            p_error(:,jj) = p_error(:,jj) / sum(p_error(:,jj));           
        end

    else % memory exceeds capacity with probability 1-K/Ns
        
        kappa1 = interp1(gvar.J_map, gvar.kappa_map, J1, 'pchip');        
        p_error = zeros(length(gvar.error_range),length(gvar.sample));
        for jj=1:length(gvar.sample) % loop error
                mu=gvar.error_range-bias(jj);
                p_error(:,jj) = (K/N)*1/(2*pi) ./ besseli0_fast(kappa1, 1).*exp(kappa1.*cos(mu)) + ...
                (1-K/N)*(1/180);
                p_error(:,jj) = p_error(:,jj) / sum(p_error(:,jj))/2;
        end        
    end
    
    p_error(p_error==0) = eps;    % avoid nan in log
    % compute probabilities of reponses, take log, and sum
    error=data.error_idx{ii}; sample=data.sample_idx{ii};
    p_resp=zeros(1,length(error));
    for t=1:length(data.error_idx{ii})
        p_resp(t) = p_error(error(t),sample(t));
    end
    
    LLH = LLH - sum(log(p_resp));  % note that this LLH is a negative value
end

% This should never happen
if isnan(LLH) || isinf(abs(LLH))
    LLH = exp(666);
end
end