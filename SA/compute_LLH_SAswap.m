function LLH = compute_LLH_SAswap(pars, data, gvar)
%
%
%
J1 = pars(1); % unit resource chunk
K = pars(2); % capacity
kappa_r = pars(3); % motor noise
s=pars(4); % swap rate

N_vec = unique(data.N);
K = floor(K); % K must be a discrete quantity

%compute model predictions and log likelihood
LLH = 0;
for ii=1:length(N_vec)
    N = N_vec(ii); % memory load
    if N<=K % memory load within capacity with probability K/N
        S1 = floor(K/N);
        S2 = floor(K/N)+1;
        p1 = 1-mod(K,N)/N;
        p2 = mod(K,N)/N;
        
        % convert to kappa parameter
        kappa1 = interp1(gvar.J_map, gvar.kappa_map, S1*J1, 'pchip');
        kappa2 = interp1(gvar.J_map, gvar.kappa_map, S2*J1, 'pchip');
        
        % calculate kc
        p_error = zeros(length(gvar.error_range),1);
        for jj=1:length(gvar.error_range)
            k_c1 = sqrt(kappa1^2+kappa_r^2+2*kappa1*kappa_r*cos(gvar.error_range(jj)));
            k_c2 = sqrt(kappa2^2+kappa_r^2+2*kappa2*kappa_r*cos(gvar.error_range(jj)));
            p_error(jj) = p2 * 1/(2*pi)*besseli0_fast(k_c2,1)/(besseli0_fast(kappa2,1)*besseli0_fast(kappa_r,1)).*exp(k_c2-kappa2-kappa_r) + ...
                p1 * 1/(2*pi)*besseli0_fast(k_c1,1)/(besseli0_fast(kappa1,1)*besseli0_fast(kappa_r,1)).*exp(k_c1-kappa1-kappa_r);
        end
        
        p_error = p_error / sum(p_error)/2;
        %compute probabilities of reponses
        p_T = (1-s)*p_error(data.error_idx{ii});
        
        error_nt=data.error_nt_idx{ii};
        p_NT=zeros(1,length(error_nt));
        for t=1:length(data.error_nt_idx{ii})
            if N_vec(ii)>1
                for j=1:N_vec(ii)-1
                    p_NT(t)=p_NT(t)+s/(N_vec(ii)-1)*p_error(error_nt(t,j));
                end
            end
        end
        p_resp=p_T+p_NT;
        
    else % memory exceeds capacity with probability 1-K/N
        
        kappa1 = interp1(gvar.J_map, gvar.kappa_map, J1, 'pchip');
        
        p_error = zeros(length(gvar.error_range),1);
        for jj=1:length(gvar.error_range)
            k_c1 = sqrt(kappa1^2+kappa_r^2+2*kappa1*kappa_r*cos(gvar.error_range(jj)));
            p_error(jj) = (K/N)*1/(2*pi) * besseli0_fast(k_c1,1) / (besseli0_fast(kappa1, 1)*besseli0_fast(kappa_r,1)).*exp(k_c1-kappa1-kappa_r) + ...
                (1-K/N)*(0.5/90);
        end
        
        %compute probabilities of reponses
        p_error = p_error / sum(p_error)/2;
        %compute probabilities of reponses
        p_T = (1-s)*p_error(data.error_idx{ii});
        
        error_nt=data.error_nt_idx{ii};
        p_NT=zeros(1,length(error_nt));
        for t=1:length(data.error_nt_idx{ii})
            if N_vec(ii)>1
                for j=1:N_vec(ii)-1
                    p_NT(t)=p_NT(t)+s/(N_vec(ii)-1)*p_error(error_nt(t,j));
                end
            end
        end
        p_resp=p_T+p_NT;
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