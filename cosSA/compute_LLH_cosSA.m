function LLH = compute_LLH_cosSA(pars, data, gvar)

Jm = pars(1); % Mean unit resource chunk
K = pars(2); % Capacity
Jf = pars(3); % Fluctitation of unit resource  
muf = pars(4); % fluctuation of center

N_vec = unique(data.N);
K = floor(K); % K must be a discrete quantity

%compute model predictions and log likelihood
LLH = 0;

J1=exp(Jm+Jf*cos(8*gvar.sample*pi/180)); % Stimulus-specific effect
        
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
            
        % calculate kc
        p_error = zeros(length(gvar.error_range),length(gvar.sample));
        for jj=1:length(gvar.error_range)
            for ss=1:length(gvar.sample)
                mu=gvar.error_range(jj);
                p_error(jj,ss) = p2 * 1/(2*pi)/besseli0_fast(kappa2(ss),1)*exp(kappa2(ss).*cos(mu-muf*sin(4*ss))) + ...
                    p1 * 1/(2*pi)/besseli0_fast(kappa1(ss),1)*exp(kappa1(ss).*cos(mu-muf*sin(4*ss)));
            end
        end
        for i=1:length(gvar.sample)
            p_error(:,i) = p_error(:,i) / sum(p_error(:,i))/2;
        end
        % compute probabilities of reponses, take log, and sum
        error=data.error_idx{ii}; sample=data.sample_idx{ii};
        p_resp=zeros(1,length(error));
        for t=1:length(data.error_idx{ii})
            p_resp(t) = p_error(error(t),sample(t));
        end

    else % memory exceeds capacity with probability 1-K/Ns
        
        kappa1 = interp1(gvar.J_map, gvar.kappa_map, J1, 'pchip');
        
        p_error = zeros(length(gvar.error_range),length(gvar.sample));
        for jj=1:length(gvar.error_range) % loop error
            for ss=1:length(gvar.sample)
                mu=gvar.error_range(jj);
                p_error(jj,ss) = (K/N)*1/(2*pi) / besseli0_fast(kappa1(ss), 1).*exp(kappa1(ss)*cos(mu-muf*sin(4*ss))) + ...
                (1-K/N)*(0.5/90);
            end
        end
        
        for i=1:length(gvar.sample)
            p_error(:,i) = p_error(:,i) / sum(p_error(:,i))/2;
        end
        % compute probabilities of reponses, take log, and sum
        error=data.error_idx{ii}; sample=data.sample_idx{ii};
        p_resp=zeros(1,length(error));
        for t=1:length(data.error_idx{ii})
            p_resp(t) = p_error(error(t),sample(t));
        end
    end
    
    LLH = LLH + sum(log(p_resp));  % note that this LLH is a negative value
end

% We output postive likelihood, maximizing negative likelihood is equivalent to
% miniziming postive likelihood.
LLH = -LLH;

% This should never happen
if isnan(LLH) || isinf(LLH)
    LLH = 1e4;
end
end