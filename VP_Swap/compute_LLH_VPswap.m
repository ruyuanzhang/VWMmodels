%-%-%-%-%-%-%-%-%-%-%
% Likelihood function
% note that this function return the positive value of likelihood, so should be
% minimized
%-%-%-%-%-%-%-%-%-%-%

function LLH = compute_LLH_VPswap(pars, data, gvar)
J1bar = pars(1);
power = pars(2);
tau = pars(3);
kappa_r = pars(4);
s=pars(5);

N_vec = unique(data.N);
if numel(N_vec)==1
    power = 0;   % such that Jbar = J1bar
end

%compute model predictions and log likelihood
LLH = 0;
for ii=1:length(N_vec)
    %draw values for J and compute corresponding kappas
    
    Jbar = ones(1,gvar.nMCSamples)*J1bar*N_vec(ii)^(-power);  % note we use negative power here
    
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
    
    
    
    %make sure p_error integrates to 0.5 (we're considering only the positive half of the pdf)
    %p_error = p_error/sum(p_error) * 1/diff(gvar.error_range(1:2))/2;
    p_error = p_error/sum(p_error) / 2;
    
    %compute probabilities of reponses, take log, and sum
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