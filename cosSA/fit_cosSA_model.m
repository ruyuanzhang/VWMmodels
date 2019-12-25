function [fitpars, neglh, neglhtrial, AIC, AICc, BIC] = fit_cosSA_model(N, probe, resp, x0, opt)
% This function is to fit cosSA VWM model
% 
% Input:
% <N>: vector of set size
% <probe>: vector of probe stimuli, range [0,180)
% <resp>: vector of subject choice, range [0,180)
% <x0>: initial parameter guess for this fit
% <opt>: options for optimization
%
% Output:
% <fitpars>: a vector of fitted parameters
% <neglh>,<AIC>,<BIC>: postive likelihood value,AIC, BIC values


%% convert data
data.N = N;
error = circulardiff(probe,resp, 180);
error = error*pi/180; % convert error to [-pi/2,pi/2];
probe = probe*pi/180;
resp = resp*pi/180;


%%
% discretization of the error space
%error_range = linspace(0,pi,91); % color exp, error_range [0, pi], input error
%range is [-pi, pi]

error_range = linspace(0, pi/2, 91); % ori exp, error_range [0, pi/2]
gvar.error_range = error_range(1:end-1)+diff(error_range(1:2))/2;
gvar.sample = linspace(pi/180, pi, 180);

% Mapping between J and kappa
gvar.kappa_max      = 700; % this number is due to matlab limit
gvar.kappa_map      = linspace(0, gvar.kappa_max, 1e5);
gvar.J_map          = gvar.kappa_map.*besseli(1,gvar.kappa_map,1)./besseli(0,gvar.kappa_map,1);

% number of parameters
gvar.n_par          = 4;                     % number of parameters (Jm, K, Jf, muf)

% get indices of errors
unique_N = unique(N);
for ii=1:length(unique_N)
    trial_idx = find(N==unique_N(ii));
    data.error_idx{ii} = interp1(gvar.error_range, 1:length(gvar.error_range), abs(error(trial_idx)),'nearest','extrap');
    data.sample_idx{ii} = interp1(gvar.sample, 1:length(gvar.sample), abs(probe(trial_idx)),'nearest','extrap');
end

%% ========= use bads to optimization ========
PLB = opt.PLB; % (Jm, K, Jf, muf)
PUB = opt.PUB;
LB = opt.LB;
UB = opt.UB;
options = opt.options;
options.MaxIter = opt.options.MaxIter;

% do it
% compute_LLH should return the postive loglikelihood
[x,fval,exitflag, output, optimState, gpstruct] = bads(@(params) compute_LLH_cosSA(params, data, gvar),x0,LB,UB,PLB,PUB);

% find ML parameter estimates
neglh = fval; % here max_lh is a positive value, a real likelihood value should be -max_lh
neglhtrial = neglh / numel(data.N);
%fitpars = X_mat(LLH==max_lh,:);
fitpars = x;
% compute AIC and BIC
n_free_pars = gvar.n_par;   % Jm, K, Jf, muf 
BIC = -2*-neglh + n_free_pars*log(numel(data.N));
AICc = -2*-neglh + 2*n_free_pars+2*n_free_pars*(n_free_pars+1)/(numel(data.N)-n_free_pars-1);
AIC = -2*-neglh + 2*n_free_pars;

