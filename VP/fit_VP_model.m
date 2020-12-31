function [fitpars, neglh, neglhtrial, AIC, AICc, BIC] = fit_VP_model(N,probe,resp,x0, opt)

%%
data.N = N;
error = circulardiff(probe,resp,180) * 2;
error = error * pi/180;

%% discretization of the error space

error_range = linspace(-pi, pi, 181); % ori exp, error_range [-pi, pi]
% input error range [-2/pi,pi/]
gvar.error_range = error_range(1:end-1)+diff(error_range(1:2))/2;

% mapping between J and kappa
gvar.kappa_max      = 700;
gvar.kappa_map      = linspace(0,gvar.kappa_max,1e5);
gvar.J_map          = gvar.kappa_map.*besseli(1,gvar.kappa_map,1)./besseli(0,gvar.kappa_map,1);

% 
gvar.nMCSamples     = 10000;                 % number of MC samples to draw when computing model predictions (Paper: 1000)
gvar.n_par          = 4;                     % number of parameters (J1bar, power, tau, kappa_r)

% get indices of errors
unique_N = unique(data.N);
for ii=1:length(unique_N)
    trial_idx = find(data.N==unique_N(ii));
    data.error_idx{ii} = interp1(gvar.error_range,1:length(gvar.error_range),error(trial_idx), 'nearest','extrap');
end

%% ========= use bads to optimization ========
PLB = opt.PLB;
PUB = opt.PUB;
LB = opt.LB;
UB = opt.UB;
options = opt.options;
options.MaxIter = opt.options.MaxIter;
% do it
% compute_LLH should return the postive loglikelihood
[x,fval,exitflag,output,optimState,gpstruct] = bads(@(params) compute_LLH_VP(params, data, gvar),x0,LB,UB,PLB,PUB);


%% find ML parameter estimates
neglh = fval; % here neglh is a positive value, a real likelihood value should be -neglh
neglhtrial = neglh / numel(data.N);
%fitpars = X_mat(LLH==neglh,:);
fitpars = x;
% compute AIC and BIC
n_free_pars = gvar.n_par - (numel(unique_N)==1);   % Jbar, tau, power, kappa_r  (3 if we have only 1 set size)
BIC = -2*-neglh + n_free_pars*log(numel(data.N));
AICc = -2*-neglh + 2*n_free_pars+2*n_free_pars*(n_free_pars+1)/(numel(data.N)-n_free_pars-1);
AIC = -2*-neglh + 2*n_free_pars;
