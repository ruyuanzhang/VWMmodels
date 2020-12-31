function [fitpars, neglh, neglhtrial, AIC, AICc, BIC] = fit_IL_model(N, probe, resp, x0, opt)
% This is the main function to fit item-limit model
%
% <N>: a vector containing set size levels, discrete
% <probe>: probe color, [1, 180], can be continous or discrete
% <resp>: resp color, [1, 180], can be continous or discrete


%%
data.N=N;
error=circulardiff(probe,resp,180) * 2; % [-180, 178]
error = error * pi/180; % convert to radians

% discretesize of the error space
error_range = linspace(-pi,pi, 181); % ori exp, error_range [-pi, pi)
gvar.error_range = error_range(1:end-1)+diff(error_range(1:2))/2;

% get indices of errors
errorRange_tmp = gvar.error_range; 

unique_N = unique(N);
for ii=1:length(unique_N)
    trial_idx = find(N==unique_N(ii));
    % this step is necessary if errors are continous...
    data.error_idx{ii} = interp1(errorRange_tmp, 1:length(errorRange_tmp), error(trial_idx),'nearest','extrap');
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
[x,fval,exitflag, output, optimState, gpstruct] = bads(@(params) compute_LLH_IL(params, data, gvar),x0,LB,UB,PLB,PUB);


%% find ML parameter estimates
neglh = fval; % here neglh is a positive value, a real likelihood value should be -neglh
neglhtrial = neglh/numel(data.N);
%fitpars = X_mat(LLH==neglh,:);
fitpars = x;
% compute AIC and BIC
n_free_pars = gvar.n_par;   %
BIC = -2*-neglh + n_free_pars*log(numel(data.N));
AICc = -2*-neglh + 2*n_free_pars+2*n_free_pars*(n_free_pars+1)/(numel(data.N)-n_free_pars-1);
AIC = -2*-neglh + 2*n_free_pars;

