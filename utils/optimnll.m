function c = optimnll(data, c)
% run the optimization of negative loglikelihood multiple times and save results and several models
% for model comparison
%
% <c>: a structure contains information for optimizations
%   opt: optimization structure, with fields
%       PLB, PUB, LB, UB: low and upper parameters boundary
%       paramsLabels: a cell, string labels for parameters
%       options: optimization options
%       nFit: int, number of optimizations to run
%       nvars: int, number of free parameters
%       x0: nFit x nvars matrix, with each row as initial values for each
%           fit       
%   negLogLiliFun: the negative log likelihood function that accept params and
%       data as input
% 
% <data>: the data structure, see details for c.negLogLikeLiFun
%
% we fit the model and save results to c.result structure with following fields
%   fitResults: a nFit x nvars matrix, with each row as the fitted parameters
%   modelMetrics: a nFit x 5 matrix, we calculate 5 model metrics,
%       currently we have neglhtrial (negative log likelihood per trial),
%       neglh, AIC, AICc, BIC.
%   modelMetricLabels: a cell with model metric names
%   bestFit: a vector, best fitted parameters with lowest neglh
 
c.result.fitResults = nan(c.opt.nFit, c.opt.nvars); % 
c.result.modelMetrics = nan(c.opt.nFit, 5); % 5 dimension include, neglhtrial, neglh, AIC, AICc, BIC
c.result.modelMetricLabels = {'neglhtrial', 'neglh', 'AIC', 'AICc', 'BIC'};

% do it
objfun = @(params) c.negLogLikeliFun(params, data); % define objective function to minimize
for iFit=1:c.opt.nFit
    % optimizatiom
    [fitpars, neglh, exitflag, output, optimState, gpstruct] = bads(objfun, c.opt.x0(iFit,:),...
        c.opt.LB, c.opt.UB, c.opt.PLB, c.opt.PUB, [], c.opt.options);
    
    % calculate model comparison variables 
    neglhtrial = neglh/numel(data.N);
    % compute AIC, AICc and BIC
    BIC = -2*-neglh + c.opt.nvars*log(numel(data.N));
    AICc = -2*-neglh + 2*c.opt.nvars+2*c.opt.nvars*(c.opt.nvars+1)/(numel(data.N)-c.opt.nvars-1);
    AIC = -2*-neglh + 2*c.opt.nvars;
    c.result.fitResults(iFit,:) = fitpars;
    c.result.modelMetrics(iFit,:) = [neglhtrial neglh AIC AICc BIC];
end

%%
% we choose best fit parameters based on negative loglilelihood
[~,idx] = min(c.result.modelMetrics(:,end-3));
c.result.bestFit = c.result.fitResults(idx,:);