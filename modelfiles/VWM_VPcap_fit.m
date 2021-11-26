function c = VWM_VPcap_fit(data, c)
% This is the main function to fit EP model
% <data>:
%   N: a vector containing set size levels
%   error: a vector of response error (probe-resp)
%   distrError: a cell, response errors with respect to distractors 
% 
% <c>: a struct, contains configuration variables, e.g. optimization

%% preprocess visual working memory data
% you can do any preprocessing here. This could be model specific
data = prepVWMdata(data);
% add more helping variable to data
data.gvar.kappa_max      = 700; % this is due to matlab limit
data.gvar.kappa_map      = linspace(0, data.gvar.kappa_max,1e5);
data.gvar.J_map          = data.gvar.kappa_map.*besseli(1, data.gvar.kappa_map,1)./besseli(0, data.gvar.kappa_map,1);
data.gvar.nMCSamples     = 10000;                 % number of MC samples to draw when computing model predictions (Paper: 1000)

%% ========= use bads to optimization ========
% below is model specific
c.result.fitResults= nan(c.opt.nFit, c.opt.nvars + 5); % 5 dimension include, neglhtrial, neglh, AIC, AICc, BIC
c.result.labels = [c.opt.paramLabels {'neglhtrial', 'neglh', 'AIC', 'AICc', 'BIC'}];

% do it
objfun = @(params) c.negLogLikeliFun(params, data); % define objective function to minimize
for iFit=1:c.opt.nFit
    % optimizatiom
    [fitpars, neglh, exitflag, output, optimState, gpstruct] = bads(objfun, c.opt.x0(iFit,:),...
        c.opt.LB, c.opt.UB, c.opt.PLB, c.opt.PUB, [], c.opt.options);
    
    % calculate some useful variable
    neglhtrial = neglh/numel(data.N);
    % compute AIC, AICc and BIC
    BIC = -2*-neglh + c.opt.nvars*log(numel(data.N));
    AICc = -2*-neglh + 2*c.opt.nvars+2*c.opt.nvars*(c.opt.nvars+1)/(numel(data.N)-c.opt.nvars-1);
    AIC = -2*-neglh + 2*c.opt.nvars;
    c.result.fitResults(iFit,:) = [fitpars neglhtrial neglh AIC AICc BIC];
end

%% we choose best fit model based on AIC
[~,idx] = min(c.result.fitResults(:,end-3));
c.result.bestFit = c.result.fitResults(idx,:);
