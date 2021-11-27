function c = VWM_ILnt_fit(data, c)
% This is the main function to fit item-limit model
% <data>:
% <c>: configuration struct includes optimization/

%% preprocess visual working memory data
% you can do any preprocessing here
data = prepVWMdatant(data);

%% ========= use bads to optimization ========
% below is model specific
c.result.fitResults = nan(c.opt.nFit, c.opt.nvars); % 
c.result.modelMetrics = nan(c.opt.nFit, 5); % 5 dimension include, neglhtrial, neglh, AIC, AICc, BIC
c.result.labels = [c.opt.paramLabels {'neglhtrial', 'neglh', 'AIC', 'AICc', 'BIC'}];

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

