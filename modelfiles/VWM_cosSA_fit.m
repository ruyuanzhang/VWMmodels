function c = VWM_cosSA_fit(data, c)
% This is the main function to fit item-limit model
% <data>:
% <c>: configuration struct includes optimization/

%% preprocess visual working memory data
data.error = data.error * 2 * pi/180; % convert error [-90, 89] to [-pi, pi];
data.probe = data.probe * 2 * pi/180 - pi; % convert probe from [1, 179] to [-pi, pi];

% discretization of the error space 
error_range = linspace(-pi, pi, 181); % to fit ori exp, error_range [0, pi/2]
data.gvar.error_range = error_range(1:end-1)+diff(error_range(1:2))/2;

data.gvar.sample = linspace(pi/180, 2*pi, 180);
data.gvar.kappa_max      = 700; % this is due to matlab limit
data.gvar.kappa_map      = linspace(0,data.gvar.kappa_max,1e5);
data.gvar.J_map          = data.gvar.kappa_map.*besseli(1,data.gvar.kappa_map,1)./besseli(0,data.gvar.kappa_map,1);

% get indices of errors
data.unique_N = unique(data.N);
for ii=1:length(data.unique_N)
    trial_idx = find(data.N==data.unique_N(ii));
    data.error_idx{ii} = interp1(data.gvar.error_range, 1:length(data.gvar.error_range), data.error(trial_idx), 'nearest','extrap');
    data.sample_idx{ii} = interp1(data.gvar.sample, 1:length(data.gvar.sample), data.probe(trial_idx), 'nearest','extrap');
end

%% ========= use bads to optimization ========
% below is model specific
c.result.fitResults = nan(c.opt.nFit, c.opt.nvars + 5); % 5 dimension include, neglhtrial, neglh, AIC, AICc, BIC
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
    c.result.fitResults(iFit,:) = [fitpars neglhtrial neglh AIC AICc BIC];
end

%%
% we choose best fit model based on AIC
[~,idx] = min(c.result.fitResults(:,end-3));
c.result.bestFit = c.result.fitResults(idx,:);
