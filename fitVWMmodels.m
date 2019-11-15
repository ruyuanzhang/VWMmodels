function results = fitVWMmodels(data, models, prefix)
% function results = fitVWMmodels(data, models, prefix)
% 
% <data>: a structure
%   data.probe: 1-180 deg, probe color
%   data.resp: 1-180 deg, response color
%   data.error: -90~89 deg, circular error
%   data.N: a vector of integer, set-size level
% <models>: a cell vector of model name to fit, or 'all' to fit all
%   avaliable models.
% <prefix>: prefix of the result file.
%
% Here we tried to fit a set of VWM models, including
%
% Limited-item(IL) model:
% Mixture (MIX) model:
% cosine Slot-plus-averaging cosine (cosSA) model:
% Slot-plus-averaging(SA) model:
% Equal precision(EP):
% Variable-precision (VP) model:
% Variable-precision-with-capacity (Vpcap) model:
%
%
% History
%   20191027 RZ wrote 1st version
%   20191029 RZ let the model output likelihood per trial to compensate 
%       different trials across subjects and groups
%
%
%


if notDefined('prefix')
    prefix='';
end

%% Settings
if ~iscell(models)
    models = {models};
end
models = upper(models);

setSize = unique(data.N);
nSetSize = length(unique(data.N));
nFit = 20;
maxJ = 700; 

results.setSize = setSize;
results.models = models;
results.nFit = nFit; % for each model, we random initialize parameters and fit <nFit> times 
results.data = data;
results.allposibModels = {'IL', 'MIX', 'cosSA', 'SA', 'EP', 'VP', 'VPcap'};
results.allposibModelParams = {
    {'Kr','cap', 'neglh', 'neglhtrial', 'AIC', 'BIC'}, ... % IL
    {'Kr for different set size levels','cap', 'neglh', 'neglhtrial','AIC', 'BIC'}, ... % MIX
    {'J1','K', 'Kr', 'neglh', 'neglhtrial','AIC', 'BIC'}, ... % SA
    {'J1bar','power', 'Kr', 'neglh', 'neglhtrial', 'AIC', 'BIC'}, ... % EP
    {'J1bar','power', 'tau', 'Kr', 'neglh', 'neglhtrial', 'AIC', 'BIC'}, ... % VP
    {'J1bar','power', 'tau', 'Kr', 'K','neglh', 'neglhtrial', 'AIC', 'BIC'}, ... % Vpcap
    {'Jm','K', 'Jf', 'muf', 'neglh', 'neglhtrial', 'AIC', 'BIC'}, ... % cosSA
    };

if numel(models)==1 && strcmp(models{1}, 'ALL')
    results.nModeltoFit = numel(results.AllModels);
else
    results.nModeltoFit = numel(models);
end

%% add BADS optimization toolbox
addpath(genpath('./'));

%% ========================================================================
% ===================== Now fit models ====================================
%% ========================================================================

%% IL model
if any(strcmp(models,'IL')) || any(strcmp(models,'ALL'))
    results.IL.allfit = nan(2+4, nFit); % 2 parameters + neglhtrial, neglh, AIC, BIC
    % initialize
    opt.x0 = [40, 2]; % initial guess
    opt.PLB = [0, 0];
    opt.PUB = [maxJ, 10];
    opt.LB = [0, 0];
    opt.UB = [maxJ, 20];
    opt.options = bads('defaults');
    opt.options.MaxIter = 'maxJ*nvars';
    % initial guess
    opt.x0 = [opt.PLB(1):(opt.PUB(1)-opt.PLB(1))/(nFit-1):opt.PUB(1);
        opt.PLB(2):(opt.PUB(2)-opt.PLB(2))/(nFit-1):opt.PUB(2);
        ];
    opt.x0=opt.x0';
    
    % fit data
    for iFit = 1:nFit
        fprintf('IL model, Fit: %d\n',iFit)
        x0_tmp = opt.x0(iFit,:);
        [fitparams,neglhtrial, neglh,AIC,BIC] = fit_IL_model(data.N, data.probe, data.resp, x0_tmp, opt);
        results.IL.allfit(:,iFit) =[fitparams neglhtrial neglh AIC BIC];
    end
    
    % find the best params
    [~,idx] = min(results.IL.allfit(end-2,:));
    results.IL.best = results.IL.allfit(:,idx); 
end

%% MIX model
if any(strcmp(models,'MIX')) || any(strcmp(models,'ALL'))
    results.MIX.allfit = nan(nSetSize + 1 + 4, nFit); %  parameters + neglhtrial neglh AIC, BIC
    % initialize
    opt.x0 = [40*ones(1, nSetSize), 2]; % initial guess
    opt.PLB = [zeros(1, nSetSize), 0];
    opt.PUB = [maxJ*ones(1, nSetSize), 10];
    opt.LB = opt.PLB;
    opt.UB = [maxJ*ones(1, nSetSize), 20];
    opt.options = bads('defaults');
    opt.options.MaxIter = 'maxJ*nvars';
    % initial guess
    x0 = nan(nSetSize+1, nFit);
    for i=1:nSetSize+1
        x0(i,:) = opt.PLB(i):(opt.PUB(i)-opt.PLB(i))/(nFit-1):opt.PUB(i);
    end
    x0=x0';
    % fit data
    for iFit = 1:nFit
        fprintf('MIX model, Fit: %d \n',iFit)
        x0_tmp = x0(iFit,:);
        [fitparams, neglhtrial, neglh,AIC,BIC] = fit_MIX_model(data.N, data.resp, data.probe,x0_tmp, opt);
        results.MIX.allfit(:,iFit) =[fitparams neglhtrial neglh AIC BIC];
    end
    
    % find the best params
    [~,idx] = min(results.MIX.allfit(end-2,:));
    results.MIX.best = results.MIX.allfit(:,idx); 
end


%% SA model
if any(strcmp(models,'SA')) || any(strcmp(models,'ALL'))
    results.SA.allfit = nan(3 + 4, nFit); % 3 parameters + neglhtrial neglh, AIC, BIC
    % initialize
    opt.PLB = [0, 0, 0];
    opt.PUB = [10, 8, maxJ];
    opt.LB = [0, 0, 0];
    opt.UB = [10, 8, maxJ];
    opt.options = bads('defaults');
    opt.options.MaxIter = 'maxJ*nvars';
    % initial guess
    x0 = [opt.PLB(1):(opt.PUB(1)-opt.PLB(1))/(nFit-1):opt.PUB(1);
        opt.PLB(2):(opt.PUB(2)-opt.PLB(2))/(nFit-1):opt.PUB(2);
        opt.PLB(3):(opt.PUB(3)-opt.PLB(3))/(nFit-1):opt.PUB(3);
        ];
    x0=x0';
    
    % fit data
    for iFit = 1:nFit
        fprintf('SA model, Fit: %d \n',iFit);
        x0_tmp = x0(iFit, :);
        [fitparams,neglhtrial, neglh,AIC,BIC] = fit_SA_model(data.N, data.probe, data.resp, x0_tmp, opt);
        results.SA.allfit(:,iFit) =[fitparams neglhtrial neglh AIC BIC];
    end
    
    % find the best params
    [~,idx] = min(results.SA.allfit(end-2,:));
    results.SA.best = results.SA.allfit(:,idx); 
end

%% EP model
if any(strcmp(models,'EP')) || any(strcmp(models,'ALL'))
    results.EP.allfit = nan(3 + 4, nFit); % 3 parameters + neglhtrial, neglh, AIC, BIC
    opt.PLB = [0, 0, 0];
    opt.PUB = [maxJ, 10, maxJ];
    opt.LB = [0, 0, 0];
    opt.UB = [maxJ, 10, maxJ];
    opt.options = bads('defaults');
    opt.options.MaxIter = 'maxJ*nvars';
    x0 = [opt.PLB(1)+eps:(opt.PUB(1)-opt.PLB(1))/(nFit-1):opt.PUB(1)-eps;
        opt.PLB(2)+eps:(opt.PUB(2)-opt.PLB(2))/(nFit-1):opt.PUB(2)-eps;
        opt.PLB(3)+eps:(opt.PUB(3)-opt.PLB(3))/(nFit-1):opt.PUB(3)-eps];
    
    x0=x0';  %x0,  nFit x 3 params
    
    % fit data
    for iFit = 1:nFit
        fprintf('EP model, Fit: %d \n',iFit);
        x0_tmp = x0(iFit, :);
        [fitparams, neglhtrial, neglh,AIC,BIC] = fit_EP_model(data.N, data.probe, data.resp, x0_tmp);
        results.EP.allfit(:,iFit) =[fitparams neglhtrial neglh AIC BIC];
    end
    
    % find the best params
    [~,idx] = min(results.EP.allfit(end-2,:));
    results.EP.best = results.EP.allfit(:,idx); 
end


%% VP model
if any(strcmp(models,'VP')) || any(strcmp(models,'ALL'))
    results.VP.allfit = nan(4 + 4, nFit); % 4 parameters + neglhtrial neglh AIC, BIC
    opt.PLB = [0, 0, 0, 0];
    opt.PUB = [maxJ, 10, maxJ, maxJ];
    opt.LB = [0, 0, 0, 0];
    opt.UB = [maxJ, 20, maxJ, maxJ];
    opt.options = bads('defaults');
    opt.options.MaxIter = 'maxJ*nvars';
    x0 = [opt.PLB(1):(opt.PUB(1)-opt.PLB(1))/(nFit-1):opt.PUB(1);
        opt.PLB(2):(opt.PUB(2)-opt.PLB(2))/(nFit-1):opt.PUB(2);
        opt.PLB(3):(opt.PUB(3)-opt.PLB(3))/(nFit-1):opt.PUB(3);
        opt.PLB(4):(opt.PUB(4)-opt.PLB(4))/(nFit-1):opt.PUB(4)];
    x0=x0';  %x0,  nFit x params
    
    % fit data
    for iFit = 1:nFit
        fprintf('VP model, Fit: %d \n',iFit);
        x0_tmp = x0(iFit, :);
        [fitparams,neglhtrial, neglh,AIC,BIC] = fit_VP_model(data.N, data.probe, data.resp, x0_tmp, opt);
        results.VP.allfit(:,iFit) =[fitparams neglhtrial neglh AIC BIC];
    end
    
    % find the best params
    [~,idx] = min(results.VP.allfit(end-2,:));
    results.VP.best = results.VP.allfit(:,idx); 
    
end


%% VPcap model
if any(strcmp(models,'VPCAP')) || any(strcmp(models,'ALL'))
    results.VPCAP.allfit = nan(5 + 4, nFit); % 5 parameters + neglhtrial neglh AIC, BIC
    
    opt.PLB = [0, 0, 0, 0, 1];
    opt.PUB = [maxJ, 10, maxJ, maxJ, 10];
    opt.LB = [0, 0, 0, 0, 1];
    opt.UB = [maxJ, 20, maxJ, maxJ, 20];
    opt.options = bads('defaults');
    opt.options.MaxIter = 'maxJ*nvars';
    x0 = [opt.PLB(1):(opt.PUB(1)-opt.PLB(1))/(nFit-1):opt.PUB(1);
        opt.PLB(2):(opt.PUB(2)-opt.PLB(2))/(nFit-1):opt.PUB(2);
        opt.PLB(3):(opt.PUB(3)-opt.PLB(3))/(nFit-1):opt.PUB(3);
        opt.PLB(4):(opt.PUB(4)-opt.PLB(4))/(nFit-1):opt.PUB(4);
        opt.PLB(5):(opt.PUB(5)-opt.PLB(5))/(nFit-1):opt.PUB(5)];
    x0=x0';  %x0,  nFit x params
    % fit data
    for iFit = 1:nFit
        fprintf('VPcap model, Fit: %d \n',iFit);
        x0_tmp = x0(iFit, :);
        [fitparams, neglhtrial, neglh,AIC,BIC] = fit_VPcap_model(data.N, data.probe, data.resp, x0_tmp, opt);
        results.VPCAP.allfit(:,iFit) =[fitparams neglhtrial neglh AIC BIC];
    end
    
    % find the best params
    [~,idx] = min(results.VPCAP.allfit(end-2,:));
    results.VPCAP.best = results.VPCAP.allfit(:,idx); 
end


%% cosSA model (fitting this model takes a long time)
if any(strcmp(models,'COSSA')) || any(strcmp(models,'ALL'))
    results.COSSA.allfit = nan(4 + 4, nFit); % 4 parameters + neglhtrial neglh, AIC, BIC
    % initialize
    opt.PLB = [0, 0, 0, 0];
    opt.PUB = [10, 8, maxJ, 2*pi];
    opt.LB = [0, 0, 0, 0];
    opt.UB = [10, 8, maxJ, 2*pi];
    opt.options = bads('defaults');
    opt.options.MaxIter = 'maxJ*nvars';
    % initial guess
    x0 = [opt.PLB(1):(opt.PUB(1)-opt.PLB(1))/(nFit-1):opt.PUB(1);
        opt.PLB(2):(opt.PUB(2)-opt.PLB(2))/(nFit-1):opt.PUB(2);
        opt.PLB(3):(opt.PUB(3)-opt.PLB(3))/(nFit-1):opt.PUB(3);
        opt.PLB(4):(opt.PUB(4)-opt.PLB(4))/(nFit-1):opt.PUB(4);
        ];
    x0=x0';
    % fit data
    for iFit = 1:nFit
        fprintf('cosSA model, Fit: %d \n',iFit);
        x0_tmp = x0(iFit, :);
        [fitparams,neglhtrial, neglh, AIC, BIC] = fit_cosSA_model(data.N, data.probe, data.resp, x0_tmp, opt);
        results.COSSA.allfit(:,iFit) =[fitparams neglhtrial neglh AIC BIC];
    end
    
    % find the best params
    [~,idx] = min(results.COSSA.allfit(end-2,:));
    results.COSSA.best = results.COSSA.allfit(:,idx); 
end

%% We do some simple analysis
% find the best-fitting model
AIC = zeros(1, results.nModeltoFit);
BIC = zeros(1, results.nModeltoFit);
for iModel = 1:results.nModeltoFit
    tmp = getfield(results, results.models{iModel});
    AIC(iModel) = tmp.best(end-1); % We extract the AIC; 
    BIC(iModel) = tmp.best(end); % We extract the BIC; 
end
[~,idx] = min(AIC);
results.bestModel_AIC = results.models{idx};
[~,idx] = min(BIC);
results.bestModel_BIC = results.models{idx};

results.AIC = AIC;
results.BIC = BIC;

%% save results and clean up
filename = strcat(prefix,'results_',datestr(now,'yymmddHHMM'),'.mat');
if exist(filename,'file')
    error('data file name exists');
end
% save(filename);

rmpath(genpath('./'));







