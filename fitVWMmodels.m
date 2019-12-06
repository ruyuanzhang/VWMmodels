function results = fitVWMmodels(data, models, prefix)
% function results = fitVWMmodels(data, models, prefix)
% 
% <data>: a structure
%   data.probe: 1-180 deg, probe color
%   data.resp: 1-180 deg, response color
%   data.error: -90~89 deg, circular error
%   data.N: a vector of integer, set-size level
%   data.squares: a cell vector, each element is a vector containing color index of targets other than probe.  
% <models>: a cell vector of model name to fit, or 'all' to fit all
%   avaliable models.
% <prefix>: prefix of the result file.
%
% Here we tried to fit a set of VWM models, including
%
% Limited-item(IL) model:
% Mixture (MIX) model:
% Slot-plus-averaging with cosine fluctuation (cosSA) model:
% Slot-plus-averaging(SA) model:
% Equal precision(EP):
% Variable-precision (VP) model:
% Variable-precision-with-capacity (Vpcap) model:
%
%
% History
%   20191116 TY revised cosSA model and added swap variants
%   20191117 TY added AICc and converted the information criteria into
%       posterior probability (i.e. p(Model | Data, Model Space))
%   20191029 RZ let the model output likelihood per trial to compensate 
%       different trials across subjects and groups
%   20191027 RZ wrote 1st version

if notDefined('prefix')
    prefix='';
end

%% Settings
if ~iscell(models)
    models = {models};
end
models = upper(models); % convert all model names to upper case
if any(strcmp(models,'ALL')) % fit all models
    models = {'IL', 'MIX', 'COSSA', 'SA', 'EP', 'VP', 'VPCAP'};
end

setSize = unique(data.N);
nSetSize = length(unique(data.N));
nFit = 20; % fit how many times for each model
maxJ1bar = 700; % upper bound of J1bar parameters 

results.setSize = setSize;
results.models = models;
results.nFit = nFit; % for each model, we random initialize parameters and fit <nFit> times 
results.data = data;
results.allposibModels = {'IL', 'MIX', 'cosSA', 'SA', 'EP', 'VP', 'VPcap',...
    'IL_swap','MIX_swap','cosSA_swap','SA_swap','EP_swap','VP_swap','VPcap_swap'};
results.allposibModelParams = {
    {'Kr','cap', 'neglh', 'neglhtrial', 'AIC', 'AICc', 'BIC'}, ... % IL
    {'Kr for different set size levels','cap', 'neglh', 'neglhtrial','AIC', 'AICc', 'BIC'}, ... % MIX
    {'J1','K', 'Kr', 'neglh', 'neglhtrial','AIC', 'AICc', 'BIC'}, ... % SA
    {'J1bar','power', 'Kr', 'neglh', 'neglhtrial', 'AIC', 'AICc', 'BIC'}, ... % EP
    {'J1bar','power', 'tau', 'Kr', 'neglh', 'neglhtrial', 'AIC', 'AICc', 'BIC'}, ... % VP
    {'J1bar','power', 'tau', 'Kr', 'K','neglh', 'neglhtrial', 'AIC', 'AICc', 'BIC'}, ... % Vpcap
    {'Jm','K', 'Jf', 'muf', 'neglh', 'neglhtrial', 'AIC', 'AICc', 'BIC'}, ... % cosSA
    {'Kr','cap', 's', 'neglh', 'neglhtrial', 'AIC', 'AICc', 'BIC'}, ... % IL + swap
    {'Kr for different set size levels','cap', 's', 'neglh', 'neglhtrial', 'AIC', 'AICc', 'BIC'}, ... % MIX + swap
    {'J1','K', 'Kr', 's', 'neglh', 'neglhtrial','AIC', 'AICc', 'BIC'}, ... % SA + swap
    {'J1bar','power', 'Kr', 's', 'neglh', 'neglhtrial', 'AIC', 'AICc', 'BIC'}, ... % EP + swap
    {'J1bar','power', 'tau', 'Kr', 's', 'neglh', 'neglhtrial', 'AIC', 'AICc', 'BIC'}, ... % VP + swap
    {'J1bar','power', 'tau', 'Kr', 'K', 's', 'neglh', 'neglhtrial', 'AIC', 'AICc', 'BIC'}, ... % Vpcap + swap
    {'Jm','K', 'Jf', 'muf', 's', 'neglh', 'neglhtrial', 'AIC', 'AICc', 'BIC'}, ... % cosSA + swap
    };
results.nModeltoFit = numel(models);


%% add BADS optimization toolbox
addpath(genpath('./'));

%% ========================================================================
% ===================== Now fit models ====================================
%% ========================================================================

%% IL model
if any(strcmp(models,'IL'))
    results.IL.allfit = nan(2+4, nFit); % 2 parameters + neglhtrial, neglh, AIC, BIC
    % initialize
    opt.x0 = [40, 2]; % initial guess
    opt.PLB = [0, 0];
    opt.PUB = [maxJ1bar, 10];
    opt.LB = [0, 0];
    opt.UB = [maxJ1bar, 20];
    opt.options = bads('defaults');
    opt.options.MaxIter = 'maxJ1bar*nvars';
    % initial guess
    opt.x0 = [opt.PLB(1):(opt.PUB(1)-opt.PLB(1))/(nFit-1):opt.PUB(1);
        opt.PLB(2):(opt.PUB(2)-opt.PLB(2))/(nFit-1):opt.PUB(2);
        ];
    opt.x0=opt.x0';
    
    % fit data
    for iFit = 1:nFit
        fprintf('IL model, Fit: %d\n',iFit)
        x0_tmp = opt.x0(iFit,:);
        [fitparams,neglhtrial, neglh,AIC,AICc,BIC] = fit_IL_model(data.N, data.probe, data.nt, data.resp, x0_tmp, opt);
        results.IL.allfit(:,iFit) =[fitparams neglhtrial neglh AIC AICc BIC];
    end
    
    % find the best params
    [~,idx] = min(results.IL.allfit(end-3,:));
    results.IL.best = results.IL.allfit(:,idx); 
end

%% IL model with swap rate
if any(strcmp(models,'IL_swap'))
    results.IL.allfit = nan(3+4, nFit); % 3 parameters + neglhtrial, neglh, AIC, BIC
    % initialize
    opt.x0 = [40, 2, 0.1]; % initial guess
    opt.PLB = [0, 0, 0];
    opt.PUB = [maxJ, 10, 0.9];
    opt.LB = [0, 0, 0];
    opt.UB = [maxJ, 20, 1];
    opt.options = bads('defaults');
    opt.options.MaxIter = 'maxJ*nvars';
    % initial guess
    opt.x0 = [opt.PLB(1):(opt.PUB(1)-opt.PLB(1))/(nFit-1):opt.PUB(1);
        opt.PLB(2):(opt.PUB(2)-opt.PLB(2))/(nFit-1):opt.PUB(2);
        opt.PLB(3):(opt.PUB(2)-opt.PLB(3))/(nFit-1):opt.PUB(3);
        ];
    opt.x0=opt.x0';
    
    % fit data
    for iFit = 1:nFit
        fprintf('IL model + swap, Fit: %d\n',iFit)
        x0_tmp = opt.x0(iFit,:);
        [fitparams,neglhtrial, neglh,AIC,AICc,BIC] = fit_ILswap_model(data.N, data.probe, data.nt, data.resp, x0_tmp, opt);
        results.ILswap.allfit(:,iFit) =[fitparams neglhtrial neglh AIC AICc BIC];
    end
    
    % find the best params
    [~,idx] = min(results.ILswap.allfit(end-3,:));
    results.ILswap.best = results.ILswap.allfit(:,idx); 
end

%% MIX model
if any(strcmp(models,'MIX'))
    results.MIX.allfit = nan(nSetSize + 1 + 4, nFit); %  parameters + neglhtrial neglh AIC, BIC
    % initialize
    opt.x0 = [40*ones(1, nSetSize), 2]; % initial guess
    opt.PLB = [zeros(1, nSetSize), 0];
    opt.PUB = [maxJ1bar*ones(1, nSetSize), 10];
    opt.LB = opt.PLB;
    opt.UB = [maxJ1bar*ones(1, nSetSize), 20];
    opt.options = bads('defaults');
    opt.options.MaxIter = 'maxJ1bar*nvars';
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
        [fitparams, neglhtrial, neglh,AIC,AICc,BIC] = fit_MIX_model(data.N, data.probe, data.resp,x0_tmp, opt);
        results.MIX.allfit(:,iFit) =[fitparams neglhtrial neglh AIC AICc BIC];
    end
    
    % find the best params
    [~,idx] = min(results.MIX.allfit(end-3,:));
    results.MIX.best = results.MIX.allfit(:,idx); 
end

%% MIX model with swap rate
if any(strcmp(models,'MIX_swap'))
    results.MIX.allfit = nan(nSetSize + 2 + 4, nFit); %  parameters + neglhtrial neglh AIC, BIC
    % initialize
    opt.x0 = [40*ones(1, nSetSize), 2]; % initial guess
    opt.PLB = [zeros(1, nSetSize), 0, 0];
    opt.PUB = [maxJ*ones(1, nSetSize), 10, 0.9];
    opt.LB = opt.PLB;
    opt.UB = [maxJ*ones(1, nSetSize), 20, 1];
    opt.options = bads('defaults');
    opt.options.MaxIter = 'maxJ*nvars';
    % initial guess
    x0 = nan(nSetSize+1, nFit);
    for i=1:nSetSize+2
        x0(i,:) = opt.PLB(i):(opt.PUB(i)-opt.PLB(i))/(nFit-1):opt.PUB(i);
    end
    x0=x0';
    % fit data
    for iFit = 1:nFit
        fprintf('MIX model + swap, Fit: %d \n',iFit)
        x0_tmp = x0(iFit,:);
        [fitparams, neglhtrial, neglh,AIC,AICc,BIC] = fit_MIXswap_model(data.N, data.probe, data.nt, data.resp,x0_tmp, opt);
        results.MIXswap.allfit(:,iFit) =[fitparams neglhtrial neglh AIC AICc BIC];
    end
    
    % find the best params
    [~,idx] = min(results.MIXswap.allfit(end-3,:));
    results.MIXswap.best = results.MIXswap.allfit(:,idx); 
end

%% SA model
if any(strcmp(models,'SA'))
    results.SA.allfit = nan(3 + 4, nFit); % 3 parameters + neglhtrial neglh, AIC, BIC
    % initialize
    opt.PLB = [0, 0, 0];
    opt.PUB = [10, 8, maxJ1bar];
    opt.LB = [0, 0, 0];
    opt.UB = [10, 8, maxJ1bar];
    opt.options = bads('defaults');
    opt.options.MaxIter = 'maxJ1bar*nvars';
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
        [fitparams,neglhtrial, neglh,AIC,AICc,BIC] = fit_SA_model(data.N, data.probe, data.resp, x0_tmp, opt);
        results.SA.allfit(:,iFit) =[fitparams neglhtrial neglh AIC AICc BIC];
    end
    
    % find the best params
    [~,idx] = min(results.SA.allfit(end-3,:));
    results.SA.best = results.SA.allfit(:,idx); 
end

%% SA model with swap rate
if any(strcmp(models,'SA_swap'))
    results.SA.allfit = nan(4 + 4, nFit); % 4 parameters + neglhtrial neglh, AIC, BIC
    % initialize
    opt.PLB = [0, 0, 0, 0];
    opt.PUB = [10, 8, maxJ, 0.9];
    opt.LB = [0, 0, 0, 0];
    opt.UB = [10, 8, maxJ, 1];
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
        fprintf('SA model + swap, Fit: %d \n',iFit);
        x0_tmp = x0(iFit, :);
        [fitparams,neglhtrial, neglh,AIC,AICc,BIC] = fit_SAswap_model(data.N, data.probe, data.nt, data.resp, x0_tmp, opt);
        results.SAswap.allfit(:,iFit) =[fitparams neglhtrial neglh AIC AICc BIC];
    end
    
    % find the best params
    [~,idx] = min(results.SAswap.allfit(end-3,:));
    results.SAswap.best = results.SAswap.allfit(:,idx); 
end

%% EP model
if any(strcmp(models,'EP'))
    results.EP.allfit = nan(3 + 4, nFit); % 3 parameters + neglhtrial, neglh, AIC, BIC
    opt.PLB = [0, 0, 0];
    opt.PUB = [maxJ1bar, 10, maxJ1bar];
    opt.LB = [0, 0, 0];
    opt.UB = [maxJ1bar, 10, maxJ1bar];
    opt.options = bads('defaults');
    opt.options.MaxIter = 'maxJ1bar*nvars';
    x0 = [opt.PLB(1)+eps:(opt.PUB(1)-opt.PLB(1))/(nFit-1):opt.PUB(1)-eps;
        opt.PLB(2)+eps:(opt.PUB(2)-opt.PLB(2))/(nFit-1):opt.PUB(2)-eps;
        opt.PLB(3)+eps:(opt.PUB(3)-opt.PLB(3))/(nFit-1):opt.PUB(3)-eps];
    
    x0=x0';  %x0,  nFit x 3 params
    
    % fit data
    for iFit = 1:nFit
        fprintf('EP model, Fit: %d \n',iFit);
        x0_tmp = x0(iFit, :);
        [fitparams, neglhtrial, neglh,AIC,AICc,BIC] = fit_EP_model(data.N, data.probe, data.resp, x0_tmp);
        results.EP.allfit(:,iFit) =[fitparams neglhtrial neglh AIC AICc BIC];
    end
    
    % find the best params
    [~,idx] = min(results.EP.allfit(end-3,:));
    results.EP.best = results.EP.allfit(:,idx); 
end

%% EP model with swap rate
if any(strcmp(models,'EP_swap'))
    results.EP.allfit = nan(4 + 4, nFit); % 4 parameters + neglhtrial, neglh, AIC, BIC
    opt.PLB = [0, 0, 0, 0];
    opt.PUB = [maxJ, 10, maxJ, 0.9];
    opt.LB = [0, 0, 0, 0];
    opt.UB = [maxJ, 10, maxJ, 1];
    opt.options = bads('defaults');
    opt.options.MaxIter = 'maxJ*nvars';
    x0 = [opt.PLB(1)+eps:(opt.PUB(1)-opt.PLB(1))/(nFit-1):opt.PUB(1)-eps;
        opt.PLB(2)+eps:(opt.PUB(2)-opt.PLB(2))/(nFit-1):opt.PUB(2)-eps;
        opt.PLB(3)+eps:(opt.PUB(3)-opt.PLB(3))/(nFit-1):opt.PUB(3)-eps;
        opt.PLB(4)+eps:(opt.PUB(4)-opt.PLB(4))/(nFit-1):opt.PUB(4)-eps];
    
    x0=x0';  %x0,  nFit x 4 params
    
    % fit data
    for iFit = 1:nFit
        fprintf('EP model + swap, Fit: %d \n',iFit);
        x0_tmp = x0(iFit, :);
        [fitparams, neglhtrial, neglh,AIC,AICc,BIC] = fit_EPswap_model(data.N, data.probe, data.nt, data.resp, x0_tmp);
        results.EPswap.allfit(:,iFit) =[fitparams neglhtrial neglh AIC AICc BIC];
    end
    
    % find the best params
    [~,idx] = min(results.EPswap.allfit(end-3,:));
    results.EPswap.best = results.EPswap.allfit(:,idx); 
end

%% VP model
if any(strcmp(models,'VP'))
    results.VP.allfit = nan(4 + 4, nFit); % 4 parameters + neglhtrial neglh AIC, BIC
    opt.PLB = [0, 0, 0, 0];
    opt.PUB = [maxJ1bar, 10, maxJ1bar, maxJ1bar];
    opt.LB = [0, 0, 0, 0];
    opt.UB = [maxJ1bar, 20, maxJ1bar, maxJ1bar];
    opt.options = bads('defaults');
    opt.options.MaxIter = 'maxJ1bar*nvars';
    x0 = [opt.PLB(1):(opt.PUB(1)-opt.PLB(1))/(nFit-1):opt.PUB(1);
        opt.PLB(2):(opt.PUB(2)-opt.PLB(2))/(nFit-1):opt.PUB(2);
        opt.PLB(3):(opt.PUB(3)-opt.PLB(3))/(nFit-1):opt.PUB(3);
        opt.PLB(4):(opt.PUB(4)-opt.PLB(4))/(nFit-1):opt.PUB(4)];
    x0=x0';  %x0,  nFit x params
    
    % fit data
    for iFit = 1:nFit
        fprintf('VP model, Fit: %d \n',iFit);
        x0_tmp = x0(iFit, :);
        [fitparams,neglhtrial, neglh,AIC,AICc,BIC] = fit_VP_model(data.N, data.probe, data.resp, x0_tmp, opt);
        results.VP.allfit(:,iFit) =[fitparams neglhtrial neglh AIC AICc BIC];
    end
    
    % find the best params
    [~,idx] = min(results.VP.allfit(end-3,:));
    results.VP.best = results.VP.allfit(:,idx); 
    
end

%% VP model with swap rate
if any(strcmp(models,'VP_swap')) || any(strcmp(models,'ALL'))
    results.VP.allfit = nan(5 + 4, nFit); % 5 parameters + neglhtrial neglh AIC, BIC
    opt.PLB = [0, 0, 0, 0, 0];
    opt.PUB = [maxJ, 10, maxJ, maxJ, .9];
    opt.LB = [0, 0, 0, 0, 0];
    opt.UB = [maxJ, 20, maxJ, maxJ, 1];
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
        fprintf('VP model + swap, Fit: %d \n',iFit);
        x0_tmp = x0(iFit, :);
        [fitparams,neglhtrial, neglh,AIC,AICc,BIC] = fit_VPswap_model(data.N, data.probe, data.nt, data.resp, x0_tmp, opt);
        results.VPswap.allfit(:,iFit) =[fitparams neglhtrial neglh AIC AICc BIC];
    end
    
    % find the best params
    [~,idx] = min(results.VPswap.allfit(end-3,:));
    results.VPswap.best = results.VPswap.allfit(:,idx); 
    
end

%% VPcap model
if any(strcmp(models,'VPCAP')) || any(strcmp(models,'ALL'))
    results.VPCAP.allfit = nan(5 + 5, nFit); % 5 parameters + neglhtrial neglh AIC, AICc, BIC
    
    opt.PLB = [0, 0, 0, 0, 1];
    opt.PUB = [maxJ1bar, 10, maxJ1bar, maxJ1bar, 10];
    opt.LB = [0, 0, 0, 0, 1];
    opt.UB = [maxJ1bar, 20, maxJ1bar, maxJ1bar, 20];
    opt.options = bads('defaults');
    opt.options.MaxIter = 'maxJ1bar*nvars';
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
        [fitparams, neglhtrial, neglh,AIC,AICc,BIC] = fit_VPcap_model(data.N, data.probe, data.resp, x0_tmp, opt);
        results.VPCAP.allfit(:,iFit) =[fitparams neglhtrial neglh AIC AICc BIC];
    end
    
    % find the best params
    [~,idx] = min(results.VPCAP.allfit(end-3,:));
    results.VPCAP.best = results.VPCAP.allfit(:,idx); 
end

%% VPcap model with swap rate
if any(strcmp(models,'VPCAP_swap')) || any(strcmp(models,'ALL'))
    results.VPCAP.allfit = nan(6 + 5, nFit); % 6 parameters + neglhtrial neglh AIC, AICc, BIC
    
    opt.PLB = [0, 0, 0, 0, 1, 0];
    opt.PUB = [maxJ, 10, maxJ, maxJ, 10, 0.9];
    opt.LB = [0, 0, 0, 0, 1, 0];
    opt.UB = [maxJ, 20, maxJ, maxJ, 20, 1];
    opt.options = bads('defaults');
    opt.options.MaxIter = 'maxJ*nvars';
    x0 = [opt.PLB(1):(opt.PUB(1)-opt.PLB(1))/(nFit-1):opt.PUB(1);
        opt.PLB(2):(opt.PUB(2)-opt.PLB(2))/(nFit-1):opt.PUB(2);
        opt.PLB(3):(opt.PUB(3)-opt.PLB(3))/(nFit-1):opt.PUB(3);
        opt.PLB(4):(opt.PUB(4)-opt.PLB(4))/(nFit-1):opt.PUB(4);
        opt.PLB(5):(opt.PUB(5)-opt.PLB(5))/(nFit-1):opt.PUB(5);
        opt.PLB(6):(opt.PUB(6)-opt.PLB(6))/(nFit-1):opt.PUB(6)];
    x0=x0';  %x0,  nFit x params
    % fit data
    for iFit = 1:nFit
        fprintf('VPcap model + swap, Fit: %d \n',iFit);
        x0_tmp = x0(iFit, :);
        [fitparams, neglhtrial, neglh,AIC,AICc,BIC] = fit_VPcapSwap_model(data.N, data.probe, data.resp, x0_tmp, opt);
        results.VPCAPswap.allfit(:,iFit) =[fitparams neglhtrial neglh AIC AICc BIC];
    end
    
    % find the best params
    [~,idx] = min(results.VPCAPswap.allfit(end-3,:));
    results.VPCAPswap.best = results.VPCAPswap.allfit(:,idx); 
end

%% cosSA model (fitting this model relatively takes a long time)
if any(strcmp(models,'COSSA')) || any(strcmp(models,'ALL'))
    results.COSSA.allfit = nan(4 + 5, nFit); % 4 parameters + neglhtrial neglh, AIC, AICc, BIC
    % initialize
    opt.PLB = [0, 0, 0, 0];
    opt.PUB = [10, 8, maxJ1bar, 2*pi];
    opt.LB = [0, 0, 0, 0];
    opt.UB = [10, 8, maxJ1bar, 2*pi];
    opt.options = bads('defaults');
    opt.options.MaxIter = 'maxJ1bar*nvars';
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
        [fitparams,neglhtrial, neglh, AIC, AICc, BIC] = fit_cosSA_model(data.N, data.probe, data.resp, x0_tmp, opt);
        results.COSSA.allfit(:,iFit) =[fitparams neglhtrial neglh AIC AICc BIC];
    end
    
    % find the best params
    [~,idx] = min(results.COSSA.allfit(end-3,:));
    results.COSSA.best = results.COSSA.allfit(:,idx); 
end

%% cosSA model with swap rate (fitting this model relatively takes a long time)
if any(strcmp(models,'COSSA_swap')) || any(strcmp(models,'ALL'))
    results.COSSA.allfit = nan(5 + 5, nFit); % 4 parameters + neglhtrial neglh, AIC, AICc, BIC
    % initialize
    opt.PLB = [0, 0, 0, 0, 0];
    opt.PUB = [10, 8, maxJ, 2*pi, .9];
    opt.LB = [0, 0, 0, 0, 0];
    opt.UB = [10, 8, maxJ, 2*pi, 1];
    opt.options = bads('defaults');
    opt.options.MaxIter = 'maxJ*nvars';
    % initial guess
    x0 = [opt.PLB(1):(opt.PUB(1)-opt.PLB(1))/(nFit-1):opt.PUB(1);
        opt.PLB(2):(opt.PUB(2)-opt.PLB(2))/(nFit-1):opt.PUB(2);
        opt.PLB(3):(opt.PUB(3)-opt.PLB(3))/(nFit-1):opt.PUB(3);
        opt.PLB(4):(opt.PUB(4)-opt.PLB(4))/(nFit-1):opt.PUB(4);
        opt.PLB(5):(opt.PUB(5)-opt.PLB(5))/(nFit-1):opt.PUB(5);
        ];
    x0=x0';
    % fit data
    for iFit = 1:nFit
        fprintf('cosSA model + swap, Fit: %d \n',iFit);
        x0_tmp = x0(iFit, :);
        [fitparams,neglhtrial, neglh, AIC, AICc, BIC] = fit_cosSAswap_model(data.N, data.probe, data.resp, x0_tmp, opt);
        results.COSSA.allfit(:,iFit) =[fitparams neglhtrial neglh AIC AICc BIC];
    end
    
    % find the best params
    [~,idx] = min(results.COSSAswap.allfit(end-3,:));
    results.COSSAswap.best = results.COSSAswap.allfit(:,idx); 
end

%% We do some simple analysis
% find the best-fitting model
AIC = zeros(1, results.nModeltoFit);
AICc = zeros(1, results.nModeltoFit);
BIC = zeros(1, results.nModeltoFit);
for iModel = 1:results.nModeltoFit
    tmp = getfield(results, results.models{iModel});
    AIC(iModel) = tmp.best(end-2); % We extract the AIC;
    AICc(iModel) = tmp.best(end-1); % We extract the AICc;
    BIC(iModel) = tmp.best(end); % We extract the BIC; 
end
[~,idx] = min(AIC);
results.bestModel_AIC = results.models{idx};
[~,idx] = min(AICc);
results.bestModel_AICc = results.models{idx};
[~,idx] = min(BIC);
results.bestModel_BIC = results.models{idx};

% Posterior probability
results.PP_AIC=IC2PP(AIC);
results.PP_BIC=IC2PP(BIC);
results.PP_AICc=IC2PP(AICc);

results.AIC = AIC;
results.AICc=AICc;
results.BIC = BIC;

%% save results and clean up
filename = strcat(prefix,'results_',datestr(now,'yymmddHHMM'),'.mat');
if exist(filename,'file')
    error('data file name exists');
end
% save(filename);

rmpath(genpath('./'));

end

%% Calculate posterior probability
function PP=IC2PP(IC)
    BestModel=min(IC);
    IC=IC-BestModel;
    PP=exp(-IC/2)/sum(exp(-IC/2));
end

