function results = fitVWMmodels(data, models, prefix)
% function results = fitVWMmodels(data, models, prefix)
% 
% <data>: a structure
%   data.probe: 1-180 deg, probe color
%   data.resp: 1-180 deg, response color
%   data.error: -90~89 deg, circular error
%   data.N: a vector of integer, set-size level
%   data.distr: a cell vector, each element is a vector containing color index of targets other than probe.
%   data.distrError: a cell vector, each element contains the response
%       errors with respect to other distractors
%
% <models>: a cell vector of model name to fit, or 'all' to fit all
%   avaliable models.
%
% <prefix>: prefix of the result file.
%
% Here we try to fit a set of VWM models, including
%
% Limited-item(IL) model:
% Mixture (MIX) model:
% Slot-plus-averaging with cosine fluctuation (cosSA) model:
% Slot-plus-averaging(SA) model:
% Equal precision(EP):
% Variable-precision (VP) model:
% Variable-precision-with-capacity (Vpcap) model:
%
% We also include a version with non-target swap errors.
%
% History
%   2021-11-25 RZ massively reorganized the data and add the versions with
%       non-target swap error.
%   20191217 RZ seperate swap variants into a different branch and focus on
%       non swap version
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
    models = {
        'IL', 'MIX', 'COSSA', 'SA', 'EP', 'VP', 'VPCAP',...
        'ILNT', 'MIXNT', 'COSSANT', 'SANT', 'EPNT', 'VPNT', 'VPCAPNT'        
        };
end

setSize = unique(data.N); % set size numbers
nSetSize = length(unique(data.N)); % levels of sit sizes
nFit = 20; % fit how many times for each model

results.setSize = setSize;
results.models = models;
results.nFit = nFit; % for each model, we random initialize parameters and fit <nFit> times 
%results.data = data;
results.nModeltoFit = numel(models);


%% add BADS optimization toolbox
addpath(genpath('./'));

%% ========================================================================
% ===================== Now fit models ====================================
%% ========================================================================

%% IL model
results.fitResults = cell(1, results.nModeltoFit);
for iModel=1:results.nModeltoFit % loop model
    models{iModel}
    if ismember(models{iModel}, results.models)
        switch models{iModel}
            case 'IL'
                c = VWM_IL_config;
            case 'SA'
                c = VWM_SA_config;
            case 'MIX'
                c = VWM_MIX_config(nSetSize);
            case 'EP'
                c = VWM_EP_config;
            case 'VP'
                c = VWM_VP_config;
            case 'VPCAP'
                c = VWM_VPcap_config;
            case 'COSSA'
                c = VWM_cosSA_config;
            case 'ILNT'
                c = VWM_ILnt_config;
            case 'SANT'
                c = VWM_SAnt_config;
            case 'MIXNT'
                c = VWM_MIXnt_config(nSetSize); % note MIX model requires number of set size levels
            case 'EPNT'
                c = VWM_EPnt_config;
            case 'VPNT'
                c = VWM_VPnt_config;
            case 'VPCAPNT'
                c = VWM_VPcapnt_config;
            case 'COSSANT'
                c = VWM_cosSAnt_config;     
        end        
        %%
        % do it, fit the model
        c = c.fitFun(data, c);
        results.fitResults{iModel} = {models{iModel}, c};
    else
        error('We have no such model!');
    end
end

% %% We do some simple analysis
% % find the best-fitting model
% AIC = zeros(1, results.nModeltoFit);
% AICc = zeros(1, results.nModeltoFit);
% BIC = zeros(1, results.nModeltoFit);
% for iModel = 1:results.nModeltoFit
%     tmp = results{};
%     AIC(iModel) = tmp.best(end-2); % We extract the AIC;
%     AICc(iModel) = tmp.best(end-1); % We extract the AICc;
%     BIC(iModel) = tmp.best(end); % We extract the BIC; 
% end
% [~,idx] = min(AIC);
% results.bestModel_AIC = results.models{idx};
% [~,idx] = min(AICc);
% results.bestModel_AICc = results.models{idx};
% [~,idx] = min(BIC);
% results.bestModel_BIC = results.models{idx};
% 
% % Posterior probability
% results.PP_AIC=IC2PP(AIC);
% results.PP_BIC=IC2PP(BIC);
% results.PP_AICc=IC2PP(AICc);
% 
% results.AIC = AIC;
% results.AICc=AICc;
% results.BIC = BIC;

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

