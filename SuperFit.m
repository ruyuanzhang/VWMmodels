%% Super Fit
% Ma, Tianye
% 11.24.2019


%% Prologue
clear all;close all;clc;
clear variables
addpath(genpath('~/Dropbox/stonesync/fitVWMmodels/VWMmodels'));
ModelSpace={'IL', 'MIX', 'COSSA', 'SA', 'EP', 'VP', 'VPCAP'};


% %% Fit 62 normal subjects on color recall
% data = load('normdata2.mat'); % 62 normal subjects (color recall)
% data_norm=data.data;
% Nsubj_norm=length(data_norm);
% 
% Results_Norm=cell(1, Nsubj_norm);
% 
% % prepare the data
% data_subj = cell(1,Nsubj_norm);
% for iSubj = 1:Nsubj_norm
%     Data=data_norm{iSubj};
%     Data=Data(Data(:,3)==0,:);
%     data_subj{iSubj}.probe=Data(:,5);
%     data_subj{iSubj}.resp=Data(:,6);
%     data_subj{iSubj}.error=Data(:,1);
%     data_subj{iSubj}.N=Data(:,2);
% end
% 
% % Run it
% try
%     if isempty(gcp)
%         pobj = parpool(16);
%     end
%     
%     starttime = gettimestr;
%     parfor iSubj=1:Nsubj_norm
%         iSubj
%         Results_Norm{iSubj}=fitVWMmodels(data_subj{iSubj},ModelSpace);
%     end
%     endtime = gettimestr;
%     
%     % notify me when it is done
%     mailcontent = sprintf('start: %s \n end: %s', starttime, endtime);
%     mailme('62normalfitdone','mailcontent','miao0607','ruyuanzhang@gmail.com');
%     
%     % save data
%     filename = ['../WMschizo/' gettimestr 'results62normori.mat'];
%     save(filename);
% catch ME
%     % notify me when error occurs
%     mailme('62normalfitERROR','something wrong','miao0607','ruyuanzhang@gmail.com');
%     rethrow(ME) % show the error message
% end 


%% test one subject
%test = fitVWMmodels(data_subj{10}, ModelSpace);


%% fit 9 schizo and 26 norm subjects on orientation recall
data = load('../WMschizo/oritaskschizo9normal26.mat'); % 9 schizophrenia subjects (orientation recall)
data_Schi = data.schizoData;
data_Norm = data.controlData;
Nsubj_schi=length(data_Schi);
Nsubj_norm=length(data_Norm);

%% Run it
try
    if isempty(gcp)
        pobj = parpool(9);
    end
    
    starttime = datestr(now);
    
    % Run schizo data
    Results_Schi=cell(1,Nsubj_schi);
    starttime = gettimestr;
    parfor iSubj=1:Nsubj_schi
        iSubj
        Results_Schi{iSubj}=fitVWMmodels(data_Schi{iSubj}, ModelSpace);
    end
    
    pobj=gcp;
    if isempty(gcp)
        pobj = parpool(13);
    end
    % Run normal data
    Results_norm=cell(1,Nsubj_norm);
    parfor iSubj=1:Nsubj_norm
        iSubj
        Results_Norm{iSubj}=fitVWMmodels(data_Norm{iSubj}, ModelSpace);
    end
    
    endtime = datestr(now);
    
    % notify me when it is done
    mailcontent = sprintf('start: %s \n end: %s', starttime, endtime);
    mailme('ori9schi26normdone',mailcontent,'miao0607','ruyuanzhang@gmail.com');
    
    % save data
    filename = ['../WMschizo/' gettimestr 'ori9schi26norm.mat'];
    save(filename);
catch ME
    % notify me when error occurs
    mailme('ori9schi26normERROR','something wrong','miao0607','ruyuanzhang@gmail.com');
    rethrow(ME) % show the error message
end 



%% Epilogue
rmpath(genpath('~/Dropbox/stonesync/fitVWMmodels/VWMmodels'))
