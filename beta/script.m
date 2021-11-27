%% Prologue
clear all;close all;clc;
clear variables
addpath(genpath('~/Dropbox/stonesync/fitVWMmodels/VWMmodels'));
ModelSpace={'IL', 'MIX', 'COSSA', 'SA', 'EP', 'VP', 'VPCAP'};


%% fit 9 schizo and 26 norm subjects on orientation recall
data = load('../data.mat'); % 9 schizophrenia subjects (orientation recall)

data_Schi = data.schizoData;
data_Norm = data.controlData;
Nsubj_schi=length(data_Schi);
Nsubj_norm=length(data_Norm);

%% Run it
try
    if isempty(gcp)
        pobj = parpool(40);
    end
    
    % Run schizo data
    Results_Schi=cell(1,Nsubj_schi);
    
    starttime = gettimestr;
    parfor iSubj=1:Nsubj_schi % loop parallel computing
        iSubj
        Results_Schi{iSubj}=fitVWMmodels(data_Schi{iSubj}, ModelSpace);
    end
    
    endtime = gettimestr;
    
    % notify me when it is done
    % mailcontent = sprintf('start: %s \n end: %s', starttime, endtime);
    % mailme('ori9schi26normdone',mailcontent,'miao0607','ruyuanzhang@gmail.com');
    
    % save data
    filename = ['../' gettimestr '78child.mat'];
    save(filename);
catch ME
    % notify me when error occurs
    %mailme('ori9schi26normERROR','something wrong','miao0607','ruyuanzhang@gmail.com');
    rethrow(ME) % show the error message
end 



%% Epilogue
rmpath(genpath('~/Dropbox/stonesync/fitVWMmodels/VWMmodels'))
