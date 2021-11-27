%% Super Fit
% Here is an example script to fit multiple subject
%% Prologue
clear all;close all;clc;
clear variables
addpath(genpath('~/Dropbox/stonesync/fitVWMmodels/VWMmodels')); % you should change the path
ModelSpace={'IL', 'MIX', 'COSSA', 'SA', 'EP', 'VP', 'VPCAP'};

%% load the data of all subject
% load data

%% Run it
try
    if isempty(gcp)
        pobj = parpool(9); % open 9 cores for parrallel computing
    end
     
    % Run schizo data
    Results_Schi=cell(1,Nsubj_schi);
    starttime = gettimestr;
    parfor iSubj=1:Nsubj_schi % loop subject
        iSubj
        Results_Schi{iSubj}=fitVWMmodels(data_Schi{iSubj}, ModelSpace);
    end
        
    endtime = datestr(now);
    
    % notify me when it is done
    %mailcontent = sprintf('start: %s \n end: %s', starttime, endtime);
    
    % save data
    filename = ['./' gettimestr 'ori9schi26norm.mat'];
    save(filename);
catch ME
    % notify me when error occurs
    %mailme('');
    rethrow(ME) % show the error message
end 



%% Epilogue
rmpath(genpath('~/Dropbox/stonesync/fitVWMmodels/VWMmodels'))
