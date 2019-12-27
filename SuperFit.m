%% Super Fit
% Ma, Tianye
% 11.24.2019

%% Prologue
clc
clear variables
addpath('Your Datapath')
load('schidata_ori.mat') % 9 schizophrenia subjects (orientation recall)
data_schi=data;
load('normdata2.mat') % 62 normal subjects (color recall)
data_norm=data;
Nsubj_schi=length(data_schi);
Nsubj_norm=length(data_norm);
ModelSpace={'IL', 'MIX', 'COSSA', 'SA', 'EP', 'VP', 'VPCAP'};
% ModelSpace={'cosSA'};

%% HaHaHa
Results_Schi=cell(1,Nsubj_schi);
Results_Norm=cell(1,Nsubj_norm);
for i=1:Nsubj_schi
    Data=data_schi{i};
    Data=Data(Data(:,3)==0,:);
    data_subj.probe=Data(:,5);
    data_subj.resp=Data(:,6);
    data_subj.error=Data(:,1);
    data_subj.N=Data(:,2);
    Results_Schi{i}=fitVWMmodels(data_subj,ModelSpace);
end

for i=1:Nsubj_norm
    Data=data_norm{i};
    Data=Data(Data(:,3)==0,:);
    data_subj.probe=Data(:,5);
    data_subj.resp=Data(:,6);
    data_subj.error=Data(:,1);
    data_subj.N=Data(:,2);
    Results_Norm{i}=fitVWMmodels(data_subj,ModelSpace);
end

%% Epilogue
rmpath('Your Datapath')
