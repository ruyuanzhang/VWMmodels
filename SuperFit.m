%% Super Fit
% Ma, Tianye
% 11.24.2019

%% Prologue
clc
clear variables
addpath('C:\Users\LabLogin\Desktop\Tianye\WM_Models\Data_GZ')
load('schidata_ori.mat')
Nsubj=length(data);
ModelSpace={'IL', 'MIX', 'COSSA', 'SA', 'EP', 'VP', 'VPCAP'};
% ModelSpace={'cosSA'};

%% HaHaHa
for i=2
    Data=data{i};
    data_subj.probe=Data(:,5);
    data_subj.resp=Data(:,6);
    data_subj.error=Data(:,1);
    data_subj.N=Data(:,2);
    Results{i}=fitVWMmodels(data_subj,ModelSpace);
end

%% Epilogue
rmpath('C:\Users\LabLogin\Desktop\Tianye\WM_Models\Data_GZ')
