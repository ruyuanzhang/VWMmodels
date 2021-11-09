%% fit one example subject
clear all;close all;clc;
clear variables
addpath(genpath('~/Dropbox/stonesync/fitVWMmodels/VWMmodels')); % you should change the path
ModelSpace={'IL', 'MIX', 'COSSA', 'SA', 'EP', 'VP', 'VPCAP'};
data = load('sampledata.mat'); 
fitVWMmodels(data,ModelSpace);