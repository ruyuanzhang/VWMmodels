%% fit one example subject
clear all;close all;clc;
clear variables

data = load('sampledata.mat'); 
data = data.data;
c = VWM_cosSAnt_config;
c = c.fitFun(data,c);
