function c = VWM_VP_fit(data, c)
% This is the main function to fit EP model
% <data>:
%   N: a vector containing set size levels
%   error: a vector of response error (probe-resp)
%   distrError: a cell, response errors with respect to distractors 
% 
% <c>: a struct, contains configuration variables, e.g. optimization

%% preprocess visual working memory data
% you can do any preprocessing here. This could be model specific
data = prepVWMdata(data);
% add more helping variable to data
data.gvar.kappa_max      = 700; % this is due to matlab limit
data.gvar.kappa_map      = linspace(0, data.gvar.kappa_max,1e5);
data.gvar.J_map          = data.gvar.kappa_map.*besseli(1, data.gvar.kappa_map,1)./besseli(0, data.gvar.kappa_map,1);
data.gvar.nMCSamples     = 10000;                 % number of MC samples to draw when computing model predictions (Paper: 1000)

%% ========= use bads to optimization ========
c = optimnll(data, c);
