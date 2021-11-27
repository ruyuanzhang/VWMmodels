function c = VWM_MIX_fit(data, c)
% This is the main function to fit item-limit model
% <data>:
% <c>: configuration struct includes optimization/

%% preprocess visual working memory data
% you can do any preprocessing here
data = prepVWMdata(data);

%% ========= use bads to optimization ========
c = optimnll(data, c);
