function c = VWM_SA_fit(data, c)
% This is the main function to fit item-limit model
% <data>:
% <c>: configuration struct includes optimization/

%% preprocess visual working memory data
% you can do any preprocessing here
data = prepVWMdata(data);
data.gvar.kappa_max      = 700; % this is due to matlab limit
data.gvar.kappa_map      = linspace(0,data.gvar.kappa_max,1e5);
data.gvar.J_map          = data.gvar.kappa_map.*besseli(1,data.gvar.kappa_map,1)./besseli(0,data.gvar.kappa_map,1);

%% ========= use bads to optimization ========
c = optimnll(data, c);

