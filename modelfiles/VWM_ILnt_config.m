function c = VWM_ILnt_config(varparams)
 
% varparams is a struct that accept input from outside functions

nFit = varparams.nFit;
optimizer = varparams.optimizer;

% some setting
maxJ1bar = 700;

%% define optimization variable
opt.PLB = [0, 0, 0];
opt.PUB = [maxJ1bar, 20, 1];
opt.LB = [0, 0, 0];
opt.UB = [maxJ1bar, 20, 1];
opt.paramLabels = {'kappa_r', 'K', 'pNT'};

nvars = length(opt.PLB);
opt.options = bads('defaults');
opt.options.MaxIter = maxJ1bar*nvars;
% initial guess
opt.x0 = [opt.PLB(1):(opt.PUB(1)-opt.PLB(1))/(nFit-1):opt.PUB(1);
    opt.PLB(2):(opt.PUB(2)-opt.PLB(2))/(nFit-1):opt.PUB(2);
    opt.PLB(3):(opt.PUB(3)-opt.PLB(3))/(nFit-1):opt.PUB(3);
    ];
opt.x0=opt.x0';

opt.nFit = nFit;
opt.nvars = nvars;
opt.optimizer = optimizer;

c.opt=opt;
%% define fitting function
c.fitFun = @VWM_ILnt_fit;

%% define negative loglikelihood function
c.negLogLikeliFun=@VWM_ILnt_nll;


end