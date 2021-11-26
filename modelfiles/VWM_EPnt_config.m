function c = VWM_EPnt_config(nFit)

if notDefined('nFit')
    nFit = 20;
end

% some setting
maxJ1bar = 700;

%% define optimization variable
opt.PLB = [0, 0, 0, 0];
opt.PUB = [maxJ1bar, 10, maxJ1bar, 1];
opt.LB = [0, 0, 0, 0];
opt.UB = [maxJ1bar, 10, maxJ1bar, 1];
opt.paramLabels = {'J1bar', 'power', 'kappa_r', 'pNT'};

% 
nvars = length(opt.PLB);
opt.options = bads('defaults');
opt.options.MaxIter = maxJ1bar*nvars;

x0 = [opt.PLB(1):(opt.PUB(1)-opt.PLB(1))/(nFit-1):opt.PUB(1);
    opt.PLB(2):(opt.PUB(2)-opt.PLB(2))/(nFit-1):opt.PUB(2);
    opt.PLB(3):(opt.PUB(3)-opt.PLB(3))/(nFit-1):opt.PUB(3);
    opt.PLB(4):(opt.PUB(4)-opt.PLB(4))/(nFit-1):opt.PUB(4);
    ];
opt.x0=x0';  %x0,  nFit x n params

opt.nFit = nFit;
opt.nvars = nvars;

c.opt=opt;
%% define fitting function
c.fitFun = @VWM_EPnt_fit;

%% define negative loglikelihood function
c.negLogLikeliFun=@VWM_EPnt_nll;

end