function c = VWM_EP_config(nFit)

if notDefined('nFit')
    nFit = 20;
end

% some setting
maxJ1bar = 700;

%% define optimization variable
opt.PLB = [0, 0, 0];
opt.PUB = [maxJ1bar, 10, maxJ1bar];
opt.LB = [0, 0, 0];
opt.UB = [maxJ1bar, 10, maxJ1bar];
opt.paramLabels = {'J1bar', 'power', 'kappa_r'};

% 
nvars = length(opt.PLB);
opt.options = bads('defaults');
opt.options.MaxIter = maxJ1bar*nvars;

x0 = [opt.PLB(1)+eps:(opt.PUB(1)-opt.PLB(1))/(nFit-1):opt.PUB(1)-eps;
    opt.PLB(2)+eps:(opt.PUB(2)-opt.PLB(2))/(nFit-1):opt.PUB(2)-eps;
    opt.PLB(3)+eps:(opt.PUB(3)-opt.PLB(3))/(nFit-1):opt.PUB(3)-eps];
opt.x0=x0';  %x0,  nFit x n params

opt.nFit = nFit;
opt.nvars = nvars;

c.opt=opt;
%% define fitting function
c.fitFun = @VWM_EP_fit;

%% define negative loglikelihood function
c.negLogLikeliFun=@VWM_EP_nll;

end