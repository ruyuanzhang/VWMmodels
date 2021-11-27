# VWMmodels
2021-11-26 by Ru-Yuan Zhang

A set of computational models (mostly probablistic models) of visual working memory



## Installnation

To run this, you need to add repositories below to your matlab path

* Ru-Yuan Zhang's matlab utility functions, https://github.com/ruyuanzhang/RZutil
* bads optimization toolbox, https://github.com/lacerbi/bads



## Model specifications

Currently we have 7 variants of models:

* IL: item-limit model
* MIX: mixture model
* SA: slot-plus-averaging model
* EP: equal presicion model
* VP: variable precision model
* VPcap: variable-precision-plus-capacity model
* cosSA: slots-plus-averaging-with-cosine-precision model

We also consider whether a model has non-target swap errors.  The total number of models is thus 7 (variants above) x 2 (non-target swap errors or not) = 14.



All model-related files are located in the 'modelfiles' directory. We label a model as 'VWM\_<model><nt>_\*\*' format.  <model> is the model name, <nt> indicates whether to include non-target swap errors. For example,  VWM\_SA\_* indicates the slot-plus-averaging model without considering non-target swaping errors, and VWM\_SAnt\_* indicates vice versa.



Each model has three files:

* VWM\_<model><nt>_config. Configurations for model fiting
* VWM\_<model><nt>_fit. The main function to preprocess data and conduct fit multiple times
* VWM\_<model><nt>_nll.  Negative log likelihood function (this is the key function of a model)



## Data preparation

Suppose you have run an experiment including <nTrials> trials, Basically you need to construct a data structure that contains the fields a below:

* N: a 1 x nTrials vector, the set size in each trial
* probe: a 1 x nTrials vector, with range of [1, 180], the color of probe in each trial
* resp: a 1 x nTrials vector, with range of [1, 180], the color of response in each trial
* error: a 1 x nTrials vector, with range of [-90, 89], the circular error between probe and resp with mode at 180
* distrError: a 1 x nTrials cell, each element is a vector of errors with respect to other distractors besides resp



Note that <distrError> is only useful when fitting non-target swapping(i.e., 'nt') models, not necessary for other models

# Fit single subject data

Once you have the data structure, you can fit a single subject data as:

```matlab
# fit a subset of models to a subject's data
result = fitVWMmodels(data, {'SA','VP', 'EP'}) % fit 'SA', 'VP', and 'EP'
result = fitVWMmodels(data, 'ALL') % fit all models
```



# Intepret the result

The variable `result` have following fields:

* setSize: a vector of set size used on this subject, e.g., [1 3 5 8]
* models: a cell of model names
* nFit: int, how many times to fit a model with random initialized parameters
* nModeltoFit: int, number of total models to fot
* fitResults: a cell of fitting results, with each element for each model. 

```matlab
% get the result for a model
c = result.fitResults{1}{2}
% c has fields:
%	 opt: options for model fitting
%	 fitFun: main fitting function to use
%	 negLogLikeliFun: negative log likelihood function to optimize
%  result: structure for fitting result, which will be further analyze as below
c = result.fitResults{1}{2}.result
%  fitResults: nFit x nvars matrix, fitted parameters for all nFit times optimization
%  modelMetrics: nFit x nmodelmetrics, calculated model metrics
%  modelMetricLabels: a cell of names of model metrics.
%	 bestFit: best fitting parameters with the lowest neglh value among all nFit optimizations 
```

