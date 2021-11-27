# VWMmodels
2021-11-26 by Ru-Yuan Zhang

A set of computational models (mostly probablistic models) of visual working memory



## Installnation

To run this, you need to add below files to your matlab path

* Ru-Yuan Zhang's matlab utility functions, https://github.com/ruyuanzhang/RZutil
* bads optimization toolbox, https://github.com/lacerbi/bads



## Model specifications

Currently we have 7 variants of models

* IL: item-limit model
* MIX: mixture model
* SA: slot-plus-averaging model
* EP: equal presicion model
* VP: variable precision model
* VPcap: variable-precision-plus-capacity model
* cosSA: slots-plus-averaging-with-cosine-precision model

 The total number of model is thus 7 x 2 = 14.



All model related files are located in 'modelfiles' directory. We label the model as 'VWM\_<model><nt>_config' format.  <model> is the model name, <nt> indicates whether to includes non-target swap error. For example,  VWM\_SA\_* indicates the slot-plus-averaging model without non-target swaping error, and VWM\_SAnt\_* indicates vice versa.



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
fitVWMmodels(data, {'SA','VP', 'EP'}) % fit 'SA', 'VP', and 'EP'
fitVWMmodels(data, 'ALL') % fit 'SA', 'VP', and 'EP'
```



# Intepret the result

